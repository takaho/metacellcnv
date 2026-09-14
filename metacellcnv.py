#!/usr/bin/env python3
"""
metacellcnv.py -- CNV-based malignant metacell identification for tumor scRNA-seq

Two-branch pipeline: builds metacells from single-cell RNA-seq data (native
implementation, no SEACells/infercnvpy dependency), estimates copy-number
variation (CNV) per metacell, calls malignant vs. normal metacells from the
CNV signal, and optionally runs differential expression (pyDESeq2) between
malignant and normal groups.

Input: one CellRanger `filtered_feature_bc_matrix/` directory (MTX format)
per sample, plus a gene-annotation GTF used to map genes to genomic
coordinates for the CNV step.

Output: QC'd/preprocessed AnnData, metacell assignments, CNV estimates,
malignant/normal calls, optional DE results, and supporting QC and
interpretation reports (see README.md for the full file list).

Usage:
    conda create -n scrna-cnv python=3.10 -y && conda activate scrna-cnv
    pip install scanpy pydeseq2 decoupler scrublet
    python metacellcnv.py --help

    # Reference implementations (SEACells / infercnvpy) are optional and
    # only needed to cross-check the native metacell/CNV implementations:
    #   pip install SEACells infercnvpy
    # then run verify_native_metacells.py / verify_native_cnv.py.

See README.md for design rationale, validation results, and the
full version history.

References:
- pyDESeq2:   https://pydeseq2.readthedocs.io/
- Scrublet (scanpy wrapper): https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html
- MAD-based adaptive QC: https://bioconductor.org/books/release/OSCA.basic/quality-control.html
- scDblFinder (design inspiration): https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
- Numbat (design inspiration for CNV-consistency-based clone/doublet discrimination):
  https://www.biorxiv.org/content/10.1101/2022.02.07.479314

Design invariants (do not violate without updating README/paper):
1. High-count QC outliers are not excluded outright (not included in qc_fail).
2. Doublet detection runs per cluster (batch_key=coarse_cluster).
3. Doublet-rescue CNV-consistency checks run at single-cell level.
4. Final CNV estimation and DE analysis run at metacell level, with
   sample_id included as a covariate.
5. coarse_cluster is kept independent of the malignant/normal call itself.
6. infercnvpy is experimental; cross-check important conclusions with other
   methods.
7. GTF attribute keys and chromosome naming conventions must always be
   verified against the real data.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# ---------------------------------------------------------------------------
# Shared building blocks live in scrna_common.py (split out so the
# preprocessing steps can be reused standalone, outside the CNV analysis).
# Behavior is unchanged; re-exported here under the same names so existing
# code importing this module keeps working.
# ---------------------------------------------------------------------------
from scrna_common import (  # noqa: F401  (re-export)
    BUILTIN_MARKER_SPECIES,
    MarkerTable,
    load_marker_table,
    resolve_marker_path,
    CANDIDATE_GENE_ID_KEYS,
    CELLRANGER_DIRS,
    CELLS_PER_METACELL,
    CHROMOSOME_MAP_PATH,
    GTF_GENE_ID_ATTR,
    GTF_PATH,
    MAX_CELLS_WARN,
    MITO_CHROMOSOMES,
    MITO_MAX_LENGTH,
    MITO_MAX_NEAR_ZERO,
    MITO_MIN_GENES,
    MITO_MIN_SPEARMAN,
    MITO_ND6_RANK_ALERT,
    MITO_NEAR_ZERO_PCT,
    MITO_PC_CORR_ALERT,
    MITO_PREFIXES,
    MITO_REFERENCE_PROFILE,
    MITO_SYMBOLS,
    NUCLEAR_MITO_GENES,
    NUCLEAR_MITO_PREFIXES,
    N_HVG,
    N_PCS,
    OUTPUT_DIR,
    RIBO_PREFIXES,
    _infer_sample_id,
    _leiden_supports_igraph,
    _open_maybe_gzip,
    _to_dense,
    _warn_if_stale_cache,
    _write_h5ad,
    check_embedding_mito_influence,
    check_exclude_chromosomes,
    check_gene_overlap,
    check_mito_composition,
    detect_doublets_cluster_aware,
    detect_mito_chromosomes,
    ensure_sample_id,
    flag_mito_genes,
    inspect_gtf,
    load_and_preprocess,
    load_chromosome_map,
    load_mito_reference_profile,
    log,
    looks_like_raw_counts,
    mad_outlier_mask,
    normalize_chromosome_names,
    parametric_qc_filter,
    parse_gtf_gene_positions,
    validate_dimensionality,
    warn,
)


# sample_id for each CellRanger output. If None, it is auto-generated from
# the directory name; obs['sample_id'] is always set explicitly, even for a
# single sample.
SAMPLE_IDS: list[str] | None = None

# Gene coordinate annotation (required for CNV estimation).
# Example: CellRanger reference's refdata-gex-.../genes/genes.gtf
__version__ = "3.3"



# Chromosome names to exclude from CNV estimation.
# GENCODE/UCSC style: ("chrX", "chrY").
# For NCBI/RefSeq accession-style names (e.g. "NC_051844.1"), chrX/chrY do
# not exist and specifying them here has no effect, so the default is an
# empty tuple. Set this once the X/Y accessions have been identified via the
# assembly report.
EXCLUDE_CHROMOSOMES: tuple[str, ...] = ()



# Marker genes used for cell-type annotation.
# Coverage note: with only 3 sets (T/NK, B, Myeloid) the major tumor
# microenvironment components have no markers to score against, so those
# clusters all fall into "Other" -- not because they can't be classified,
# but because they were never asked about. Sets deliberately do not share
# genes with each other (shared genes would shrink the 1st-vs-2nd margin
# used for classification below).
MARKER_SETS: dict[str, list[str]] = {
    "T/NK": ["CD3D", "CD3E", "CD3G", "CD2", "LCK", "ITK", "THEMIS", "SKAP1", "CD247",
             "IL7R", "CD28", "GZMA", "GZMB", "NKG7", "GNLY", "KLRD1", "PRF1", "EOMES"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19", "PAX5", "BANK1", "EBF1", "BLK", "CD22", "FCRL1"],
    "Plasma": ["JCHAIN", "MZB1", "XBP1", "DERL3", "PRDM1", "TNFRSF17", "SDC1"],
    "Myeloid": ["CD68", "CSF1R", "MRC1", "CD163", "C1QA", "C1QB", "C1QC", "MSR1",
                "AIF1", "TYROBP", "LYZ", "CD14", "ITGAM", "FCGR3A", "MARCO", "VSIG4",
                "FCER1G", "MNDA", "IRF8", "ZBTB46", "BATF3", "FLT3", "S100A8", "S100A9",
                "MMP9", "CSF3R"],
    "Mast": ["KIT", "CPA3", "MS4A2", "TPSAB1", "CMA1", "GATA2", "HDC"],
    "Endothelial": ["PECAM1", "CDH5", "VWF", "KDR", "CLDN5", "EGFL7", "TEK", "ERG",
                    "FLT1", "ESAM"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FBN1", "POSTN",
                   "THY1", "FAP", "PDGFRB", "ACTA2", "RGS5", "TAGLN", "NOTCH3", "MYH11"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "SFN", "KRT5", "KRT14"],
}

# Fallback used if every set above scores zero: a conservative pan-leukocyte
# panel. This sacrifices cell-type resolution, but the goal here is only to
# secure a "normal reference" for CNV estimation.
FALLBACK_MARKER_SETS: dict[str, list[str]] = {
    "Leukocyte": ["PTPRC", "LAPTM5", "CORO1A", "CD52", "SRGN"],
}

# Known normal cell types usable as the malignant/normal CNV reference
# (reference_cat). Fibroblast / Epithelial are kept in MARKER_SETS but are
# NOT used as a reference by default: in carcinomas the tumor itself can be
# epithelial, and in fibrous/mesenchymal tumors it can look fibroblast-like,
# so including them would cancel out the CNV signal. Empirically, "restrict
# to safe types, then take the top-scoring set" backfires (on real data 4 of
# 5 malignant clusters were misassigned as Plasma); the correct order is
# "take the top-scoring set from all sets, and only use it as a reference if
# that top set is a safe type." If the tissue-of-origin type is known and
# should be included deliberately, add it via --normal-celltype.
KNOWN_NORMAL_CELLTYPES = ["T/NK", "B", "Plasma", "Myeloid", "Mast", "Endothelial", "Leukocyte"]

#: Types that can themselves be the tumor, so excluded from the normal
#: reference by default.
MALIGNANCY_CAPABLE_CELLTYPES = ("Fibroblast", "Epithelial")


def apply_marker_table(spec: str | None) -> "object":
    """Load the marker table given via --markers and replace this module's
    default marker sets with it.

    `MARKER_SETS` / `FALLBACK_MARKER_SETS` / `KNOWN_NORMAL_CELLTYPES` are
    names referenced elsewhere in this module, so their contents are updated
    in place rather than rebound.
    """
    tbl = load_marker_table(spec)
    MARKER_SETS.clear(); MARKER_SETS.update(tbl.marker_sets)
    FALLBACK_MARKER_SETS.clear(); FALLBACK_MARKER_SETS.update(tbl.fallback_sets)
    KNOWN_NORMAL_CELLTYPES[:] = [c for c in tbl.normal_celltypes
                                 if c not in tbl.fallback_sets] + list(tbl.fallback_sets)
    globals()["MARKER_TABLE"] = tbl
    return tbl


def set_normal_celltypes(extra: list[str] | None) -> None:
    """Add extra cell types to the normal reference set, via --normal-celltype."""
    if not extra:
        return
    for name in extra:
        name = str(name).strip()
        if not name or name in KNOWN_NORMAL_CELLTYPES:
            continue
        KNOWN_NORMAL_CELLTYPES.append(name)
        if name in MALIGNANCY_CAPABLE_CELLTYPES:
            warn(
                f"Added '{name}' to the normal reference. This type can itself be "
                "the tumor, so including malignant clusters in the reference will "
                "cancel out the CNV signal. Use this only when the tissue type is known."
            )
        else:
            log(f"Added '{name}' to the normal reference")

# ---------------------------------------------------------------------------
# Step -1. Validate and record the runtime environment
# ---------------------------------------------------------------------------

def check_runtime_environment(strict: bool = True) -> None:
    """Validate the Python version, interpreter, ssl, and required packages.

    Fail-fast check for an old/wrong Python version, a missing ssl module,
    or missing required packages.
    """
    log("=== Runtime environment check ===")
    problems: list[str] = []

    v = sys.version_info
    log(f"Python {v.major}.{v.minor}.{v.micro} ({sys.executable})")
    if (v.major, v.minor) < (3, 10):
        problems.append(
            f"Python {v.major}.{v.minor} is older than 3.10. Older scanpy versions "
            "(e.g. 1.7.2) will be resolved instead, causing incompatibilities "
            "throughout the pipeline. Create and activate a dedicated conda "
            "environment with python=3.10+ (e.g. scrna-cnv) before running this."
        )
    if "seurat" in sys.executable.lower():
        warn(
            "The running Python's environment name looks like it is meant for "
            "Seurat/R. Please check that the intended dedicated environment is "
            "activated."
        )

    try:
        import ssl  # noqa: F401

        log(f"ssl available: {ssl.OPENSSL_VERSION}")
    except ImportError as exc:  # pragma: no cover - environment dependent
        problems.append(f"ssl module is not available ({exc}). Using a conda environment is recommended.")

    for pkg in ("numpy", "pandas", "scipy", "anndata", "scanpy", "sklearn"):
        _report_version(pkg, problems, required=True)
    for pkg in ("pydeseq2", "scrublet"):
        _report_version(pkg, problems, required=True)
    for pkg in ("decoupler",):
        _report_version(pkg, problems, required=False)

    if problems:
        for p in problems:
            print(f"[pipeline][FAIL] {p}", file=sys.stderr)
        if strict:
            raise SystemExit(
                "Environment check failed. Resolve the issues above and re-run "
                "(or pass --skip-env-check to force continuation)."
            )


def _report_version(pkg: str, problems: list[str], required: bool) -> None:
    import importlib

    try:
        mod = importlib.import_module(pkg)
        log(f"{pkg} {getattr(mod, '__version__', 'unknown')}")
    except ImportError as exc:
        msg = f"Cannot import {pkg}: {exc}"
        if required:
            problems.append(msg)
        else:
            warn(f"(optional) {msg}")


def write_environment_lock(out_dir: Path) -> Path | None:
    """Save `pip list` output to environment.lock.txt."""
    out_dir.mkdir(parents=True, exist_ok=True)
    lock_path = out_dir / "environment.lock.txt"
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "pip", "list", "--format=freeze"],
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
        )
    except Exception as exc:  # pragma: no cover - environment dependent
        warn(f"Failed to create environment.lock.txt: {exc}")
        return None

    if proc.returncode != 0:
        warn(f"pip list failed (returncode={proc.returncode}): {proc.stderr.strip()[:200]}")
        return None

    header = (
        f"# python: {sys.version.replace(chr(10), ' ')}\n"
        f"# executable: {sys.executable}\n"
    )
    lock_path.write_text(header + proc.stdout, encoding="utf-8")
    log(f"Saved environment.lock.txt: {lock_path}")
    return lock_path


def check_marker_availability(
    var_names, marker_sets: dict[str, list[str]], fallback: dict[str, list[str]] | None = None
) -> None:
    """Report which marker genes are present in var_names for each set."""
    log("=== Marker gene check ===")
    var_set = set(map(str, var_names))
    any_hit = False
    for label, genes in marker_sets.items():
        present = [g for g in genes if g in var_set]
        missing = [g for g in genes if g not in var_set]
        if present:
            any_hit = True
            log(f"'{label}': {len(present)}/{len(genes)} matched {present} (missing: {missing})")
        else:
            warn(f"'{label}': no matching genes (candidates: {genes})")
    if not any_hit and fallback:
        for label, genes in fallback.items():
            present = [g for g in genes if g in var_set]
            if present:
                any_hit = True
                log(f"(fallback) '{label}': {present} available")
            else:
                warn(f"(fallback) '{label}': no matching genes (candidates: {genes})")
    if not any_hit:
        raise ValueError(
            "No matching genes in any marker set (including fallback). "
            "annotate_coarse_celltype() will fail. Adjust MARKER_SETS to match "
            "this species' gene symbols or LOC IDs."
        )


def preflight(
    cellranger_dirs: list[str],
    gtf_path: str,
    gtf_gene_id_attr: str | None,
    chromosome_map_path: str | None = None,
    mito_chromosomes: tuple[str, ...] | list[str] | None = None,
) -> dict:
    """Run all cheap validation checks before the heavy pipeline steps.

    Verifies that CellRanger outputs and required files exist, checks the
    GTF format/attribute keys/chromosome naming, computes the overlap
    between features.tsv gene names and the GTF, and checks marker gene
    availability. The goal is to fail before the expensive mtx load.
    """
    log("=== CellRanger output check ===")
    if not cellranger_dirs:
        raise ValueError("No CellRanger output directories were given")
    if len(cellranger_dirs) == 1:
        log("Single sample: obs['sample_id'] will be set explicitly by this script")

    feature_names: set[str] = set()
    for d in cellranger_dirs:
        path = Path(d)
        if not path.exists():
            raise FileNotFoundError(f"Directory does not exist: {d}")
        matrix = next((p for p in (path / "matrix.mtx.gz", path / "matrix.mtx") if p.exists()),
                      path / "matrix.mtx.gz")
        barcodes = next((p for p in (path / "barcodes.tsv.gz", path / "barcodes.tsv") if p.exists()),
                        path / "barcodes.tsv.gz")
        features = next((p for p in (path / "features.tsv.gz", path / "features.tsv",
                                     path / "genes.tsv.gz", path / "genes.tsv") if p.exists()),
                        path / "features.tsv.gz")
        missing = [p.name for p in (matrix, barcodes, features) if not p.exists()]
        if missing:
            raise FileNotFoundError(f"{d}: required files not found (gzipped or plain is fine): {missing}")

        genes: list[str] = []
        with _open_maybe_gzip(features) as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    genes.append(parts[1])
                elif parts:
                    genes.append(parts[0])
        with _open_maybe_gzip(barcodes) as handle:
            n_barcodes = sum(1 for _ in handle)
        log(f"{d}: {len(genes)} genes / {n_barcodes} barcodes")
        if n_barcodes > MAX_CELLS_WARN:
            warn(f"{d}: {n_barcodes:,} cells is large; SEACells-style optimization may be unstable")
        feature_names |= set(genes)

    gtf_info = inspect_gtf(
        gtf_path, gtf_gene_id_attr, extra_mito_chromosomes=mito_chromosomes
    )
    log("=== var_names vs. GTF overlap ===")
    check_gene_overlap(feature_names, gtf_info["gene_names"])

    log("=== Mitochondrial genes ===")
    mito_in_data = gtf_info["mito_genes"] & feature_names
    if gtf_info["mito_chromosomes"]:
        log(
            f"Of the {gtf_info['mito_chromosomes']} mtDNA sequence(s)' "
            f"{len(gtf_info['mito_genes'])} genes, {len(mito_in_data)} are present in "
            "CellRanger features"
        )
        if mito_in_data:
            log(f"  e.g.: {sorted(mito_in_data)[:13]}")
        else:
            warn(
                "No mtDNA genes are present in features.tsv. The CellRanger "
                "reference may not include mtDNA. pctMT-based QC will not work."
            )
    gtf_info["mito_genes_in_data"] = mito_in_data
    check_marker_availability(feature_names, MARKER_SETS, FALLBACK_MARKER_SETS)

    # Validate exclude_chromosomes against chromosome names normalized the
    # same way as in run_cnv_branch.
    chromosome_map = load_chromosome_map(chromosome_map_path)
    normalized = normalize_chromosome_names(
        pd.DataFrame({"chromosome": sorted(gtf_info["chromosomes"])}),
        chromosome_map,
        mito_chromosomes=gtf_info["mito_chromosomes"],
    )["chromosome"]
    gtf_info["chromosome_map"] = chromosome_map
    gtf_info["normalized_chromosomes"] = set(normalized)
    gtf_info["exclude_chromosomes"] = check_exclude_chromosomes(
        EXCLUDE_CHROMOSOMES, gtf_info["normalized_chromosomes"]
    )
    return gtf_info


def is_known_normal(labels: pd.Series | np.ndarray) -> np.ndarray:
    """Check whether a label is a known normal cell type (also matches
    lineage-suffixed labels like 'Myeloid_3')."""
    arr = np.asarray(pd.Series(labels).astype(str))
    known = tuple(KNOWN_NORMAL_CELLTYPES)
    return np.array(
        [v in known or any(v.startswith(k + "_") for k in known) for v in arr], dtype=bool
    )


def annotate_coarse_celltype(
    adata: ad.AnnData,
    cluster_key: str = "coarse_cluster",
    marker_sets: dict[str, list[str]] | None = None,
    fallback_marker_sets: dict[str, list[str]] | None = None,
    min_score: float | None = None,
    min_zscore: float = 1.0,
    normal_clusters: list[str] | None = None,
    manual_label: str = "Leukocyte",
    per_cluster_labels: bool = True,
    rule: str = "margin",
    min_margin: float = 0.02,
    min_within_z: float = 2.0,
) -> pd.Series:
    """Assign a coarse cell-type label ("T/NK", "B", "Myeloid", "Other", ...)
    to each coarse cluster based on canonical marker-gene scores.

    This is not a substitute for rigorous cell-type annotation; it exists
    only to secure a "normal reference" for CNV estimation. Scoring uses
    only genes actually present in var_names, and cell types with zero
    matching genes are dropped. If every set scores zero, falls back to a
    conservative pan-leukocyte panel (fallback_marker_sets); if that also
    fails, raises ValueError (fail-fast).

    Classification rule (default changed to "margin"):

    Earlier versions used a **relative** threshold: z-score marker scores
    across clusters and keep a cluster if z >= min_zscore. This has a
    structural flaw -- with k clusters, only about the top 1-2 clusters can
    ever satisfy z >= 1 regardless of how uniformly strong the true signal
    is, so most correctly-positive clusters get labeled "Other".

    The current default, "margin", judges **each cluster independently**
    instead of competing clusters against each other:

      margin   = (top score) - (2nd-highest score)   (within the same cluster)
      within_z = (top score - mean of the rest) / (SD of the rest)  (within
                 the same cluster)

    A cluster is kept if: top score > 0 AND margin >= min_margin AND
    within_z >= min_within_z. No cluster's result depends on any other
    cluster's scores.

    `rule` selects the mode:
      "margin"     (default) as above
      "relative-z"           legacy behavior using min_zscore (kept for
                              reproducibility)
      "absolute"             absolute threshold via min_score (min_score
                              required)
    If min_score is given explicitly, "absolute" is used regardless of `rule`.

    Because a relative rule always makes *some* cluster the top scorer, a
    sample with few immune cells risks mislabeling a malignant cluster as
    the normal reference. This function therefore always prints a
    cluster x marker-set score table. To decide manually from that table,
    pass cluster IDs via normal_clusters (bypasses scoring entirely; only
    the given clusters get manual_label).
    """
    if cluster_key not in adata.obs:
        raise KeyError(
            f"obs['{cluster_key}'] does not exist. Run coarse clustering "
            "(Leiden) before annotate_coarse_celltype()."
        )
    if marker_sets is None:
        marker_sets = MARKER_SETS
    if fallback_marker_sets is None:
        fallback_marker_sets = FALLBACK_MARKER_SETS

    clusters = adata.obs[cluster_key].astype(str)

    # --- Manual-override mode (escape hatch when scoring isn't trusted) ---
    if normal_clusters:
        wanted = {str(c) for c in normal_clusters}
        unknown = wanted - set(clusters.unique())
        if unknown:
            raise ValueError(
                f"normal_clusters contains unknown cluster IDs: {sorted(unknown)}"
                f" / actual clusters: {sorted(clusters.unique())}"
            )
        result = pd.Series(
            [f"{manual_label}_{c}" if c in wanted and per_cluster_labels
             else (manual_label if c in wanted else "Other") for c in clusters.values],
            index=adata.obs_names,
            dtype=object,
        ).astype(str)
        log(
            f"Manually assigned clusters {sorted(wanted)} to '{manual_label}*'"
            f" (normal reference): {int(is_known_normal(result).sum())} cells"
        )
        return result

    scores = _score_marker_sets(adata, marker_sets)
    if scores.empty and fallback_marker_sets:
        warn(
            "No primary marker set matched var_names. Falling back to the "
            f"conservative marker set {list(fallback_marker_sets)}."
        )
        scores = _score_marker_sets(adata, fallback_marker_sets)

    if scores.empty:
        raise ValueError(
            "No marker genes from any set were found in var_names. Adjust "
            "marker_sets to match this species' gene symbols (or LOC IDs)."
        )

    # --- Cluster x marker-set score table (always printed for review) ---
    cluster_scores = scores.groupby(clusters, observed=True).mean()
    cluster_sizes = clusters.value_counts().reindex(cluster_scores.index)
    table = cluster_scores.copy()
    table.insert(0, "n_cells", cluster_sizes)
    log("Marker scores (per-cluster mean):")
    print(table.round(4).to_string())

    best_label = cluster_scores.idxmax(axis=1)
    best_score = cluster_scores.max(axis=1)

    # --- Within-cluster comparison (margin / within_z): independent of other clusters ---
    if cluster_scores.shape[1] >= 2:
        second = cluster_scores.apply(lambda r: r.nlargest(2).iloc[1], axis=1)
    else:
        second = pd.Series(0.0, index=cluster_scores.index)
    margin = best_score - second
    rest_mean, rest_sd = [], []
    for idx, lab in best_label.items():
        rest = cluster_scores.loc[idx].drop(labels=[lab])
        rest_mean.append(rest.mean() if len(rest) else 0.0)
        rest_sd.append(rest.std(ddof=0) if len(rest) > 1 else np.nan)
    rest_mean = pd.Series(rest_mean, index=cluster_scores.index)
    rest_sd = pd.Series(rest_sd, index=cluster_scores.index)
    within_z = ((best_score - rest_mean) / rest_sd.replace(0, np.nan)).fillna(0.0)

    # --- Across-cluster comparison (legacy relative z; always shown) ---
    sd = cluster_scores.std(ddof=0).replace(0, np.nan)
    zscores = ((cluster_scores - cluster_scores.mean()) / sd).fillna(0.0)
    best_z = pd.Series(
        [zscores.loc[idx, lab] for idx, lab in best_label.items()], index=cluster_scores.index
    )

    effective_rule = "absolute" if min_score is not None else rule
    if effective_rule == "absolute":
        if min_score is None:
            raise ValueError("rule='absolute' requires min_score to be set.")
        keep = best_score >= min_score
        log(f"Classification mode: absolute threshold (min_score={min_score})")
    elif effective_rule == "relative-z":
        keep = (best_score > 0) & (best_z >= min_zscore)
        log(f"Classification mode: across-cluster relative (score>0 and cross-cluster z>={min_zscore})")
        warn(
            "rule='relative-z' is a legacy compatibility mode. A relative "
            "threshold only passes about the top fifth of k clusters, so most "
            "clusters fall to Other in samples where all clusters share the "
            "same true cell type (the default 'margin' is recommended)."
        )
    elif effective_rule == "margin":
        keep = (best_score > 0) & (margin >= min_margin) & (within_z >= min_within_z)
        log(
            f"Classification mode: within-cluster comparison (score>0, margin>={min_margin}"
            f", within_z>={min_within_z})"
        )
    else:
        raise ValueError(
            f"Invalid rule='{rule}'. Must be one of margin / relative-z / absolute."
        )

    cluster_label = best_label.where(keep, "Other")
    summary = pd.DataFrame(
        {
            "label": cluster_label,
            "best_score": best_score.round(4),
            "margin": margin.round(4),
            "within_z": within_z.round(2),
            "cluster_z": best_z.round(2),
            "reference": np.where(
                is_known_normal(cluster_label), "normal-ref", "-"
            ),
        }
    )
    print(summary.to_string())
    # Flag clusters whose top label can't be used as a reference, instead of
    # silently dropping them.
    blocked = cluster_label[
        cluster_label.isin(MALIGNANCY_CAPABLE_CELLTYPES)
        & ~pd.Series(is_known_normal(cluster_label), index=cluster_label.index)
    ]
    if len(blocked):
        log(
            f"Clusters labeled but not used as normal reference: {dict(blocked)}"
            " (type can itself be the tumor; use --normal-celltype to include it deliberately)"
        )
    if len(set(best_label)) == 1 and best_label.iloc[0] in MALIGNANCY_CAPABLE_CELLTYPES:
        warn(
            f"All {len(best_label)} clusters' top label is '{best_label.iloc[0]}'. "
            "The tumor itself likely belongs to this lineage, and markers alone "
            "cannot distinguish tumor from normal stroma. Specify clusters "
            "explicitly with --normal-clusters, or use --cnv-refine to use "
            "CNV-flat clones as the reference."
        )

    result = clusters.map(cluster_label).astype(str)
    if per_cluster_labels:
        # Distinguish coarse clusters that were assigned the same cell type.
        # Rationale: infercnvpy switches to a "bounded" mode when there are
        # 2+ reference categories (_infercnv_chunk: only counts logFC outside
        # the reference's min/max), which suppresses false positives from
        # cell-type-specific expression bias. A single category collapses to
        # a plain difference and this mechanism does not engage.
        result = pd.Series(
            [f"{lab}_{cl}" if lab in KNOWN_NORMAL_CELLTYPES else lab
             for lab, cl in zip(result.values, clusters.values)],
            index=result.index, dtype=object,
        ).astype(str)
        n_cat = result[result.str.startswith(tuple(KNOWN_NORMAL_CELLTYPES))].nunique()
        log(f"Split normal reference by lineage: {n_cat} categories")
        if n_cat >= 2:
            log("  -> bounded mode is active, suppressing cell-type bias")
        else:
            log(
                "  -> only one category, so bounded mode will not engage. "
                "Consider --cnv-refine to supply multiple CNV-flat clones as references."
            )
    n_normal = int(is_known_normal(result).sum())
    log(f"Cells assigned to normal-reference candidates: {n_normal} / {len(result)}")
    if n_normal == 0:
        warn(
            "Zero cells were assigned a known normal cell type. Check the score "
            "table above to determine which applies:\n"
            "  (a) Labels were assigned, but only Fibroblast / Epithelial ->"
            " the tumor is likely of that lineage. Specify stromal clusters with"
            " --normal-clusters, or add the type with --normal-celltype if the"
            " tissue type is known.\n"
            "  (b) No cluster meets the margin / within_z thresholds ->"
            " lower --marker-min-margin / --marker-min-within-z.\n"
            "  (c) Marker scores are near zero for all sets ->"
            " gene names likely don't match this species; check GTF gene names.\n"
            "  (d) None of the above -> run with --cnv-refine (uses CNV-flat"
            " clones as reference) or --cnv-reference none (uses the mean over"
            " all cells as reference)."
        )
    elif n_normal < 0.01 * len(result):
        warn(
            f"Normal reference is only {n_normal / len(result):.2%} of cells; the"
            " CNV baseline may be unstable. Check the score table for validity."
        )
    elif n_normal > 0.5 * len(result):
        warn(
            f"Normal reference is {n_normal / len(result):.1%} of all cells. This"
            " may be too large for a tumor sample -- including malignant clusters"
            " in the normal reference would cancel out the CNV signal. Check the"
            " score table to confirm each cluster is truly immune/endothelial"
            " (raise --marker-min-margin / --marker-min-within-z, or specify"
            " clusters explicitly with --normal-clusters)."
        )
    return result


def _score_marker_sets(adata: ad.AnnData, marker_sets: dict[str, list[str]]) -> pd.DataFrame:
    scores = pd.DataFrame(index=adata.obs_names)
    for label, genes in marker_sets.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            print(f"[annotate_coarse_celltype] No marker genes found for '{label}': {genes}")
            continue
        score_col = f"_score_{label}"
        sc.tl.score_genes(adata, present, score_name=score_col)
        scores[label] = adata.obs[score_col]
        del adata.obs[score_col]  # don't keep the temporary column
    return scores
def rescue_doublets_by_cnv_consistency(
    doublet_calls: pd.DataFrame,
    cnv_matrix: np.ndarray,
    clone_labels: pd.Series,
    corr_threshold: float = 0.7,
    margin_threshold: float = 0.15,
) -> pd.Series:
    """Rescue doublet-flagged cells whose CNV profile clearly matches a
    single clone, relabeling them as high-RNA malignant singlets rather than
    true doublets.

    Logic:
    - Compute each clone's (e.g. cnv_adata.obs['cnv_leiden']) CNV profile centroid.
    - Correlate each doublet-flagged cell's CNV profile with each clone centroid.
    - If the best-matching clone's correlation is >= corr_threshold, and its
      margin over the second-best clone is >= margin_threshold, treat the
      cell as clearly deriving from a single clone and clear its doublet flag.
    - Conversely, if a cell is similarly close to two clones (small margin),
      keep the doublet call -- likely a true mixture of two clones' RNA.

    A simplified, downstream-only application of the idea behind Numbat's
    joint phylogeny/CNV reconstruction, applied here just to doublet rescue.

    This check must run at single-cell level: after metacell aggregation the
    sum over multiple cells mixes together, and doublet-specific mixture
    signal is lost.
    """
    import scipy.spatial.distance as ssd

    cnv_matrix = _to_dense(cnv_matrix)
    n_cells = len(doublet_calls)
    if cnv_matrix.shape[0] != n_cells or len(clone_labels) != n_cells:
        raise ValueError(
            "doublet_calls / cnv_matrix / clone_labels have mismatched row counts "
            f"({n_cells}, {cnv_matrix.shape[0]}, {len(clone_labels)}). Generate "
            "them from the same single-cell adata with matching row order."
        )
    if not doublet_calls.index.equals(pd.Index(clone_labels.index)):
        raise ValueError(
            "doublet_calls and clone_labels indices do not match. Pass them in "
            "the same adata.obs_names order."
        )

    rescued = doublet_calls["predicted_doublet"].astype(bool).copy()
    clone_ids = pd.Index(pd.Series(clone_labels).astype(str).unique())
    if len(clone_ids) < 2:
        warn(
            f"Only {len(clone_ids)} CNV clone(s) detected; cannot compare between "
            "clones. Skipping doublet rescue."
        )
        return rescued

    labels = pd.Series(clone_labels).astype(str).values
    clone_centroids = np.vstack(
        [np.asarray(cnv_matrix[labels == c]).mean(axis=0) for c in clone_ids]
    )

    doublet_idx = np.where(rescued.values)[0]
    n_rescued = 0
    for i in doublet_idx:
        profile = np.asarray(cnv_matrix[i]).ravel()
        if not np.isfinite(profile).all() or np.std(profile) == 0:
            continue  # correlation is undefined for this cell; don't rescue it
        corrs = np.array(
            [1 - ssd.correlation(profile, centroid) for centroid in clone_centroids]
        )
        corrs = corrs[np.isfinite(corrs)]
        if corrs.size == 0:
            continue
        sorted_corrs = np.sort(corrs)[::-1]
        best = sorted_corrs[0]
        second = sorted_corrs[1] if sorted_corrs.size > 1 else -1.0

        if best >= corr_threshold and (best - second) >= margin_threshold:
            rescued.iloc[i] = False  # rescued: treated as a single-clone-origin singlet
            n_rescued += 1

    log(f"Doublet candidates rescued via CNV consistency: {n_rescued} / {len(doublet_idx)} cells")
    return rescued
# ---------------------------------------------------------------------------
# Step 2. Metacell construction (native implementation, no SEACells dependency)
# ---------------------------------------------------------------------------
# Metacell construction and CNV estimation are implemented natively (see
# metacells_native.py and cnv_native.py) rather than relying on SEACells /
# infercnvpy, so the pipeline no longer breaks when those upstream packages
# change their internals.
#
# Two behavior notes based on measurement:
#  * The line-search step size for the underlying Frank-Wolfe solver was
#    tried and reverted: it lowers RSS but distorts the metacell size
#    distribution (in one real dataset, the largest metacell grew to consume
#    57.8% of cells vs. 3.3% with the fixed step). Neither RSS nor cell-type
#    purity detects this degradation.
#  * Consequently `fw_iters` should be understood as a regularization
#    parameter rather than a speed knob -- see the table in
#    metacells_native.py's docstring for details.
def load_seacell_assignments(
    adata: ad.AnnData,
    source: str,
    metacell_key: str = "SEACell",
    min_coverage: float = 0.80,
    fill_unmatched: bool = True,
) -> int:
    """Load an existing SEACell/metacell assignment, skipping optimization.

    Metacell optimization can take tens of minutes to hours over tens of
    thousands of cells x hundreds of metacells. When only re-running with
    different QC thresholds, the assignment often does not need to be
    rebuilt, so a previous result can be reused.

    `source` is either cell_to_metacell.csv (barcode column + SEACell
    column) or singlecells_qc.h5ad (obs['SEACell']).

    Important: matching is done by **barcode**, not row position -- changing
    QC thresholds changes the cell set, so positional alignment would break
    silently. Cells in the current data with no matching assignment are, if
    fill_unmatched=True, assigned to the metacell of their nearest
    already-assigned neighbor in PCA space (e.g. when relaxing a pctMT
    filter brings cells back).

    Returns: number of cells that received a carried-over assignment.
    """
    path = Path(source)
    if not path.exists():
        raise FileNotFoundError(f"SEACell assignment file does not exist: {source}")

    if path.suffix in (".csv", ".tsv", ".txt"):
        sep = "\t" if path.suffix == ".tsv" else ","
        table = pd.read_csv(path, sep=sep)
        cols = {c.lower(): c for c in table.columns}
        bc_col = cols.get("barcode") or table.columns[0]
        mc_col = cols.get(metacell_key.lower()) or cols.get("seacell")
        if mc_col is None:
            raise ValueError(
                f"No SEACell column found in {source}. Columns: {list(table.columns)}"
            )
        mapping = pd.Series(
            table[mc_col].astype(str).values, index=table[bc_col].astype(str).values
        )
    elif path.suffix == ".h5ad":
        import h5py

        with h5py.File(path, "r") as f:
            o = f["obs"]
            barcodes = [
                b.decode() if isinstance(b, bytes) else str(b)
                for b in o[o.attrs["_index"]][()]
            ]
            if metacell_key not in o:
                raise ValueError(f"{source}'s obs has no '{metacell_key}' column")
            node = o[metacell_key]
            if isinstance(node, h5py.Group):
                cats = [
                    c.decode() if isinstance(c, bytes) else str(c)
                    for c in node["categories"][()]
                ]
                codes = node["codes"][()]
                values = [cats[c] if c >= 0 else None for c in codes]
            else:
                arr = node[()]
                values = [
                    v.decode() if isinstance(v, bytes) else str(v) for v in arr
                ]
        mapping = pd.Series(values, index=barcodes)
    else:
        raise ValueError(f"Unsupported file extension: {path.suffix} (.csv / .h5ad)")

    mapping = mapping[~mapping.index.duplicated()].dropna()
    labels = mapping.reindex(adata.obs_names.astype(str))
    matched = labels.notna().to_numpy()
    coverage = matched.mean() if adata.n_obs else 0.0
    log(
        f"Reusing SEACell assignment: {source}"
        f" ({int(matched.sum()):,}/{adata.n_obs:,} cells matched, {coverage:.1%})"
    )
    if coverage < min_coverage:
        raise ValueError(
            f"Barcode match rate is {coverage:.1%}, below the minimum {min_coverage:.0%}. "
            "This may be an assignment from a different run (different sample or "
            "CellRanger output, or a mix of sample_id-prefixed and raw barcodes). "
            "Stop reusing it and re-run metacell construction."
        )
    if not matched.all():
        n_missing = int((~matched).sum())
        if not fill_unmatched:
            raise ValueError(
                f"{n_missing} cells have no assignment. Enable "
                "--seacell-fill-unmatched, or re-run metacell construction."
            )
        if "X_pca" not in adata.obsm:
            raise KeyError("obsm['X_pca'] is required to fill in unmatched cells")
        from sklearn.neighbors import NearestNeighbors

        emb = adata.obsm["X_pca"]
        nn = NearestNeighbors(n_neighbors=1).fit(emb[matched])
        nearest = nn.kneighbors(emb[~matched], return_distance=False).ravel()
        donor = labels.to_numpy()[matched][nearest]
        filled = labels.to_numpy().copy()
        filled[~matched] = donor
        labels = pd.Series(filled, index=adata.obs_names)
        warn(
            f"Assigned {n_missing} unmatched cells ({n_missing / adata.n_obs:.1%}) to "
            "the metacell of their nearest already-assigned neighbor in PCA space. "
            "This is expected when QC thresholds change the cell set. Metacell "
            "composition will not exactly match the previous run -- re-run metacell "
            "construction to confirm any borderline results."
        )

    adata.obs[metacell_key] = pd.Categorical(labels.astype(str).values)
    n_mc = int(adata.obs[metacell_key].nunique())
    log(f"Number of metacells: {n_mc} (median cells per metacell: {int(labels.value_counts().median())})")
    return int(matched.sum())


def find_cached_assignments(out_dir: Path, extra: list[Path] | None = None) -> Path | None:
    """Look for a reusable SEACell assignment in the output directory (or elsewhere)."""
    candidates = [
        out_dir / "cell_to_metacell.csv",
        out_dir / "singlecells_qc.h5ad",
    ]
    for e in extra or []:
        candidates.extend([e / "cell_to_metacell.csv", e / "singlecells_qc.h5ad"])
    for c in candidates:
        if c.exists():
            return c
    return None


def build_seacells(
    adata: ad.AnnData,
    cells_per_metacell: int = CELLS_PER_METACELL,
    min_iter: int = 10,
    max_iter: int = 50,
    convergence_epsilon: float = 1e-5,
    verbose_iterations: bool = True,
    fw_iters: int | None = None,
    line_search: bool = False,
    init: str = "cssp",
    seed: int = 0,
    **_ignored,
):
    """Build metacells and assign labels to adata.obs['SEACell'].

    Implemented on top of metacells_native.MetacellModel (no SEACells
    dependency). The function name and obs key are kept for compatibility
    with existing outputs and downstream scripts.

    `line_search=True` degrades metacell quality and is disabled if passed
    (see metacells_native.py's docstring for measurements); a warning is
    logged instead.
    """
    import time
    import metacells_native as MC

    if "X_pca" not in adata.obsm:
        raise KeyError(
            "obsm['X_pca'] does not exist. Run PCA before build_seacells()."
        )
    if line_search:
        warn("Disabling exact line search: it distorts the metacell size distribution "
             "(measured: the largest metacell consumed 57.8% of cells). RSS improves "
             "but the result is not usable as metacells")
        line_search = False
    if max_iter < min_iter:
        warn(f"max_iter={max_iter} is smaller than min_iter={min_iter}; "
             f"lowering min_iter to {max_iter}")
        min_iter = max_iter

    n_metacells = max(2, int(np.floor(adata.n_obs / cells_per_metacell)))
    if n_metacells >= adata.n_obs:
        raise ValueError(
            f"Number of metacells {n_metacells} is >= number of cells {adata.n_obs}. "
            "Increase CELLS_PER_METACELL."
        )
    log(f"Building metacells: n_cells={adata.n_obs:,}, n_metacells={n_metacells},"
        f" min_iter={min_iter}, max_iter={max_iter},"
        f" fw_iters={fw_iters or MC.DEFAULT_FW_ITERS}, init={init}")

    model = MC.MetacellModel(
        n_metacells=n_metacells,
        use_rep="X_pca",
        fw_iters=int(fw_iters) if fw_iters else MC.DEFAULT_FW_ITERS,
        line_search=False,
        convergence_epsilon=convergence_epsilon,
        init=init, seed=seed, verbose=True,
    )
    t0 = time.time()
    model.fit(adata, max_iter=max_iter, min_iter=min_iter,
              log_every=1 if verbose_iterations else 10)
    log(f"Optimization complete: {model.n_iter_} iterations / {(time.time() - t0) / 60:.1f} min")
    log(f"Metacell construction complete: {adata.obs['SEACell'].nunique()} metacells")
    return model


def evaluate_metacells(
    adata: ad.AnnData, model=None, celltype_key: str | None = None
) -> pd.DataFrame:
    """Compute compactness / separation / (if celltype_key given) purity.
    Low compactness, high separation, and high purity indicate good metacells.

    Diffusion components are computed natively (no palantir/SEACells
    dependency). Absolute values are not on the same scale as the palantir
    version, so compare across runs by rank rather than absolute value.
    """
    import metacells_native as MC

    metrics = MC.compactness(adata, use_rep="X_pca", key="SEACell").join(
        MC.separation(adata, use_rep="X_pca", key="SEACell", nth_nbr=1), how="outer")
    metrics.index.name = "SEACell"

    if celltype_key is not None and celltype_key in adata.obs:
        try:
            metrics = metrics.join(MC.celltype_purity(adata, celltype_key, key="SEACell"),
                                   how="outer")
        except Exception as exc:
            warn(f"Failed to compute purity (skipping): {exc}")

    bal = MC.size_balance(adata.obs["SEACell"].astype(str).values)
    log(f"Metacell sizes: median {bal['size_median']:.0f}"
        f" [{bal['size_min']}, {bal['size_max']}] / Gini {bal['gini']:.3f}"
        f" / largest share {bal['max_share']:.2%} / singletons {bal['n_singleton']}")
    return metrics.reset_index()


def find_low_quality_metacells(
    metrics: pd.DataFrame,
    compactness_quantile: float = 0.95,
    separation_quantile: float = 0.5,
) -> list[str]:
    """Select candidate low-quality metacells using two dimensions
    (compactness x separation).

    Thresholding compactness alone is risky: compactness and separation are
    often highly rank-correlated in real data, because a metacell in a
    sparse region of the manifold is both "loose internally" (high
    compactness) and "far from its neighbors" (high separation) -- this
    is often simply a rare subpopulation collapsed into one metacell, and
    excluding it risks dropping a rare malignant subclone.

    Truly low-quality metacells are ones that are internally heterogeneous
    *and* still overlap their neighbors, so this selects on the AND of:
      compactness > compactness_quantile quantile, AND
      separation  < separation_quantile quantile
    """
    if metrics.empty or "compactness" not in metrics or "separation" not in metrics:
        return []
    comp_cut = metrics["compactness"].quantile(compactness_quantile)
    sep_cut = metrics["separation"].quantile(separation_quantile)
    mask = (metrics["compactness"] > comp_cut) & (metrics["separation"] < sep_cut)
    bad = metrics.loc[mask, "SEACell"].astype(str).tolist()
    n_comp_only = int((metrics["compactness"] > comp_cut).sum())
    log(
        f"Low-quality metacell candidates: {len(bad)}"
        f" (compactness > {comp_cut:.3g} and separation < {sep_cut:.3g})."
        f" compactness alone would flag {n_comp_only}, but high-separation ones"
        " are kept since they may be rare subpopulations"
    )
    if bad:
        print(
            metrics.loc[mask, ["SEACell", "compactness", "separation"]]
            .sort_values("compactness", ascending=False)
            .to_string(index=False)
        )
    return bad


def _as_seacell_frame(obj) -> pd.DataFrame:
    """Normalize the return value of a metacell evaluation function
    (DataFrame / Series, sometimes indexed by SEACell) into a DataFrame with
    a 'SEACell' column (absorbs version differences)."""
    if isinstance(obj, pd.Series):
        obj = obj.to_frame()
    obj = obj.copy()
    if "SEACell" not in obj.columns:
        obj.index.name = "SEACell"
        obj = obj.reset_index()
    return obj

# ---------------------------------------------------------------------------
# Step 3. Aggregate raw counts to metacells (sum)
# ---------------------------------------------------------------------------

def aggregate_to_metacells(
    adata: ad.AnnData,
    metacell_key: str = "SEACell",
    sample_key: str | None = "sample_id",
    celltype_key: str | None = None,
) -> ad.AnnData:
    """Sum-aggregate raw counts per metacell and return a metacell x gene AnnData.
    obs carries the constituent cell count (n_cells), sample ID (majority
    vote), and cell type (majority vote, if available).

    Implemented with scanpy/pandas only, to avoid a decoupler dependency
    (an alternative would be decoupler's get_pseudobulk).
    """
    import scipy.sparse as sp

    if metacell_key not in adata.obs:
        raise KeyError(f"obs['{metacell_key}'] does not exist. Run build_seacells() first.")
    if "counts" not in adata.layers:
        raise KeyError(
            "layers['counts'] does not exist. load_and_preprocess() is expected to "
            "keep raw counts around (required for aggregation)."
        )

    counts = adata.layers["counts"]
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)

    groups = adata.obs[metacell_key].astype(str)
    metacell_ids = pd.Index(groups.unique())

    agg_matrix = sp.lil_matrix((len(metacell_ids), adata.n_vars))
    obs_rows = []

    for i, mc in enumerate(metacell_ids):
        mask = (groups == mc).values
        agg_matrix[i, :] = counts[mask, :].sum(axis=0)

        row = {"SEACell": mc, "n_cells": int(mask.sum())}

        if sample_key is not None and sample_key in adata.obs:
            # Majority vote within the metacell (metacells are ideally built
            # to stay within a single sample).
            row[sample_key] = str(adata.obs.loc[mask, sample_key].mode().iat[0])
            n_samples_in_mc = adata.obs.loc[mask, sample_key].nunique()
            row["n_samples_in_metacell"] = int(n_samples_in_mc)

        if celltype_key is not None and celltype_key in adata.obs:
            row[celltype_key] = str(adata.obs.loc[mask, celltype_key].mode().iat[0])

        obs_rows.append(row)

    obs_df = pd.DataFrame(obs_rows).set_index("SEACell")
    if "n_samples_in_metacell" in obs_df.columns:
        n_mixed = int((obs_df["n_samples_in_metacell"] > 1).sum())
        if n_mixed:
            warn(
                f"{n_mixed} metacell(s) contain cells from multiple samples. Since "
                "sample_id is assigned by majority vote, DE covariate correction may "
                "be weaker for these."
            )

    mc_adata = ad.AnnData(
        X=agg_matrix.tocsr(),
        obs=obs_df,
        var=adata.var.copy(),
    )
    mc_adata.layers["counts"] = mc_adata.X.copy()
    log(f"Metacell aggregation: {mc_adata.n_obs} metacells x {mc_adata.n_vars} genes")
    return mc_adata


# ---------------------------------------------------------------------------
# Step 3b. Detect abnormal metacells via mitochondrial/nuclear expression balance
# ---------------------------------------------------------------------------

def find_nuclear_mito_genes(adata: ad.AnnData) -> list[str]:
    """Find nuclear-encoded mitochondrial-machinery genes in var_names.

    Always excludes mtDNA-encoded genes (var['mt'] True). COX1-COX3 are
    mtDNA-encoded while COX4I1/COX5A etc. are nuclear-encoded but similarly
    named, so a simple prefix match alone would mix the two; using the
    sequence-ID-derived var['mt'] flag keeps them separated reliably.
    """
    upper = pd.Index(adata.var_names.str.upper())
    is_mito = (
        np.asarray(adata.var["mt"], dtype=bool)
        if "mt" in adata.var
        else np.zeros(adata.n_vars, dtype=bool)
    )
    mask = np.zeros(adata.n_vars, dtype=bool)
    for prefix in NUCLEAR_MITO_PREFIXES:
        mask |= np.array([str(n).startswith(prefix) for n in upper], dtype=bool)
    mask |= np.isin(np.asarray(upper, dtype=object), list(NUCLEAR_MITO_GENES))
    mask &= ~is_mito
    return list(adata.var_names[mask])


def detect_expression_imbalance(
    mc_adata: ad.AnnData,
    sc_adata: ad.AnnData | None = None,
    metacell_key: str = "SEACell",
    nmads: float = 3.5,
) -> pd.DataFrame:
    """Evaluate mtDNA/nuclear expression balance per metacell and flag anomalies.

    Excluding mitochondrial genes from the CNV signal discards the
    information that something unusual happened in that metacell. This
    compares each metacell to the rest along four axes, all using the same
    MAD-based adaptive thresholding as the rest of the pipeline (no fixed
    cutoffs):

    1. High mt_frac (fraction of counts from mtDNA)
       -> a metacell enriched for dying/apoptotic cells. Cells that passed
          single-cell QC individually can still concentrate into one metacell.
    2. Low mt_frac
       -> candidate for nucleus-only capture, or ambient-RNA-dominated metacell.
    3. mtDNA-encoded / nuclear-encoded mitochondrial-machinery count ratio
       deviates from the trend across metacells
       -> OXPHOS complexes combine mtDNA- and nuclear-encoded subunits in
          fixed stoichiometry, so this ratio should be roughly constant;
          a deviation suggests mitochondrial dysfunction or altered mtDNA
          copy number, which plain mt_frac would miss. Computed as the
          direct log-ratio between the two count totals (log2((mito+1)/(nuc+1))),
          NOT via regressing log2(mt_frac) on log2(nuclear fraction) -- the
          latter is wrong because both are fractions of the same total count
          and are therefore structurally anti-correlated, which lets the
          regression line absorb the very deviation this is meant to detect.
    4. High per-cell pctMT variance within a metacell (only if sc_adata is given)
       -> aggregation accidentally mixed healthy and dying cells; invisible
          from the metacell-level mean alone, so single-cell-level variance
          is checked explicitly.

    Returns: a DataFrame of metrics and flags, indexed like mc_adata.obs.
    """
    import scipy.sparse as sp

    if "mt" not in mc_adata.var:
        warn(
            "var['mt'] is missing. Skipping mitochondrial/nuclear balance detection "
            "(specify mtDNA via --mito-chromosome)."
        )
        return pd.DataFrame(index=mc_adata.obs_names)

    counts = mc_adata.layers["counts"] if "counts" in mc_adata.layers else mc_adata.X
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)

    is_mito = np.asarray(mc_adata.var["mt"], dtype=bool)
    nuc_genes = find_nuclear_mito_genes(mc_adata)
    is_nuc = np.asarray(mc_adata.var_names.isin(nuc_genes), dtype=bool)

    log("=== Mitochondrial/nuclear expression balance check ===")
    log(f"mtDNA-encoded genes: {int(is_mito.sum())} / nuclear-encoded machinery genes: {len(nuc_genes)}")
    if int(is_mito.sum()) == 0:
        warn("No mtDNA-encoded genes found; cannot run this check.")
        return pd.DataFrame(index=mc_adata.obs_names)

    total = np.asarray(counts.sum(axis=1)).ravel()
    mito = np.asarray(counts[:, is_mito].sum(axis=1)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        mt_frac = np.where(total > 0, mito / total, np.nan)

    res = pd.DataFrame(index=mc_adata.obs_names)
    res["mt_frac"] = mt_frac
    res["pct_counts_mt_metacell"] = mt_frac * 100

    # --- Axis 1/2: mt_frac outliers ---
    finite = np.isfinite(mt_frac)
    res["mito_high"] = False
    res["mito_low"] = False
    if finite.sum() >= 10:
        res.loc[finite, "mito_high"] = mad_outlier_mask(
            mt_frac[finite], nmads=nmads, log=True, direction="higher"
        )
        res.loc[finite, "mito_low"] = mad_outlier_mask(
            mt_frac[finite], nmads=nmads, log=True, direction="lower"
        )
    else:
        warn("Fewer than 10 metacells; skipping mt_frac outlier detection")

    # --- Axis 3: outliers in the mtDNA-encoded / nuclear-encoded log2 ratio ---
    #
    # Note: this compares the two count totals directly via a log-ratio
    # (compositional-data style), not by regressing one fraction against the
    # other -- see the docstring above for why the naive regression fails to
    # detect a real deviation.
    res["nuc_mito_frac"] = np.nan
    res["mito_nuc_log2ratio"] = np.nan
    res["ratio_outlier"] = False
    if len(nuc_genes) >= 10:
        nuc = np.asarray(counts[:, is_nuc].sum(axis=1)).ravel()
        with np.errstate(divide="ignore", invalid="ignore"):
            nuc_frac = np.where(total > 0, nuc / total, np.nan)
        res["nuc_mito_frac"] = nuc_frac
        # +1 avoids log blowing up when one side is 0 for a metacell
        log_ratio = np.log2((mito + 1.0) / (nuc + 1.0))
        ok = np.isfinite(log_ratio)
        if ok.sum() >= 10:
            res.loc[ok, "mito_nuc_log2ratio"] = log_ratio[ok]
            res.loc[ok, "ratio_outlier"] = mad_outlier_mask(
                log_ratio[ok], nmads=nmads, log=False, direction="both"
            )
            med = float(np.median(log_ratio[ok]))
            mad = float(np.median(np.abs(log_ratio[ok] - med)) * 1.4826)
            log(
                f"log2(mtDNA-encoded / nuclear-encoded): median {med:.3f}, MAD {mad:.3f}"
                f" -> thresholds {med - nmads * mad:.3f} to {med + nmads * mad:.3f}"
            )
        else:
            warn("Fewer than 10 valid metacells; skipping ratio outlier detection")
    else:
        warn(
            f"Only {len(nuc_genes)} nuclear-encoded mitochondrial-machinery genes "
            "found (fewer than 10); skipping ratio-based detection. Check whether "
            "this species uses NDUF*/MRPL*-style symbols."
        )

    # --- Axis 4: per-cell pctMT variance within a metacell ---
    res["mt_frac_sd_within"] = np.nan
    res["mito_heterogeneous"] = False
    if (
        sc_adata is not None
        and metacell_key in sc_adata.obs
        and "pct_counts_mt" in sc_adata.obs
    ):
        grp = sc_adata.obs.groupby(sc_adata.obs[metacell_key].astype(str), observed=True)[
            "pct_counts_mt"
        ]
        sd = grp.std().reindex(res.index.astype(str))
        res["mt_frac_sd_within"] = sd.values
        ok = np.isfinite(res["mt_frac_sd_within"].values)
        if ok.sum() >= 10:
            res.loc[ok, "mito_heterogeneous"] = mad_outlier_mask(
                res["mt_frac_sd_within"].values[ok], nmads=nmads, log=True, direction="higher"
            )
    else:
        log("No single-cell data given; skipping within-metacell variance check")

    # --- Combine ---
    flags = ["mito_high", "mito_low", "ratio_outlier", "mito_heterogeneous"]
    res["mito_anomaly"] = res[flags].any(axis=1)
    reason_map = {
        "mito_high": "high mtDNA fraction (possible dying cells)",
        "mito_low": "low mtDNA fraction (possible nucleus-only/ambient)",
        "ratio_outlier": "mtDNA/nuclear ratio deviates",
        "mito_heterogeneous": "pctMT heterogeneous within metacell",
    }
    res["mito_anomaly_reason"] = [
        "; ".join(reason_map[f] for f in flags if row[f]) if row["mito_anomaly"] else ""
        for _, row in res.iterrows()
    ]

    n_flag = int(res["mito_anomaly"].sum())
    log(
        f"Anomalous metacells: {n_flag} / {len(res)} ({n_flag / max(len(res), 1):.1%})"
        f" -- breakdown {{{', '.join(f'{f}: {int(res[f].sum())}' for f in flags)}}}"
    )
    if n_flag:
        show = res.loc[
            res["mito_anomaly"],
            ["pct_counts_mt_metacell", "mito_nuc_log2ratio", "mt_frac_sd_within",
             "mito_anomaly_reason"],
        ].sort_values("pct_counts_mt_metacell", ascending=False)
        print(show.head(20).round(3).to_string())
    if n_flag / max(len(res), 1) > 0.3:
        warn(
            "Over 30% of metacells were flagged. This can happen with the MAD-based "
            "relative thresholds if the underlying distribution is itself bimodal "
            "(e.g. the sample contains two distinct cell states). Consider raising "
            "nmads, or inspect the distribution visually."
        )
    return res


def report_mito_confounding(
    mc_obs: pd.DataFrame,
    anomaly_key: str = "mito_anomaly",
    condition_key: str = "putative_malignant",
) -> None:
    """Test whether mitochondrial anomalies are confounded with the malignant label.

    This is the most important use of this check: if anomalous metacells
    are skewed toward the malignant (or normal) side, DE may end up
    capturing "dying vs. healthy cells" rather than "tumor vs. normal" --
    and since DE runs over all genes, excluding mitochondrial genes from CNV
    does not protect it from this confound.
    """
    if anomaly_key not in mc_obs or condition_key not in mc_obs:
        return
    tab = pd.crosstab(mc_obs[anomaly_key], mc_obs[condition_key])
    log("=== Mitochondrial-anomaly / malignant-label confounding check ===")
    print(tab.to_string())
    if tab.shape != (2, 2):
        log("Not a 2x2 table; skipping the test")
        return
    try:
        from scipy.stats import fisher_exact

        odds, pval = fisher_exact(tab.values)
        log(f"Fisher's exact test: odds ratio {odds:.3g}, p = {pval:.3g}")
        if pval < 0.05:
            warn(
                "Mitochondrial anomalies are significantly skewed with respect to the "
                "malignant label. DE differences may reflect cell state (dying/stressed) "
                "rather than tumor biology. Consider re-running DE with anomalous "
                "metacells excluded (via the mito_anomaly column) and checking whether "
                "the result holds."
            )
        else:
            log("The skew is not significant; confounding of DE appears limited.")
    except Exception as exc:  # pragma: no cover
        log(f"Skipped the test: {exc}")


def add_genomic_positions(
    cnv_adata: ad.AnnData, gtf_path: str, gtf_gene_id: str = GTF_GENE_ID_ATTR
) -> None:
    """Attach chromosome/start/end to adata.var.

    Uses a native GTF attribute parser rather than infercnvpy's
    GENCODE-oriented one, which fails to extract attribute keys from
    NCBI/RefSeq or CellRanger-filtered GTFs.
    """
    if not Path(gtf_path).exists():
        raise FileNotFoundError(f"GTF file does not exist: {gtf_path}")

    gene_pos = parse_gtf_gene_positions(gtf_path, gene_id_type=gtf_gene_id)
    for col in ("chromosome", "start", "end"):
        if col in cnv_adata.var.columns:
            del cnv_adata.var[col]
    cnv_adata.var = cnv_adata.var.join(gene_pos, how="left")
    n = int(cnv_adata.var["chromosome"].notna().sum())
    log(f"Attached coordinates: {n:,}/{cnv_adata.n_vars:,} genes"
        f" / {cnv_adata.var['chromosome'].nunique()} chromosomes")

def run_cnv_branch(
    mc_adata: ad.AnnData,
    gtf_path: str = GTF_PATH,
    gtf_gene_id: str = GTF_GENE_ID_ATTR,
    reference_key: str | None = "cell_type",
    reference_cat: list[str] | None = None,
    window_size: int = 100,
    step: int = 10,
    exclude_chromosomes: tuple[str, ...] = EXCLUDE_CHROMOSOMES,
    chromosome_map: dict[str, str] | None = None,
    mito_chromosomes: list[str] | None = None,
    max_missing_rate: float = 0.8,
    method: str = "infercnv",
    bin_genes: int = 100,
) -> ad.AnnData:
    """Run CNV estimation on a metacell x gene (or single-cell x gene) AnnData.

    Notes:
    - Regardless of whether the input adata.X is already normalized/log
      transformed, normalization and log-transform are always redone from
      layers['counts'] (raw counts), so behavior does not depend on the
      caller's preprocessing state.
    - Requires gene chromosome coordinates (chromosome, start, end) in
      adata.var; add via add_genomic_positions() from the GTF.
    - If reference_key is given, that obs column and the reference_cat
      values are validated to exist before use.
    - Warns if exclude_chromosomes does not match the real chromosome names.
    - `method="infercnv"` (default) is a windowed running-mean approach, a
      native re-implementation matching infercnvpy's formula (verified:
      correlation 1.000000, mean absolute difference 5.4e-09 against
      infercnvpy on real data).
    - `method="bins"` is an alternative that tests genomic-bin count
      composition with a beta-binomial model. Reference-group values are
      not floored at 0, so tumor contamination in the reference remains
      visible; always reports the ratio against a null control (genes
      randomly assigned to bins). On real-data validation this correlated
      slightly less well with verified karyotypes than the windowed
      approach (0.929 vs. 0.969), so it is not the default.
    - Both methods estimate CNV from expression, so cross-checking important
      conclusions against an independent method (e.g. CopyKAT/SCEVAN) is
      recommended.
    """
    # Validate inputs before importing anything heavy, so config mistakes surface fast.
    if "counts" not in mc_adata.layers:
        raise KeyError(
            "layers['counts'] does not exist. run_cnv_branch always renormalizes "
            "from raw counts."
        )

    # --- Validate reference_key ---
    if reference_key is not None:
        if reference_key not in mc_adata.obs:
            raise KeyError(
                f"reference_key='{reference_key}' is not in obs. Run "
                "annotate_coarse_celltype() before the CNV branch. Current obs "
                f"columns: {list(mc_adata.obs.columns)}"
            )
        available = set(map(str, mc_adata.obs[reference_key].unique()))
        if reference_cat is not None:
            # Expand by prefix to also match lineage-suffixed labels (e.g. 'Myeloid_3')
            expanded: list[str] = []
            for c in reference_cat:
                if c in available:
                    expanded.append(c)
                expanded += sorted(v for v in available if v.startswith(str(c) + "_"))
            reference_cat = list(dict.fromkeys(expanded)) or list(reference_cat)
            present_cat = [c for c in reference_cat if c in available]
            if not present_cat:
                raise ValueError(
                    f"None of reference_cat={reference_cat} exist in "
                    f"obs['{reference_key}'] (actual values: {sorted(available)}). "
                    "Cannot secure a normal reference. Based on the score table from "
                    "annotate_coarse_celltype, choose one of:\n"
                    "  - stromal/immune clusters are known -> specify explicitly with"
                    " --normal-clusters '3,7'\n"
                    f"  - tissue type is known and it's OK to use"
                    f" {'/'.join(MALIGNANCY_CAPABLE_CELLTYPES)} as reference ->"
                    " --normal-celltype Fibroblast\n"
                    "  - thresholds were too strict -> lower --marker-min-margin /"
                    " --marker-min-within-z\n"
                    "  - sample has almost no normal cells -> --cnv-refine (use"
                    " CNV-flat clones as a second-pass reference) or --cnv-reference"
                    " none (mean over all metacells)"
                )
            if len(present_cat) < len(reference_cat):
                warn(
                    "Removed reference_cat values that are not present: "
                    f"{[c for c in reference_cat if c not in available]}"
                )
            reference_cat = present_cat
            log(f"CNV normal reference categories: {reference_cat}")
            if len(reference_cat) >= 2:
                log(
                    "  2+ reference categories: using bounded mode (only counts logFC"
                    " outside the reference's min/max, suppressing cell-type-specific"
                    " expression bias)"
                )
            else:
                warn(
                    "Only one reference category. This collapses to a plain "
                    "difference, making cell-type-specific expression bias more "
                    "likely to be mistaken for CNV. Consider using --cnv-refine."
                )

    cnv_adata = mc_adata.copy()

    # Redo normalization/log-transform from raw counts (independent of the input's .X state)
    raw_counts = cnv_adata.layers["counts"]
    if not looks_like_raw_counts(raw_counts):
        warn(
            "layers['counts'] does not look like non-negative integers. Verify that "
            "raw counts were preserved correctly."
        )
    cnv_adata.X = raw_counts.copy()
    # If the caller's adata was already log1p'd, uns['log1p'] would carry over and
    # cause scanpy to falsely warn "adata.X seems to be already log-transformed".
    # Since .X is rebuilt from raw counts here, drop that history before re-applying.
    cnv_adata.uns.pop("log1p", None)
    sc.pp.normalize_total(cnv_adata)
    sc.pp.log1p(cnv_adata)

    # --- Attach gene coordinate annotation ---
    if not {"chromosome", "start", "end"}.issubset(cnv_adata.var.columns):
        add_genomic_positions(cnv_adata, gtf_path, gtf_gene_id=gtf_gene_id)

    # Genes that didn't match the GTF (chromosome is NaN) are excluded from the
    # infercnv calculation. If the match rate is low, check gtf_gene_id (typically
    # "gene" for NCBI-style GTFs, "gene_name" for GENCODE-style) and whether
    # var_names uses a matching naming convention.
    n_missing = int(cnv_adata.var["chromosome"].isna().sum())
    missing_rate = n_missing / cnv_adata.n_vars if cnv_adata.n_vars else 1.0
    log(f"Genes with no chromosome position: {n_missing}/{cnv_adata.n_vars} ({missing_rate:.1%})")
    if missing_rate > max_missing_rate:
        raise ValueError(
            f"Could not assign chromosome positions for {missing_rate:.1%} of genes. "
            f"Check the gtf_gene_id='{gtf_gene_id}' setting, or var_names naming "
            "convention."
        )

    # Normalize chromosome naming (the CNV estimator only processes chromosomes
    # whose name starts with 'chr')
    cnv_adata.var = normalize_chromosome_names(
        cnv_adata.var, chromosome_map, mito_chromosomes=mito_chromosomes
    )
    n_mito_var = int((cnv_adata.var["chromosome"].astype(str) == "chrM").sum())
    if n_mito_var:
        log(f"Excluding {n_mito_var} mitochondrial gene(s) ('chrM') from the CNV calculation")

    actual_chroms = set(map(str, cnv_adata.var["chromosome"].dropna().unique()))
    if not actual_chroms:
        raise ValueError("No genes have a chromosome position assigned. Check the GTF settings.")

    # A chromosome with fewer genes than window_size degenerates to a single window
    per_chrom = cnv_adata.var["chromosome"].value_counts()
    small_chroms = per_chrom[per_chrom < window_size]
    if len(small_chroms):
        warn(
            f"window_size={window_size} exceeds the gene count on "
            f"{len(small_chroms)} chromosome(s) (e.g. {small_chroms.head(3).to_dict()}). "
            "CNV on these chromosomes degenerates to the all-gene average. Consider "
            "lowering window_size."
        )

    # Check that exclude_chromosomes actually matches real chromosome names
    effective_exclude = check_exclude_chromosomes(tuple(exclude_chromosomes), actual_chroms)

    # CNV estimation is implemented natively in cnv_native.infercnv_scores using
    # the same formula as infercnvpy (validated: correlation 1.000000, mean
    # absolute difference 5.4e-09 on real data, ~2x faster). The one behavioral
    # difference is that the dynamic threshold's SD is computed once over all
    # cells rather than per-chunk.
    import cnv_native as _CNVN

    if method == "bins":
        ref_mask = None
        if reference_key and reference_cat and reference_key in cnv_adata.obs:
            ref_mask = cnv_adata.obs[reference_key].astype(str).isin(
                [str(c) for c in reference_cat]).values
        _res = _CNVN.bin_composition_cnv(
            cnv_adata, bin_genes=bin_genes, reference_mask=ref_mask,
            layer="counts",
            exclude_chromosomes=tuple(effective_exclude) if effective_exclude else (),
        )
        log(f"Bin-composition method: amplitude {_res['amplitude']:.4f}"
            f" / expression-program control {_res.get('control_amplitude', float('nan')):.4f}"
            f" -> ratio {_res.get('amplitude_ratio', float('nan')):.2f}x")
    else:
        _CNVN.infercnv_scores(
            cnv_adata,
            reference_key=reference_key,   # obs column name (known malignant/normal annotation)
            reference_cat=reference_cat,   # values within reference_key treated as "normal"
            window_size=window_size,
            step=step,
            exclude_chromosomes=tuple(effective_exclude) if effective_exclude else (),
        )
    # Result is stored in cnv_adata.obsm["X_cnv"]

    # Attach the L2 norm as an aneuploidy-burden CNV score for downstream use
    cnv_matrix = cnv_adata.obsm["X_cnv"]
    if hasattr(cnv_matrix, "multiply"):
        squared_sum = np.asarray(cnv_matrix.multiply(cnv_matrix).sum(axis=1))
    else:
        squared_sum = np.asarray(np.square(np.asarray(cnv_matrix)).sum(axis=1))
    cnv_adata.obs["cnv_score"] = np.sqrt(squared_sum).ravel()

    # Embedding/clustering on the CNV matrix is just scanpy applied to a
    # temporary AnnData with obsm['X_cnv'] as X, so call scanpy directly.
    _CNVN.cnv_embedding(cnv_adata, key="cnv", cluster_key="cnv_leiden")
    log(f"Candidate CNV clones (cnv_leiden): {cnv_adata.obs['cnv_leiden'].nunique()}")

    return cnv_adata


# ---------------------------------------------------------------------------
# Step 4a-1b. Clone lineage consistency (agreement with a single-origin expectation)
# ---------------------------------------------------------------------------

def chromosome_level_profiles(
    cnv_adata: ad.AnnData,
    clone_key: str = "cnv_leiden",
    exclude_unplaced: bool = True,
) -> pd.DataFrame:
    """Return per-clone chromosome-level CNV profiles (chromosome x clone).

    Window-level X_cnv is noisy, so this averages within each chromosome for
    comparison. Unplaced scaffolds (NW_*) usually have only a few windows
    and are prone to being outliers, so they are excluded by default.
    """
    X = cnv_adata.obsm["X_cnv"]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    chr_pos = cnv_adata.uns["cnv"]["chr_pos"]
    order = sorted(chr_pos.items(), key=lambda kv: kv[1])
    bounds = [
        (k, v, (order[i + 1][1] if i + 1 < len(order) else X.shape[1]))
        for i, (k, v) in enumerate(order)
    ]
    if exclude_unplaced:
        bounds = [t for t in bounds if "NW_" not in t[0]]
    labels = cnv_adata.obs[clone_key].astype(str)
    data = {}
    for g in sorted(labels.unique(), key=lambda v: (len(v), v)):
        m = (labels == g).values
        data[g] = np.array([X[m, st:en].mean() for _, st, en in bounds])
    return pd.DataFrame(data, index=[k for k, _, _ in bounds])


def write_interpretation_caveats(
    out_dir: Path,
    mito_comp: dict | None = None,
    lineage: dict | None = None,
    use_pctmt: bool = True,
    malignant_call: str = "",
    sample_ids: list[str] | None = None,
) -> Path:
    """Write run-specific interpretation caveats to INTERPRETATION_CAVEATS.txt.

    Why this file exists: data-quality issues (e.g. a rescue alignment, an
    insufficient reference) show up in the run log, but that context is lost
    for anyone who later looks only at the result files, or compares them
    against other samples. This bundles detected issues with the output so
    that data of differing quality/provenance is not misread under the same
    assumptions.
    """
    lines: list[str] = []
    lines.append("Interpretation caveats (specific to this run)")
    lines.append("=" * 60)
    lines.append(f"Generated by: metacellcnv.py / sample_id = {sample_ids}")
    lines.append("")
    lines.append("IMPORTANT: this file's content differs per run. When comparing results")
    lines.append("across multiple samples, read each sample's copy of this file and confirm")
    lines.append("they can be interpreted under the same assumptions. Do not interpret a")
    lines.append("sample with data-quality issues under the same standard as a clean one.")
    lines.append("")

    severe = False
    if mito_comp and mito_comp.get("status") == "suspect":
        severe = True
        lines.append("[SEVERE] mtDNA gene composition looks biologically implausible")
        for msg in mito_comp.get("problems", []):
            lines.append(f"  - {msg}")
        lines.append("  Possible causes: alignment against a reference excluding mtDNA,")
        lines.append("  sequence rescue from BAM, NUMT mismapping, or GTF/BAM mismatch.")
        lines.append("  Impact: pctMT does not reflect true mitochondrial content.")
        lines.append("    - pctMT-based QC is meaningless (does not remove dying cells)")
        lines.append("    - the mitochondrial/nuclear balance anomaly flag (mito_anomaly)")
        lines.append("      is picking up technical artifact level, not dying cells")
        lines.append("    - if mtDNA genes enter the HVG set, even the neighbor graph is")
        lines.append("      shaped by this technical factor (excluded from HVG by default)")
        lines.append("  -> This sample cannot be interpreted the same way as a normal sample.")
        lines.append("")
    elif mito_comp and mito_comp.get("status") == "ok":
        lines.append("[OK] mtDNA gene composition is consistent with known transcript patterns")
        lines.append("")

    if not use_pctmt:
        lines.append("[SETTING] pctMT-based cell exclusion was disabled for this run")
        lines.append("  (--no-pctmt-filter). The high_pctmt flag is recorded but not used to filter.")
        lines.append("")

    if lineage:
        st = lineage.get("status")
        if st == "consistent":
            lines.append("[OK] Malignant clones are consistent with a single-origin expectation")
            lines.append(
                f"  malignant clones: {lineage.get('n_malignant_clones')} / "
                f"affected chromosomes: {lineage.get('n_affected_chromosomes')}"
            )
            if np.isfinite(lineage.get("sign_agreement", float("nan"))):
                lines.append(
                    f"  sign agreement across clones: {lineage['sign_agreement']:.1%} / "
                    f"minimum inter-clone correlation: {lineage.get('min_correlation', float('nan')):.3f}"
                )
            nr, nt = lineage.get("normal_residual_pairs"), lineage.get("normal_residual_total")
            if nr is not None:
                lines.append(f"  residual deviation in normal clones (negative control): {nr} / {nt}")
            lines.append("")
        elif st == "suspect":
            severe = True
            lines.append("[SEVERE] Malignant clones conflict with a single-origin expectation")
            for msg in lineage.get("problems", []):
                lines.append(f"  - {msg}")
            lines.append("  Tumor cell origin is normally very limited, so multiple lineages")
            lines.append("  with mutually unrelated CNV profiles are not expected.")
            lines.append("  Some may be false positives from cell-type-specific expression bias.")
            lines.append("")

    if malignant_call:
        lines.append("[Call basis] " + malignant_call)
        lines.append("")

    lines.append("Always-applicable caveats:")
    lines.append("  - All metacells in a single sample are pseudoreplicates; DE p-values should not be used for inference")
    lines.append("  - infercnvpy's developers describe it as experimental; cross-check")
    lines.append("    important conclusions against CopyKAT / SCEVAN or similar tools")
    lines.append("  - CNV is estimated from expression, so DE between CNV-defined groups")
    lines.append("    retains some circularity in principle")
    lines.append("")
    lines.append(
        "Overall: " + ("issues found -- do not interpret alongside other samples under the"
                        " same standard" if severe else "no serious issues detected")
    )

    path = out_dir / "INTERPRETATION_CAVEATS.txt"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"Saved interpretation caveats: {path.name}"
        + (" (severe issues found)" if severe else ""))
    return path


def check_clone_lineage_consistency(
    cnv_adata: ad.AnnData,
    malignant_clones: list[str],
    normal_clones: list[str],
    clone_key: str = "cnv_leiden",
    deviation_threshold: float = 0.02,
    min_correlation: float = 0.5,
    out_path: Path | None = None,
) -> dict:
    """Test whether clones called malignant are consistent with a
    single-origin expectation.

    Biological rationale: a tumor's cells normally have a very limited
    number of origins. Even if subclones diverge in the magnitude of
    amplification/deletion, **the affected chromosomes themselves are
    expected to be shared**. Conversely, if many mutually unrelated
    "tumor lineages" are found, many are likely false positives (cell-type-
    specific expression bias mistaken for CNV).

    Two tests:
      1. Cross-correlation of chromosome-level profiles between malignant
         clones. Low correlation implies unrelated lineages coexisting,
         which conflicts with a single-origin expectation.
      2. Among chromosomes whose |deviation| exceeds the threshold, the
         fraction where the **sign** (gain/loss) agrees across all
         malignant clones. Should be near 100% under single origin.

    The residual deviation in normal clones is also reported as a negative
    control -- a large value suggests a poorly chosen reference, prone to
    false positives.
    """
    log("=== Clone lineage consistency (agreement with single-origin expectation) ===")
    if not malignant_clones:
        log("No malignant clones; skipping this test")
        return {"status": "skipped"}

    P = chromosome_level_profiles(cnv_adata, clone_key=clone_key)
    missing = [g for g in malignant_clones + normal_clones if g not in P.columns]
    if missing:
        warn(f"Could not build a profile for some clones: {missing}")
    malignant_clones = [g for g in malignant_clones if g in P.columns]
    normal_clones = [g for g in normal_clones if g in P.columns]
    if not normal_clones:
        warn("No normal clones; cannot establish a baseline. Using relative comparison only.")
        base = np.zeros(len(P))
    else:
        base = P[normal_clones].mean(axis=1).values

    D = P[malignant_clones].sub(base, axis=0)
    result: dict = {"n_malignant_clones": len(malignant_clones)}

    if len(malignant_clones) >= 2:
        corr = D.corr()
        result["min_correlation"] = float(
            corr.values[np.triu_indices_from(corr.values, k=1)].min()
        )
        log("Chromosome-level profile correlation between malignant clones:")
        print(corr.round(3).to_string())
    else:
        result["min_correlation"] = float("nan")
        log("Only one malignant clone; skipping inter-clone correlation (single-clone tumor)")

    hit = D[(D.abs() > deviation_threshold).any(axis=1)]
    result["n_affected_chromosomes"] = int(len(hit))
    if len(hit) and len(malignant_clones) >= 2:
        sign = np.sign(hit)
        agree = float(sign.eq(sign.iloc[:, 0], axis=0).all(axis=1).mean())
        result["sign_agreement"] = agree
        log(
            f"Of the chromosomes with |deviation| > {deviation_threshold} ({len(hit)} total), "
            f"gain/loss direction agrees across all malignant clones: {agree:.1%}"
        )
    else:
        result["sign_agreement"] = float("nan")

    if normal_clones:
        N = P[normal_clones].sub(base, axis=0)
        n_resid = int((N.abs() > deviation_threshold).sum().sum())
        result["normal_residual_pairs"] = n_resid
        result["normal_residual_total"] = int(N.size)
        log(
            f"Residual deviation in normal clones (negative control): "
            f"(clone, chromosome) pairs with |deviation|>{deviation_threshold} = {n_resid} / {N.size}"
        )
        if N.size and n_resid / N.size > 0.05:
            warn(
                "Normal clones still show residual deviation. The reference may be "
                "poorly chosen, and cell-type-specific expression bias may be "
                "mistaken for CNV. Consider re-running with --cnv-refine, supplying "
                "normal clones as multiple reference categories."
            )

    if len(hit):
        log("Top affected chromosomes (by max absolute deviation, top 10):")
        print(hit.reindex(hit.abs().max(axis=1).sort_values(ascending=False).index)
              .head(10).round(4).to_string())

    # Verdict
    problems = []
    if len(malignant_clones) >= 2:
        if result["min_correlation"] < min_correlation:
            problems.append(
                f"Minimum correlation between malignant clones is "
                f"{result['min_correlation']:.3f}, below threshold {min_correlation}. "
                "Lineages with mutually unrelated CNV profiles coexist, conflicting "
                "with a single-origin expectation; some may be false positives"
            )
        if np.isfinite(result["sign_agreement"]) and result["sign_agreement"] < 0.7:
            problems.append(
                f"Gain/loss sign agreement is only {result['sign_agreement']:.1%}. "
                "Clones without a shared CNV event appear to be mixed together"
            )
    result["problems"] = problems
    result["status"] = "suspect" if problems else "consistent"
    if problems:
        warn("Lineage consistency issues found:")
        for i, m in enumerate(problems, 1):
            warn(f"  ({i}) {m}")
        warn(
            "Options: (a) raise --malignant-n-sd for a stricter call, "
            "(b) improve the reference with --cnv-refine, "
            "(c) inspect low-correlation clones individually via a chromosome heatmap."
        )
    else:
        log(
            "Consistent with a single-origin expectation (malignant clones share the "
            "same direction of change on the same chromosomes)"
        )

    if out_path is not None:
        D.round(5).to_csv(out_path)
        log(f"Saved chromosome-level profiles: {out_path.name}")
    return result


# ---------------------------------------------------------------------------
# Step 4a-2. Calling malignant metacells
# ---------------------------------------------------------------------------

def call_malignant_metacells(
    cnv_obs: pd.DataFrame,
    method: str = "auto",
    clone_key: str = "cnv_leiden",
    score_key: str = "cnv_score",
    celltype_key: str | None = "cell_type",
    normal_celltypes: list[str] | None = None,
    n_sd: float = 3.0,
    min_gap_frac: float = 0.25,
) -> tuple[pd.Series, str]:
    """Call malignant metacells from cnv_score / cnv_leiden.

    Rationale: thresholding directly on cnv_score (e.g. score > its median)
    only splits the data roughly 50:50 by construction and does not detect
    tumor cells -- it does not distinguish "high vs. low aneuploidy burden"
    from "tumor vs. normal". It also arbitrarily splits a single clone if
    its scores straddle the median.

    This function instead calls malignancy per **clone** (cnv_leiden): every
    metacell in the same clone gets the same label, so there is no arbitrary
    split near a threshold.

    method:
      "reference" -- Threshold = mean + n_sd*SD of the cnv_score distribution
                    of normal-reference metacells (cell_type in
                    normal_celltypes). A clone is malignant if its
                    **clone median** exceeds this threshold. Most
                    interpretable since it's anchored to real normal cells.
                    Requires at least 3 normal-reference metacells.
      "gap"       -- Used when there's no normal reference. Sort per-clone
                    median cnv_score ascending and split into two groups at
                    the **largest gap** between adjacent medians (the
                    natural 1-D break point). If the largest gap is less
                    than min_gap_frac of the median range, there is no clear
                    bimodality and the call is abandoned (all "unassigned").
      "median"    -- Legacy median split. Kept for comparison/reproducibility
                    only; not used by default.
      "auto"      -- "reference" if a normal reference is available, else "gap".

    Returns: (label Series, description of how the call was made)

    Note: since CNV is itself estimated from expression, DE between
    CNV-defined groups retains some circularity in principle. Per-clone
    assignment and anchoring to a normal reference mitigate this but do not
    eliminate it -- cross-check important conclusions against an independent
    method (e.g. CopyKAT/SCEVAN).
    """
    if normal_celltypes is None:
        normal_celltypes = KNOWN_NORMAL_CELLTYPES
    if score_key not in cnv_obs:
        raise KeyError(f"obs['{score_key}'] does not exist.")

    scores = cnv_obs[score_key].astype(float)
    clones = (
        cnv_obs[clone_key].astype(str)
        if clone_key in cnv_obs
        else pd.Series(["0"] * len(cnv_obs), index=cnv_obs.index)
    )
    clone_median = scores.groupby(clones).median().sort_values()

    log("=== Calling malignant metacells ===")
    log(f"Per-clone median cnv_score (ascending):\n{clone_median.round(4).to_string()}")

    # --- Check whether a normal reference is available ---
    n_ref = 0
    if celltype_key is not None and celltype_key in cnv_obs:
        ref_mask = pd.Series(
            is_known_normal(cnv_obs[celltype_key]), index=cnv_obs.index
        )
        n_ref = int(ref_mask.sum())
    else:
        ref_mask = pd.Series(False, index=cnv_obs.index)

    if method == "auto":
        method = "reference" if n_ref >= 3 else "gap"
        log(f"Auto-selected method: '{method}' ({n_ref} normal-reference metacells)")

    if method == "reference":
        if n_ref < 3:
            warn(
                f"Only {n_ref} normal-reference metacell(s); the reference method "
                "cannot be used. Falling back to the gap method."
            )
            method = "gap"
        else:
            ref_scores = scores[ref_mask]
            mu, sd = float(ref_scores.mean()), float(ref_scores.std(ddof=1))
            threshold = mu + n_sd * sd if sd > 0 else float(ref_scores.max())
            malignant_clones = clone_median.index[clone_median > threshold].tolist()
            desc = (
                f"reference method: {n_ref} normal-reference metacells' cnv_score = "
                f"{mu:.4g} +/- {sd:.4g}, threshold = mean + {n_sd}SD = {threshold:.4g}. "
                f"{len(malignant_clones)} clone(s) with median above this were called malignant"
            )
            labels = pd.Series(
                np.where(clones.isin(malignant_clones), "malignant", "normal"),
                index=cnv_obs.index,
            )
            _report_call(labels, clones, cnv_obs, celltype_key, desc)
            return labels, desc

    if method == "gap":
        if len(clone_median) < 2:
            warn(
                "Only one CNV clone detected; cannot split into malignant/normal. "
                "Abandoning the call (DE will be skipped). Consider lowering "
                "--cells-per-metacell for higher resolution, or revisiting the CNV "
                "window size."
            )
            return pd.Series("unassigned", index=cnv_obs.index), "Cannot call (only one clone)"

        values = clone_median.values.astype(float)
        gaps = np.diff(values)
        i = int(np.argmax(gaps))
        max_gap = float(gaps[i])
        value_range = float(values[-1] - values[0])
        rel = max_gap / value_range if value_range > 0 else 0.0
        log(
            f"Largest gap between clone medians: {max_gap:.4g}"
            f" (range {value_range:.4g}, {rel:.0%})"
            f" / between {clone_median.index[i]} and {clone_median.index[i + 1]}"
        )
        if rel < min_gap_frac:
            warn(
                f"The largest gap is only {rel:.0%} of the range, i.e. there is no "
                f"clear bimodality across clones in cnv_score (threshold {min_gap_frac:.0%}). "
                "Not applying a mechanical malignant/normal split (all 'unassigned', "
                "DE will be skipped). Check the cnv_score distribution and chromosome "
                "heatmap visually, and define groups using biological knowledge "
                "(e.g. --normal-clusters)."
            )
            return (
                pd.Series("unassigned", index=cnv_obs.index),
                f"Cannot call (largest gap is only {rel:.0%} of range)",
            )
        malignant_clones = clone_median.index[i + 1 :].tolist()
        desc = (
            f"gap method: split at the largest gap between clone medians "
            f"({max_gap:.4g}, {rel:.0%} of range). "
            f"{len(malignant_clones)} higher-burden clone(s) called malignant"
        )
        labels = pd.Series(
            np.where(clones.isin(malignant_clones), "malignant", "normal"), index=cnv_obs.index
        )
        _report_call(labels, clones, cnv_obs, celltype_key, desc)
        return labels, desc

    if method == "median":
        warn(
            "The median method splits at the cnv_score median, so it produces a "
            "roughly 50:50 split by construction. It does not detect tumor cells -- "
            "use only for comparison/reproducibility, not for real analysis."
        )
        threshold = float(np.nanmedian(scores))
        labels = pd.Series(
            np.where(scores > threshold, "malignant", "normal"), index=cnv_obs.index
        )
        desc = f"median method (not recommended): cnv_score > {threshold:.4g}"
        _report_call(labels, clones, cnv_obs, celltype_key, desc)
        return labels, desc

    raise ValueError(f"Unknown method: {method}")


def _report_call(
    labels: pd.Series,
    clones: pd.Series,
    cnv_obs: pd.DataFrame,
    celltype_key: str | None,
    desc: str,
) -> None:
    """Print a breakdown and sanity check of the malignancy call."""
    log(desc)
    vc = labels.value_counts().to_dict()
    frac = labels.eq("malignant").mean()
    log(f"Call result: {vc} (malignant fraction {frac:.1%})")
    log("Clone x label:")
    print(pd.crosstab(clones, labels).to_string())
    if celltype_key is not None and celltype_key in cnv_obs:
        if cnv_obs[celltype_key].nunique() > 1:
            log(f"{celltype_key} x label (consistency check):")
            print(pd.crosstab(cnv_obs[celltype_key].astype(str), labels).to_string())
            log(
                "  -> if many immune cells (T/NK, B, Myeloid) fall on the malignant "
                "side, question the call."
            )
    if 0.45 <= frac <= 0.55:
        warn(
            f"malignant fraction is close to {frac:.0%}. This call is per-clone, not "
            "a mechanical median split, but the split could still coincide by chance "
            "-- visually confirm on the chromosome heatmap that the CNV patterns "
            "actually differ."
        )


# ---------------------------------------------------------------------------
# Step 4b. DE analysis branch (pyDESeq2, pseudobulk)
# ---------------------------------------------------------------------------

def run_de_branch(
    mc_adata: ad.AnnData,
    condition_key: str,
    condition_test: str,
    condition_ref: str,
    covariate_key: str | None = "sample_id",
    min_total_counts: int = 10,
    n_cpus: int = 4,
):
    """Pseudobulk DE analysis (pyDESeq2) at metacell level.

    Notes (pseudoreplication):
    - Metacells from the same patient are not independent of each other.
      Passing a sample/patient ID as covariate_key and including it in the
      design partially mitigates this.
    - condition_key is expected to already be merged into mc_adata.obs
      (e.g. the malignant/normal label from the CNV branch).
    - Reports the number of genes remaining after the min_total_counts filter.
    - A covariate with only one level, or fully confounded with condition,
      is dropped from the design (otherwise the design matrix becomes
      singular and dispersion estimation breaks).
    """
    from pydeseq2.ds import DeseqStats

    try:  # absorb pydeseq2 version differences
        from pydeseq2.dds import DeseqDataSet
    except ImportError as exc:  # pragma: no cover
        raise ImportError(f"Failed to import pydeseq2: {exc}") from exc
    try:
        from pydeseq2.default_inference import DefaultInference
    except ImportError:  # older version
        from pydeseq2.dds import DefaultInference  # type: ignore
    if condition_key not in mc_adata.obs:
        raise KeyError(
            f"obs['{condition_key}'] does not exist. Merge the malignant/normal "
            "label from the CNV branch into mc_adata.obs before calling this."
        )

    counts_layer = mc_adata.layers["counts"]
    counts_df = pd.DataFrame(
        _to_dense(counts_layer),
        index=mc_adata.obs_names,
        columns=mc_adata.var_names,
    )
    counts_df = counts_df.round().astype(int)

    # Remove low-count genes
    n_before = counts_df.shape[1]
    keep = counts_df.sum(axis=0) >= min_total_counts
    counts_df = counts_df.loc[:, keep]
    log(
        f"DE: gene count after min_total_counts={min_total_counts} filter"
        f" {counts_df.shape[1]}/{n_before}"
    )
    if counts_df.shape[1] < 500:
        warn(
            f"Only {counts_df.shape[1]} genes remain, which may make pyDESeq2's "
            "dispersion estimation unstable. Consider lowering min_total_counts, or "
            "revisiting metacell size."
        )
    if counts_df.shape[1] == 0:
        raise ValueError("No genes remain after filtering. Lower min_total_counts.")

    metadata = mc_adata.obs.copy()
    metadata[condition_key] = metadata[condition_key].astype(str)

    counts_by_level = metadata[condition_key].value_counts()
    for level in (condition_test, condition_ref):
        if counts_by_level.get(level, 0) < 2:
            raise ValueError(
                f"'{condition_key}' level '{level}' has only "
                f"{counts_by_level.get(level, 0)} metacell(s). DE analysis needs at "
                "least 2 per group (3+ recommended in practice)."
            )
    log(f"DE group sizes: {counts_by_level.to_dict()}")

    design = f"~{condition_key}"
    if covariate_key is not None and covariate_key in metadata.columns:
        metadata[covariate_key] = metadata[covariate_key].astype(str)
        n_levels = metadata[covariate_key].nunique()
        crosstab = pd.crosstab(metadata[covariate_key], metadata[condition_key])
        confounded = ((crosstab > 0).sum(axis=1) == 1).all()
        if n_levels < 2:
            warn(f"Covariate '{covariate_key}' has only one level; dropping it from the design (single sample).")
        elif confounded:
            warn(
                f"Covariate '{covariate_key}' is fully confounded with "
                f"'{condition_key}' (each sample appears in only one group). "
                "Dropping it from the design. Interpret this as a between-patient "
                "comparison with caution."
            )
        else:
            design = f"~{covariate_key} + {condition_key}"
    log(f"DE design: {design}")

    inference = DefaultInference(n_cpus=n_cpus)
    try:
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design=design,
            refit_cooks=True,
            inference=inference,
        )
    except TypeError:  # pydeseq2 < 0.4 takes design_factors instead
        design_factors = [f for f in design.lstrip("~").split("+")]
        design_factors = [f.strip() for f in design_factors]
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design_factors=design_factors,
            refit_cooks=True,
            inference=inference,
        )
    dds.deseq2()

    stat_res = DeseqStats(
        dds,
        contrast=[condition_key, condition_test, condition_ref],
        inference=inference,
    )
    stat_res.summary()

    return dds, stat_res.results_df


# ---------------------------------------------------------------------------
# Main (overall pipeline flow)
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Metacell -> CNV / DE branch pipeline."
        " Passing --scanpy as the first argument delegates to preprocessing"
        " (metacellcnv_scanpy.py); --visualize delegates to interpretation/plotting"
        " (metacellcnv_visualize.py)"
    )
    parser.add_argument(
        "--cellranger-dir",
        action="append",
        default=None,
        help="CellRanger output directory (repeatable). Defaults to CELLRANGER_DIRS if omitted",
    )
    parser.add_argument(
        "--sample-id",
        action="append",
        default=None,
        help="sample_id for each --cellranger-dir (inferred from the path if omitted)",
    )
    parser.add_argument("--gtf", default=GTF_PATH, help="Path to the gene-coordinate annotation GTF")
    parser.add_argument(
        "--gtf-gene-id",
        default=GTF_GENE_ID_ATTR,
        help="GTF gene-name attribute key (GENCODE: gene_name / NCBI-RefSeq: gene)."
        " Use 'auto' to auto-detect",
    )
    parser.add_argument(
        "--chromosome-map",
        default=CHROMOSOME_MAP_PATH,
        help="NCBI assembly report, or a 2-column chromosome conversion table, "
        "converting RefSeq accessions to chr1../chrX/chrY",
    )
    parser.add_argument("--out-dir", default=str(OUTPUT_DIR), help="Output directory for results")
    parser.add_argument(
        "--seacells-fw-iters",
        type=int,
        default=None,
        help="Number of inner Frank-Wolfe iterations (default 50). This is a "
        "regularization parameter, not just a speed knob: increasing it lowers RSS "
        "but can distort the metacell size distribution (measured on one dataset: "
        "largest metacell share grew from 1.8%% at FW=15 to 2.7%% at FW=50 to 33%% "
        "at FW=100). Use 15-50 for more even sizes.",
    )
    parser.add_argument(
        "--lineage", action="store_true",
        help="Build a tumor lineage tree from CNV (Dollo parsimony) and output a "
        "Newick file plus a per-branch event table. Runs a permutation test for "
        "structure before drawing it, and warns if unsupported",
    )
    parser.add_argument(
        "--lineage-level", choices=["chromosome", "window"], default="chromosome",
        help="Granularity of the regions used as lineage characters (default "
        "chromosome; validation data did not support sub-chromosomal structure)",
    )
    parser.add_argument(
        "--lineage-events", choices=["loss", "gain", "both"], default="loss",
        help="Event type used as the character (default loss). Deletions are "
        "irreversible while amplification can occur independently multiple times, "
        "so only loss fits the Dollo parsimony assumption",
    )
    parser.add_argument("--lineage-alpha", type=float, default=0.05,
                        help="Significance level for calling an event (Bonferroni-corrected over regions; default 0.05)")
    parser.add_argument("--lineage-perm", type=int, default=200,
                        help="Number of permutations for the lineage-structure test (default 200)")
    parser.add_argument(
        "--markers",
        default="dog",
        help="Cell-type marker table: a built-in name (dog / human / mouse) or a "
        "CSV path (columns: cell_type, gene, role, use_as_normal_reference, "
        "use_for_labels, validation). Default: dog",
    )
    parser.add_argument(
        "--metacell-init",
        choices=["cssp", "maxmin", "mix"],
        default="cssp",
        help="Archetype initialization: cssp=greedy adaptive column subset "
        "selection (default, deterministic) / maxmin=farthest-point sampling / "
        "mix=combination of both. cssp performed best on RSS and size balance "
        "in testing",
    )
    parser.add_argument(
        "--metacell-seed", type=int, default=0,
        help="Random seed for metacell construction (default 0). This "
        "implementation is deterministic (ARI 0.99 across seeds)",
    )
    parser.add_argument(
        "--cnv-method",
        choices=["infercnv", "bins"],
        default="infercnv",
        help="CNV estimation method. infercnv=windowed running mean (default, "
        "native reimplementation matching infercnvpy) / bins=beta-binomial test "
        "on genomic-bin count composition (reference not floored at 0, includes "
        "an expression-program control; experimental)",
    )
    # --- Flags kept only for backward compatibility with the earlier
    #     SEACells-patch-based implementation; accepted but ignored now that
    #     metacell construction is fully native (so old command lines don't break).
    for _dead, _why in (
        ("--seacells-benchmark", "no SEACells dependency remains to benchmark against"),
        ("--seacells-benchmark-only", "same as above"),
        ("--no-seacells-patch", "replaced by the native implementation (no patch to disable)"),
        ("--no-seacells-patch-verify", "same as above"),
        ("--seacells-float32", "float32 is now an internal setting of metacells_native"),
    ):
        parser.add_argument(_dead, action="store_true",
                            help=f"[removed, ignored] {_why}")
    parser.add_argument(
        "--seacells-line-search", action="store_true",
        help="[removed; if passed, a warning is logged and this is disabled] Exact "
        "line search lowers RSS but distorts the metacell size distribution "
        "(measured: the largest metacell consumed 57.8%% of cells)",
    )
    parser.add_argument(
        "--marker-rule",
        choices=["margin", "relative-z", "absolute"],
        default="margin",
        help="Cell-type classification rule. margin=compare the top two scores "
        "within a cluster (default) / relative-z=legacy compatibility mode "
        "(z-score across clusters; passes only ~1/5 of clusters) / "
        "absolute=fixed threshold via --marker-min-score",
    )
    parser.add_argument(
        "--marker-min-margin",
        type=float,
        default=0.02,
        help="margin mode: minimum gap between the top and 2nd-place scores (default 0.02)",
    )
    parser.add_argument(
        "--marker-min-within-z",
        type=float,
        default=2.0,
        help="margin mode: minimum (top score - mean of rest) / SD of rest (default 2.0)",
    )
    parser.add_argument(
        "--marker-min-score",
        type=float,
        default=None,
        help="Use an absolute score threshold for cell-type classification (implies --marker-rule absolute)",
    )
    parser.add_argument(
        "--marker-min-zscore",
        type=float,
        default=1.0,
        help="Cross-cluster z-score threshold for relative-z mode (default 1.0)",
    )
    parser.add_argument(
        "--normal-celltype",
        default=None,
        help="Add cell type(s) to use as the normal reference (comma-separated, "
        f"e.g. 'Fibroblast'). Default reference types: {', '.join(KNOWN_NORMAL_CELLTYPES)}. "
        f"{'/'.join(MALIGNANCY_CAPABLE_CELLTYPES)} are excluded by default since they "
        "can themselves be the tumor",
    )
    parser.add_argument(
        "--normal-clusters",
        default=None,
        help="Explicitly specify coarse_cluster IDs to use as the normal reference, "
        "comma-separated (e.g. '3,7'). Bypasses marker-score-based classification",
    )
    parser.add_argument(
        "--cnv-reference",
        choices=["celltype", "none"],
        default="celltype",
        help="CNV reference for the malignant call. celltype=normal cell types in "
        "cell_type / none=mean over all metacells",
    )
    parser.add_argument(
        "--mito-chromosome",
        action="append",
        default=None,
        help="Sequence ID(s) of the mitochondrial genome (repeatable). Defaults to "
        f"known IDs {list(MITO_CHROMOSOMES)} plus auto-detection from the GTF.",
    )
    parser.add_argument(
        "--seacell-assignments",
        default=None,
        help="Reuse an existing metacell assignment, skipping optimization. Give "
        "cell_to_metacell.csv or singlecells_qc.h5ad. 'auto' searches --out-dir",
    )
    parser.add_argument(
        "--seacell-min-coverage",
        type=float,
        default=0.80,
        help="Minimum required barcode match rate when reusing an assignment (default 0.80)",
    )
    parser.add_argument(
        "--no-seacell-fill-unmatched",
        action="store_true",
        help="Do not fill in unmatched cells by nearest-neighbor; error out instead",
    )
    parser.add_argument(
        "--no-pctmt-filter",
        action="store_true",
        help="Disable pctMT-based cell exclusion (for samples with implausible "
        "mtDNA composition). The flag itself is still recorded in qc_metrics.csv",
    )
    parser.add_argument(
        "--mito-reference-profile",
        default=None,
        help="[optional input] Reference mtDNA gene composition profile: path to a "
        "mito_gene_profile.csv produced by a run on a different sample (this file "
        "is also an auto-generated output of each run's --out-dir; the check works "
        "with a built-in reference if this is omitted). 'auto' searches --out-dir. "
        "If the given path doesn't exist, warns and continues with the built-in reference",
    )
    parser.add_argument(
        "--keep-mito-in-hvg",
        action="store_true",
        help="Keep mtDNA genes in HVG/PCA (default excludes them). Use only to "
        "restore pre-v2.2 behavior",
    )
    parser.add_argument(
        "--exclude-ribosomal-from-hvg",
        action="store_true",
        help="Also exclude cytoplasmic ribosomal genes (RPL*/RPS*) from HVG. "
        "Consider this if pctMT still correlates with a principal component",
    )
    parser.add_argument(
        "--mito-nmads",
        type=float,
        default=3.5,
        help="Number of MADs used for mitochondrial/nuclear balance anomaly detection (default 3.5)",
    )
    parser.add_argument(
        "--exclude-mito-anomaly",
        action="store_true",
        help="Exclude metacells flagged as mitochondrial anomalies from DE "
        "(for sensitivity analysis when the confounding check shows a skew)",
    )
    parser.add_argument(
        "--expected-doublet-rate",
        type=float,
        default=0.05,
        help="Expected doublet rate for Scrublet (default 0.05); also used as the "
        "fallback selection fraction",
    )
    parser.add_argument(
        "--max-doublet-rate",
        type=float,
        default=0.25,
        help="Treat automatic threshold selection as having failed if the doublet "
        "call rate exceeds this (default 0.25)",
    )
    parser.add_argument(
        "--no-celltype-per-cluster",
        action="store_true",
        help="Do not append the coarse cluster number to cell_type (default appends "
        "it). Appending it splits the normal reference by lineage, enabling the "
        "CNV estimator's bounded mode",
    )
    parser.add_argument(
        "--cnv-refine",
        action="store_true",
        help="Re-estimate CNV using clones that were CNV-flat in the first pass as "
        "multiple normal references (suppresses false positives). This is an "
        "iterative refinement, so always check how much the call changes",
    )
    parser.add_argument(
        "--malignant-call",
        choices=["auto", "reference", "gap", "median"],
        default="auto",
        help="Method for calling malignant metacells. auto=reference if a normal "
        "reference exists, else gap. median is the legacy split (not recommended, "
        "roughly 50:50 by construction)",
    )
    parser.add_argument(
        "--malignant-n-sd",
        type=float,
        default=3.0,
        help="reference method threshold = normal-reference cnv_score mean + N*SD (default 3.0)",
    )
    parser.add_argument(
        "--min-clone-gap",
        type=float,
        default=0.25,
        help="Minimum gap (as a fraction of the clone-median range) for the gap "
        "method to call bimodality (default 0.25)",
    )
    parser.add_argument(
        "--seacells-min-iter", type=int, default=10, help="Minimum number of metacell-optimization iterations (default 10)"
    )
    parser.add_argument(
        "--seacells-max-iter",
        type=int,
        default=50,
        help="Maximum number of metacell-optimization iterations (default 50)",
    )
    parser.add_argument(
        "--drop-bad-metacells",
        action="store_true",
        help="Exclude low-quality metacells (high compactness and low separation)",
    )
    parser.add_argument(
        "--bad-metacell-quantile",
        type=float,
        default=0.95,
        help="Compactness quantile used for the low-quality call (default 0.95)",
    )
    parser.add_argument(
        "--cells-per-metacell", type=int, default=CELLS_PER_METACELL, help="Target cells per metacell"
    )
    parser.add_argument("--n-hvg", type=int, default=N_HVG)
    parser.add_argument("--n-pcs", type=int, default=N_PCS)
    parser.add_argument("--leiden-resolution", type=float, default=0.5)
    parser.add_argument("--metacell-annotation", action="store_true",
                        help="Annotate metacells as typed / measured-mixture / "
                             "unclassified and write metacell_annotation.csv "
                             "(validated panel + karyotype projection + shuffle calibration)")
    parser.add_argument("--no-metacell-annotation", action="store_true",
                        help="Disable --metacell-annotation (already off by default, usually unneeded)")
    parser.add_argument("--annotation-panel-fdr", type=float, default=0.01,
                        help="FDR for the panel-score threshold (against a cell->metacell shuffle null)")
    parser.add_argument("--annotation-mixture-fdr", type=float, default=0.05,
                        help="FDR for the mixture call (binomial test)")
    parser.add_argument("--annotation-perm", type=int, default=30,
                        help="Number of shuffles used to calibrate the panel threshold")
    parser.add_argument("--min-total-counts", type=int, default=10)
    parser.add_argument("--n-cpus", type=int, default=4)
    parser.add_argument(
        "--skip-env-check", action="store_true", help="Do not stop when the environment check fails"
    )
    parser.add_argument(
        "--preflight-only", action="store_true", help="Run only the preflight checks and exit"
    )
    parser.add_argument(
        "--skip-doublet-rescue",
        action="store_true",
        help="Skip the single-cell CNV-based doublet rescue step (reduces compute cost)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    apply_marker_table(getattr(args, "markers", "dog"))

    cellranger_dirs = args.cellranger_dir or CELLRANGER_DIRS
    set_normal_celltypes(
        [c.strip() for c in args.normal_celltype.split(",") if c.strip()]
        if args.normal_celltype
        else None
    )
    normal_clusters = (
        [c.strip() for c in args.normal_clusters.split(",") if c.strip()]
        if args.normal_clusters
        else None
    )
    gtf_gene_id = None if args.gtf_gene_id == "auto" else args.gtf_gene_id

    # --- Step -1: environment check and recording ---
    check_runtime_environment(strict=not args.skip_env_check)
    write_environment_lock(out_dir)

    # Resolve optional input file paths before loading data, so a bad path is reported early
    ref_profile = load_mito_reference_profile(args.mito_reference_profile, out_dir)

    # --- Step -1b: preflight ---
    gtf_info = preflight(
        cellranger_dirs,
        args.gtf,
        gtf_gene_id,
        args.chromosome_map,
        mito_chromosomes=args.mito_chromosome,
    )
    gtf_gene_id = gtf_info["gene_id_attr"]
    exclude_chromosomes = gtf_info["exclude_chromosomes"]
    chromosome_map = gtf_info["chromosome_map"]
    mito_chromosomes = gtf_info["mito_chromosomes"]
    mito_genes = gtf_info["mito_genes"]
    if args.preflight_only:
        log("--preflight-only was given; exiting (no fatal issues detected)")
        return 0

    # --- Step 1: load and preprocess ---
    adata = load_and_preprocess(
        cellranger_dirs,
        sample_ids=args.sample_id,
        n_hvg=args.n_hvg,
        n_pcs=args.n_pcs,
        mito_genes=mito_genes,
        mito_chromosomes=mito_chromosomes,
        exclude_mito_from_hvg=not args.keep_mito_in_hvg,
        exclude_ribo_from_hvg=args.exclude_ribosomal_from_hvg,
    )
    ensure_sample_id(adata)

    # --- Step 1a: coarse clustering ---
    # Used as the unit for QC stratification and doublet detection. Kept
    # independent of the final malignant/normal call by design.
    sc.pp.neighbors(adata, use_rep="X_pca")
    leiden_kwargs = {"resolution": args.leiden_resolution, "key_added": "coarse_cluster"}
    if _leiden_supports_igraph():
        # Recommended implementation on scanpy 1.10+ (flavor='leidenalg' is being phased out)
        leiden_kwargs.update({"flavor": "igraph", "n_iterations": 2, "directed": False})
    try:
        sc.tl.leiden(adata, **leiden_kwargs)
    except (ImportError, TypeError, ValueError) as exc:
        warn(f"Recommended leiden settings are unavailable; retrying with defaults: {exc}")
        sc.tl.leiden(adata, resolution=args.leiden_resolution, key_added="coarse_cluster")
    log(f"Coarse clusters: {adata.obs['coarse_cluster'].nunique()}")

    # --- Step 1a-2: coarse cell-type annotation (secures a CNV normal reference) ---
    # Ordering matters: cell_type must exist before the CNV branch.
    adata.obs["cell_type"] = annotate_coarse_celltype(
        adata,
        cluster_key="coarse_cluster",
        min_score=args.marker_min_score,
        min_zscore=args.marker_min_zscore,
        normal_clusters=normal_clusters,
        rule=args.marker_rule,
        min_margin=args.marker_min_margin,
        min_within_z=args.marker_min_within_z,
        per_cluster_labels=not args.no_celltype_per_cluster,
    )
    print(adata.obs["cell_type"].value_counts())

    # --- Step 1b2: sanity-check mtDNA gene composition (before QC) ---
    mito_comp = check_mito_composition(
        adata, reference_profile=ref_profile, out_path=out_dir / "mito_gene_profile.csv"
    )
    use_pctmt = not args.no_pctmt_filter
    if mito_comp.get("status") == "suspect" and use_pctmt:
        warn(
            "mtDNA composition looks implausible, but the pctMT filter is still "
            "enabled. Strongly consider re-running with --no-pctmt-filter."
        )

    # --- Step 1b: parametric QC (MAD-based, stratified by sample x coarse cluster) ---
    qc = parametric_qc_filter(
        adata,
        sample_key="sample_id",
        coarse_cluster_key="coarse_cluster",
        nmads=5.0,
        use_pctmt=use_pctmt,
    )
    qc.to_csv(out_dir / "qc_metrics.csv")
    # Only low-count / low-gene-count / high-pctMT cells are excluded outright here
    # (high-count cells may be malignant, so they are not excluded at this stage)
    for col in qc.columns:
        adata.obs[col] = qc[col].values
    n_before = adata.n_obs
    adata = adata[~qc["qc_fail"].values].copy()
    log(f"Cells after QC: {adata.n_obs} (excluded {n_before - adata.n_obs})")

    # --- Step 1c: cluster-aware doublet detection ---
    doublet_calls = detect_doublets_cluster_aware(
        adata,
        coarse_cluster_key="coarse_cluster",
        expected_doublet_rate=args.expected_doublet_rate,
        max_doublet_rate=args.max_doublet_rate,
    )
    adata.obs["doublet_score"] = doublet_calls["doublet_score"].values
    adata.obs["predicted_doublet"] = doublet_calls["predicted_doublet"].values

    # --- Step 1d: rescue doublet candidates via single-cell CNV consistency ---
    # Runs a preliminary single-cell-level CNV estimate for this purpose only;
    # the final CNV estimate happens later, per metacell, via run_cnv_branch.
    # cell_type at this stage is a provisional label from coarse marker
    # scoring, so reference_key=None is used (reference = mean over all cells).
    if args.skip_doublet_rescue:
        warn("--skip-doublet-rescue given; skipping Step 1d")
        adata.obs["final_doublet"] = adata.obs["predicted_doublet"].astype(bool)
    else:
        prelim_cnv = run_cnv_branch(
            adata,
            gtf_path=args.gtf,
            gtf_gene_id=gtf_gene_id,
            reference_key=None,
            reference_cat=None,
            window_size=250,
            exclude_chromosomes=exclude_chromosomes,
            chromosome_map=chromosome_map,
            mito_chromosomes=mito_chromosomes,
        )
        rescued_calls = rescue_doublets_by_cnv_consistency(
            doublet_calls=adata.obs[["doublet_score", "predicted_doublet"]],
            cnv_matrix=prelim_cnv.obsm["X_cnv"],
            clone_labels=prelim_cnv.obs["cnv_leiden"],
        )
        adata.obs["final_doublet"] = rescued_calls.values

        n_rescued = int((adata.obs["predicted_doublet"] & ~adata.obs["final_doublet"]).sum())
        print(f"Doublet candidates rescued via CNV consistency: {n_rescued} cells")
        del prelim_cnv

    # Final doublet exclusion (rescued cells are kept)
    n_before = adata.n_obs
    adata = adata[~adata.obs["final_doublet"].values].copy()
    removed = n_before - adata.n_obs
    log(
        f"Cells after doublet exclusion: {adata.n_obs} (excluded {removed}"
        f" = {removed / max(n_before, 1):.1%} of QC-passing cells)"
    )
    if removed / max(n_before, 1) > 0.2:
        warn(
            f"{removed / max(n_before, 1):.0%} of QC-passing cells were removed as "
            "doublets, which is biologically implausible. Check the doublet_score "
            "distribution and --expected-doublet-rate."
        )

    # --- Step 2: build and evaluate metacells ---
    _, _, cells_per_metacell = validate_dimensionality(
        adata.n_obs, adata.n_vars, args.n_hvg, args.n_pcs, args.cells_per_metacell
    )
    assignment_source = args.seacell_assignments
    if assignment_source == "auto":
        found = find_cached_assignments(out_dir)
        if found is None:
            warn("--seacell-assignments auto: no reusable assignment found; building metacells")
        else:
            log(f"--seacell-assignments auto: using {found}")
        assignment_source = str(found) if found else None

    model = None
    if assignment_source:
        load_seacell_assignments(
            adata,
            assignment_source,
            min_coverage=args.seacell_min_coverage,
            fill_unmatched=not args.no_seacell_fill_unmatched,
        )
        log("Skipped metacell optimization (reused an existing assignment)")
    else:
        for _dead in ("no_seacells_patch", "no_seacells_patch_verify",
                      "seacells_float32", "seacells_benchmark",
                      "seacells_benchmark_only"):
            if getattr(args, _dead, False):
                warn(f"--{_dead.replace('_', '-')} has been removed "
                     "(metacell construction no longer depends on SEACells). Ignoring it and continuing")
        model = build_seacells(
            adata,
            cells_per_metacell=cells_per_metacell,
            min_iter=args.seacells_min_iter,
            max_iter=args.seacells_max_iter,
            fw_iters=args.seacells_fw_iters,
            line_search=args.seacells_line_search,
            init=args.metacell_init,
            seed=args.metacell_seed,
        )
    metrics = evaluate_metacells(adata, model, celltype_key="cell_type")
    metrics.to_csv(out_dir / "metacell_metrics.csv", index=False)
    print(metrics.describe())

    # Exclude low-quality metacells (off by default; enable with --drop-bad-metacells)
    bad_mc = find_low_quality_metacells(metrics, compactness_quantile=args.bad_metacell_quantile)
    if args.drop_bad_metacells and len(bad_mc):
        n_before = adata.n_obs
        adata = adata[~adata.obs["SEACell"].isin(bad_mc)].copy()
        log(
            f"Excluded {len(bad_mc)} low-quality metacell(s)"
            f" ({n_before - adata.n_obs} cells)"
        )

    # --- Step 3: aggregate to metacells ---
    mc_adata = aggregate_to_metacells(
        adata, metacell_key="SEACell", sample_key="sample_id", celltype_key="cell_type"
    )

    # --- Step 3a-2: write the cell barcode -> metacell membership table ---
    # This mapping only exists in singlecells_qc.h5ad's obs; the aggregated
    # file retains only the constituent cell count. Since the h5ad can be
    # hundreds of MB, the mapping is also saved as a lightweight standalone
    # CSV (needed downstream to look up a given metacell's constituent cells).
    membership_cols = [
        c
        for c in ("SEACell", "sample_id", "coarse_cluster", "cell_type",
                  "total_counts", "n_genes_by_counts", "pct_counts_mt",
                  "doublet_score")
        if c in adata.obs
    ]
    membership = adata.obs[membership_cols].copy()
    membership.index.name = "barcode"
    membership.to_csv(out_dir / "cell_to_metacell.csv")
    log(
        f"Saved cell->metacell membership table: cell_to_metacell.csv"
        f" ({len(membership):,} cells x {len(membership_cols)} columns)"
    )

    # --- Step 3b: mitochondrial/nuclear expression balance anomaly detection ---
    imbalance = detect_expression_imbalance(mc_adata, sc_adata=adata, nmads=args.mito_nmads)
    for col in imbalance.columns:
        mc_adata.obs[col] = imbalance[col].values
    if len(imbalance.columns):
        imbalance.to_csv(out_dir / "metacell_mito_qc.csv")

    # --- Step 4a: CNV estimation (final, per metacell) ---
    if args.cnv_reference == "none":
        warn(
            "--cnv-reference none: using the mean over all metacells as the CNV "
            "baseline. The malignant signal will be weaker than with a normal "
            "reference (still usable for relative comparison)."
        )
        cnv_reference_key, cnv_reference_cat = None, None
    else:
        cnv_reference_key, cnv_reference_cat = "cell_type", KNOWN_NORMAL_CELLTYPES
    cnv_adata = run_cnv_branch(
        mc_adata,
        gtf_path=args.gtf,
        gtf_gene_id=gtf_gene_id,
        reference_key=cnv_reference_key,
        reference_cat=cnv_reference_cat,
        exclude_chromosomes=exclude_chromosomes,
        chromosome_map=chromosome_map,
        mito_chromosomes=mito_chromosomes,
        method=args.cnv_method,
    )
    # --- Step 4a-2: call malignant metacells (per clone; no median split) ---
    cnv_obs = cnv_adata.obs.reindex(mc_adata.obs_names)
    mc_adata.obs["cnv_score"] = cnv_obs["cnv_score"].values
    mc_adata.obs["cnv_leiden"] = cnv_obs["cnv_leiden"].astype(str).values
    labels, call_desc = call_malignant_metacells(
        mc_adata.obs,
        method=args.malignant_call,
        n_sd=args.malignant_n_sd,
        min_gap_frac=args.min_clone_gap,
    )
    mc_adata.obs["putative_malignant"] = labels.values
    mc_adata.uns["malignant_call"] = call_desc

    def _clone_split(lab_series):
        mal, nor = [], []
        for g in sorted(mc_adata.obs["cnv_leiden"].astype(str).unique(), key=lambda v: (len(v), v)):
            m = (mc_adata.obs["cnv_leiden"].astype(str) == g).values
            top = pd.Series(lab_series[m]).mode()
            (mal if len(top) and top.iat[0] == "malignant" else nor).append(g)
        return mal, nor

    mal_clones, norm_clones = _clone_split(mc_adata.obs["putative_malignant"].values)
    lineage = check_clone_lineage_consistency(
        cnv_adata, mal_clones, norm_clones,
        out_path=out_dir / "clone_chromosome_profiles.csv",
    )

    # --- Step 4a-3: optionally re-estimate CNV using normal clones as multiple references ---
    if args.cnv_refine and len(norm_clones) >= 2 and args.cnv_reference != "none":
        log("=== --cnv-refine: re-estimating CNV using CNV-flat clones as multiple normal references ===")
        log(
            "Clones that were CNV-flat in the first pass are re-supplied as separate "
            "reference categories. With multiple references, the CNV estimator's "
            "bounded mode engages, suppressing false positives from cell-type-"
            "specific expression bias."
        )
        warn(
            "This is an iterative refinement: if the first-pass call was wrong, this "
            "can amplify that error. Always check how much the call changed."
        )
        mc_adata.obs["cnv_reference_group"] = np.where(
            mc_adata.obs["cnv_leiden"].astype(str).isin(norm_clones),
            "ref_" + mc_adata.obs["cnv_leiden"].astype(str),
            "query",
        )
        cnv_adata2 = run_cnv_branch(
            mc_adata,
            gtf_path=args.gtf,
            gtf_gene_id=gtf_gene_id,
            reference_key="cnv_reference_group",
            reference_cat=[f"ref_{g}" for g in norm_clones],
            exclude_chromosomes=exclude_chromosomes,
            chromosome_map=chromosome_map,
            mito_chromosomes=mito_chromosomes,
        )
        o2 = cnv_adata2.obs.reindex(mc_adata.obs_names)
        prev = mc_adata.obs["putative_malignant"].copy()
        mc_adata.obs["cnv_score"] = o2["cnv_score"].values
        mc_adata.obs["cnv_leiden"] = o2["cnv_leiden"].astype(str).values
        labels2, call_desc2 = call_malignant_metacells(
            mc_adata.obs, method=args.malignant_call,
            n_sd=args.malignant_n_sd, min_gap_frac=args.min_clone_gap,
        )
        mc_adata.obs["putative_malignant"] = labels2.values
        changed = int((prev.values != labels2.values).sum())
        log(
            f"Call changes after re-estimation: {changed} / {len(prev)} metacells"
            f" ({changed / max(len(prev), 1):.1%})"
        )
        if changed / max(len(prev), 1) > 0.2:
            warn(
                f"{changed / len(prev):.0%} of calls changed. Check the chromosome "
                "heatmap to determine which pass (first or second) is more plausible."
            )
        call_desc = call_desc2 + " (2nd pass via --cnv-refine)"
        mc_adata.uns["malignant_call"] = call_desc
        cnv_adata = cnv_adata2
        mal_clones, norm_clones = _clone_split(mc_adata.obs["putative_malignant"].values)
        lineage = check_clone_lineage_consistency(
            cnv_adata, mal_clones, norm_clones,
            out_path=out_dir / "clone_chromosome_profiles_refined.csv",
        )
    elif args.cnv_refine:
        warn(
            f"--cnv-refine was given but there are only {len(norm_clones)} normal "
            "clone(s) (2+ required); skipping re-estimation."
        )

    report_mito_confounding(mc_adata.obs)

    # --- Write interpretation caveats to the output ---
    write_interpretation_caveats(
        out_dir,
        mito_comp=mito_comp,
        lineage=lineage,
        use_pctmt=use_pctmt,
        malignant_call=call_desc,
        sample_ids=list(map(str, adata.obs["sample_id"].unique())),
    )
    _write_h5ad(cnv_adata, out_dir / "cnv_metacells.h5ad")

    # --- Step 4a-4: metacell annotation (optional) ---
    # A plain cell_type label (e.g. Myeloid_3) can be driven by genes with
    # little discriminative power. This instead scores against a validated
    # panel with a shuffle-calibrated threshold, checks mixture via a
    # single-cell-level karyotype call, and marks anything undecidable as
    # unclassified with a reason. The existing cell_type column is not overwritten.
    if args.metacell_annotation and not args.no_metacell_annotation:
        try:
            import metacell_annotation as _MA
            log("=== Step 4a-4: metacell annotation ===")

            # Chromosome coordinates only exist on the CNV metacell AnnData
            # (add_genomic_positions writes to cnv_adata.var, and chr-prefix
            # normalization also happens inside run_cnv_branch). The
            # karyotype call needs single-cell counts, so coordinates are
            # copied over to the single-cell var by gene name. Skipping this
            # makes chromosome_karyotype raise, which is caught below and
            # silently skips the whole annotation step.
            _karyo = None
            if "chromosome" in cnv_adata.var.columns:
                _chrom = (cnv_adata.var["chromosome"].astype(object)
                          .reindex(adata.var_names))
                adata.var["chromosome"] = _chrom.values
                _n_pos = int(pd.notna(adata.var["chromosome"]).sum())
                log(f"Copied {_n_pos}/{adata.n_vars} genes with coordinates to the single-cell data")
                _grp = (mc_adata.obs["putative_malignant"].astype(str)
                        .reindex(adata.obs["SEACell"].astype(str)).values)
                _karyo = _MA.chromosome_karyotype(adata, _grp, "malignant", "normal")
                log("Karyotype amplitude (mean |per-chromosome log2 ratio|): "
                    f"{_karyo.abs().mean():.4f}")
            else:
                warn("var['chromosome'] is missing; cannot compute a karyotype. "
                     "Skipping the mixture call and running only the type call")

            _ann = _MA.annotate_metacells(
                adata, "SEACell",
                group_key=None, tumor_groups=None,
                karyotype=_karyo,
                markers=getattr(args, "markers", "dog"),
                panel_fdr=args.annotation_panel_fdr,
                mixture_fdr=args.annotation_mixture_fdr,
                n_perm=args.annotation_perm,
                out_dir=out_dir)
            for _c in ("label", "label_kind", "evidence",
                       "cnv_frac_tumor_cells", "cnv_mixture_q", "karyotype_proj"):
                if _c in _ann.columns:
                    mc_adata.obs["anno_" + _c] = _ann[_c].reindex(
                        mc_adata.obs_names.astype(str)).values
        except Exception as exc:
            warn(f"Metacell annotation failed (continuing without it): {exc}")

    # --- Step 4a-5: tumor lineage from CNV (Dollo parsimony, optional) ---
    # Builds a tree under the constraint that a once-lost chromosomal region
    # cannot be regained. **Structure is tested for first; if unsupported,
    # has_structure=False is set so the tree is not read as a real lineage.**
    if getattr(args, "lineage", False):
        try:
            import cnv_lineage as _LIN
            log("=== Step 4a-5: CNV-based tumor lineage ===")
            _prof = _LIN.profiles_from_adata(cnv_adata, aggregate=args.lineage_level)
            _norm = (mc_adata.obs["putative_malignant"].astype(str)
                     .reindex(_prof.index).values == "normal")
            _lin = _LIN.build_lineage(_prof, _norm, kind=args.lineage_events,
                                      alpha=args.lineage_alpha,
                                      n_perm=args.lineage_perm,
                                      out_dir=out_dir)
            if not _lin["has_structure"]:
                warn("Lineage structure was not supported. cnv_lineage.nwk is still "
                     "written, but should not be interpreted as a real lineage")
        except Exception as exc:
            warn(f"Failed to build the lineage (continuing without it): {exc}")

    mc_adata.obs.to_csv(out_dir / "metacell_obs.csv")
    (out_dir / "malignant_call.txt").write_text(call_desc + "\n", encoding="utf-8")

    # --- Step 4b: DE analysis (malignant vs. normal metacells) ---
    de_input = mc_adata
    if args.exclude_mito_anomaly and "mito_anomaly" in mc_adata.obs:
        n_drop = int(mc_adata.obs["mito_anomaly"].sum())
        if n_drop:
            de_input = mc_adata[~mc_adata.obs["mito_anomaly"].values].copy()
            log(f"--exclude-mito-anomaly: excluded {n_drop} anomalous metacell(s) from DE")
    n_mal = int((de_input.obs["putative_malignant"] == "malignant").sum())
    n_norm = int((de_input.obs["putative_malignant"] == "normal").sum())
    if n_mal < 3 or n_norm < 3:
        warn(
            f"Cannot run DE analysis with {n_mal} malignant / {n_norm} normal "
            "metacells (3+ per group required). CNV results (cnv_metacells.h5ad, "
            "metacell_obs.csv) have already been saved -- inspect the plots, "
            "redefine the groups, and re-run just Step 4b."
        )
    else:
        dds, de_results = run_de_branch(
            de_input,
            condition_key="putative_malignant",
            condition_test="malignant",
            condition_ref="normal",
            covariate_key="sample_id",
            min_total_counts=args.min_total_counts,
            n_cpus=args.n_cpus,
        )
        de_path = out_dir / "de_malignant_vs_normal.csv"
        de_results.sort_values("padj").to_csv(de_path)
        log(f"Saved DE results: {de_path}")
        print(de_results.query("padj < 0.05").sort_values("padj").head(20))

    _write_h5ad(mc_adata, out_dir / "metacells.h5ad")
    _write_h5ad(adata, out_dir / "singlecells_qc.h5ad")
    log("Done. Outputs: " + ", ".join(sorted(p.name for p in out_dir.iterdir())))
    return 0


# ---------------------------------------------------------------------------
# Subcommands: preprocessing-only / visualization-only, via the same entry point
# ---------------------------------------------------------------------------
#   metacellcnv.py             metacell -> CNV -> malignant call -> DE (this file)
#   metacellcnv_scanpy.py      shared preprocessing (also usable standalone)
#   metacellcnv_visualize.py   interpretation report and plots (also usable standalone)
# --scanpy / --visualize, given as the first argument to this script, delegate
# to the corresponding module with the remaining arguments passed through.

def _dispatch(argv: list[str]) -> int | None:
    """Delegate to the corresponding module if --scanpy / --visualize is the first argument."""
    if not argv:
        return None
    head, rest = argv[0], argv[1:]
    if head == "--scanpy":
        import metacellcnv_scanpy as _m
        log("=== metacellcnv --scanpy: running preprocessing ===")
        return int(_m.main(rest) or 0)
    if head == "--visualize":
        import metacellcnv_visualize as _m
        log("=== metacellcnv --visualize: generating interpretation report and plots ===")
        return int(_m.main(rest) or 0)
    return None


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=FutureWarning)
    _rc = _dispatch(sys.argv[1:])
    sys.exit(_rc if _rc is not None else main())
