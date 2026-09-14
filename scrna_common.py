#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scrna_common.py -- Shared building blocks for scRNA-seq preprocessing.

Common ground used by both `metacellcnv_scanpy.py` (unified preprocessing)
and `metacellcnv.py` (CNV/DE analysis). Split out so the preprocessing
steps can be used standalone, independent of CNV analysis; no behavior
changes were made in the split -- functions and constants were moved
as-is from the pipeline.

Contents:
  logging              log / warn
  GTF                  inspect_gtf / parse_gtf_gene_positions / detect_mito_chromosomes
  chromosome names     load_chromosome_map / normalize_chromosome_names / check_exclude_chromosomes
  loading & QC prep    load_and_preprocess / looks_like_raw_counts / _warn_if_stale_cache
  mitochondria         flag_mito_genes / check_mito_composition / load_mito_reference_profile /
                        check_embedding_mito_influence
  QC                   parametric_qc_filter / mad_outlier_mask
  doublets             detect_doublets_cluster_aware
  misc                 ensure_sample_id / validate_dimensionality / _write_h5ad /
                        _leiden_supports_igraph / check_gene_overlap

This module does not depend on CNV estimation (infercnvpy) or DE analysis
(pyDESeq2). See README.md for setup, usage, and version history.
"""
from __future__ import annotations

import gzip
import inspect
import warnings
from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr

__version__ = "1.0"


# ===========================================================================
# Constants
# ===========================================================================

# ---------------------------------------------------------------------------
# Step 0. Configuration (dataset-specific defaults)
# ---------------------------------------------------------------------------

# CellRanger output directory. List multiple samples as a list.
CELLRANGER_DIRS: list[str] = ["filtered_feature_bc_matrix/"]

GTF_PATH = "genes.gtf"

# Attribute key that holds the gene name in the GTF attribute column.
# "gene_name" for GENCODE-format files (the standard human/mouse CellRanger
# reference). "gene" for NCBI/RefSeq-format files (e.g. the dog genome
# Dog10K_Boxer_Tasha used in this project).
# The pipeline auto-detects this from the actual data during preflight and
# stops if it doesn't match.
GTF_GENE_ID_ATTR = "gene"

# Path to an NCBI assembly report (or a 2-column mapping table). When set,
# RefSeq accessions are converted to chr1..chrN / chrX / chrY, which makes
# EXCLUDE_CHROMOSOMES=("chrX","chrY") work as intended. For the dog
# Dog10K_Boxer_Tasha reference this is
# GCF_000002285.5_Dog10K_Boxer_Tasha_assembly_report.txt.
CHROMOSOME_MAP_PATH: str | None = None

CELLS_PER_METACELL = 75  # recommended value from the official SEACells tutorial

N_HVG = 1500

N_PCS = 50

# Sequence IDs for the mitochondrial genome.
# Determined from the GTF's first column (chromosome/sequence ID) rather
# than from gene-name patterns, which is more reliable.
# For dog (Canis lupus familiaris, CanFam6 / Dog10K_Boxer_Tasha), mtDNA is
# NC_002008.4. GENCODE/UCSC-style references use chrM, Ensembl-style use MT.
# Can be extended/overridden with --mito-chromosome, and is also
# auto-detected from the GTF (see detect_mito_chromosomes: a sequence
# where all genes fit within 25 kb is treated as mtDNA).
MITO_CHROMOSOMES: tuple[str, ...] = ("NC_002008.4", "chrM", "chrMT", "MT", "M")

# Gene-name fallback, used only when the above cannot identify mtDNA.
MITO_PREFIXES = ("MT-", "MT.", "MT_")

MITO_SYMBOLS = {
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB",
}

# Upper bound on sequence length for calling a contig mtDNA (mammalian
# mtDNA is roughly 16.6 kb).
MITO_MAX_LENGTH = 25_000

MITO_MIN_GENES = 5

# Nuclear-encoded mitochondrial machinery genes (matched by prefix and by
# exact symbol). OXPHOS complexes combine mtDNA-encoded and nuclear-encoded
# subunits in a fixed stoichiometry, so this ratio should stay roughly
# constant across metacells. A deviation is a candidate signal for
# mitochondrial dysfunction, altered mtDNA copy number, dead-cell
# contamination, or ambient RNA.
# Note: COX1/COX2/COX3 etc. are mtDNA-encoded, so the actual check always
# excludes genes already flagged in var['mt'] from this set.
NUCLEAR_MITO_PREFIXES = (
    "NDUF",    # nuclear-encoded subunits of complex I
    "SDH",     # complex II (entirely nuclear-encoded)
    "UQCR",    # complex III
    "COX",     # nuclear-encoded subunits/assembly factors of complex IV (COX4I1, COX5A, COX10, etc.)
    "ATP5",    # nuclear-encoded subunits of ATP synthase
    "MRPL",    # mitochondrial ribosome, large subunit
    "MRPS",    # mitochondrial ribosome, small subunit
    "TIMM",    # inner-membrane import machinery
    "TOMM",    # outer-membrane import machinery
)

NUCLEAR_MITO_GENES = ("CYC1", "CYCS", "POLRMT", "TFAM", "TFB2M", "SSBP1", "OPA1", "MFN2")

# Prefixes for cytoplasmic ribosomal genes (used when excluding them from HVGs).
# MRPL*/MRPS* (mitochondrial ribosome, nuclear-encoded) start with M so they
# are not matched by these.
RIBO_PREFIXES = ("RPL", "RPS")

# Warn that "the embedding is being dragged by mitochondrial content" when
# the correlation between pctMT and a principal component exceeds this.
MITO_PC_CORR_ALERT = 0.30

# Warning threshold for large datasets.
MAX_CELLS_WARN = 100_000

OUTPUT_DIR = Path("results")

# ---------------------------------------------------------------------------
# Step -1b. GTF inspection
# ---------------------------------------------------------------------------

CANDIDATE_GENE_ID_KEYS = ["gene_name", "gene", "Name", "gene_symbol", "gene_id"]

# ---------------------------------------------------------------------------
# Step 1-b2. mtDNA gene composition sanity check (detects e.g. rescue alignments)
# ---------------------------------------------------------------------------

# mtDNA 13-gene composition (%) measured on a normally-aligned dog sample
# (Case1, Dog10K_Boxer_Tasha). Because mtDNA is transcribed as a single
# polycistronic transcript, this composition is expected to be fairly
# stable across samples. ND6, transcribed from the light strand, is known
# to be the lowest.
MITO_REFERENCE_PROFILE: dict[str, float] = {
    "ND1": 2.3, "ND2": 2.2, "ND3": 0.7, "ND4": 11.7, "ND4L": 1.0, "ND5": 0.8,
    "ND6": 0.3, "COX1": 10.4, "COX2": 15.6, "COX3": 27.1, "ATP6": 24.6,
    "ATP8": 0.3, "CYTB": 3.1,
}

# Flag as abnormal if ND6 ranks within this many places from the top (it is
# normally near the bottom).
MITO_ND6_RANK_ALERT = 3

# Genes with a composition share below this value are considered "near
# undetected".
MITO_NEAR_ZERO_PCT = 0.5

# Flag as abnormal if more than this many genes are "near undetected".
MITO_MAX_NEAR_ZERO = 4

# Flag as abnormal if the rank correlation with the reference profile falls
# below this.
MITO_MIN_SPEARMAN = 0.3


# ===========================================================================
# Functions
# ===========================================================================


# ---------------------------------------------------------------------------
# Common utilities
# ---------------------------------------------------------------------------

def log(msg: str) -> None:
    print(f"[pipeline] {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"[pipeline][WARN] {msg}", flush=True)


def _open_maybe_gzip(path: Path | str):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _to_dense(matrix):
    """Return an ndarray whether the input is sparse or dense."""
    if hasattr(matrix, "toarray"):
        return matrix.toarray()
    return np.asarray(matrix)


def looks_like_raw_counts(matrix, n_sample: int = 200_000) -> bool:
    """Check whether `.X` looks like raw counts (non-negative integers)."""
    import scipy.sparse as sp

    if sp.issparse(matrix):
        values = matrix.data
    else:
        values = np.asarray(matrix).ravel()
    if values.size == 0:
        return True
    if values.size > n_sample:
        rng = np.random.default_rng(0)
        values = rng.choice(values, size=n_sample, replace=False)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return True
    if (values < 0).any():
        return False
    return bool(np.allclose(values, np.round(values)))


def inspect_gtf(
    gtf_path: str,
    gene_id_attr: str | None = None,
    sample_limit: int = 5000,
    extra_mito_chromosomes: tuple[str, ...] | list[str] | None = None,
) -> dict:
    """Inspect a GTF's format, attribute keys, and chromosome naming, and
    return its gene-name set.

    Returns a dict with:
        gene_id_attr   : the attribute key to actually use
        gene_names     : the set of gene names extracted using that key
        chromosomes    : the set of chromosome names
        has_chr_prefix : whether names use a 'chr' prefix
    """
    log("=== GTF check ===")
    path = Path(gtf_path)
    if not path.exists():
        raise FileNotFoundError(
            f"GTF file does not exist: {gtf_path}\n"
            "Per-gene chromosome/start/end coordinates are required for CNV estimation."
        )
    log(f"GTF: {path}")

    from collections import Counter

    key_counter: Counter[str] = Counter()
    gff3_like = 0
    chromosomes: set[str] = set()
    n_gene_rows = 0

    with _open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chromosomes.add(fields[0])
            attr_str = fields[8].strip()
            if "=" in attr_str and '"' not in attr_str:
                gff3_like += 1
            for kv in attr_str.split(";"):
                kv = kv.strip()
                if not kv or " " not in kv:
                    continue
                key, _ = kv.split(" ", 1)
                key_counter[key] += 1
            n_gene_rows += 1
            if n_gene_rows >= sample_limit:
                break

    if n_gene_rows == 0:
        raise ValueError(
            "Could not find any 'gene' rows in the GTF. Confirm the file is"
            " tab-separated with feature type in column 3"
            " (a GFF3 file needs to be converted to GTF first)."
        )
    if gff3_like > n_gene_rows * 0.5:
        raise ValueError(
            "The attribute column looks like `key=value` format, which means"
            " this is GFF3 rather than GTF. Convert it to GTF (e.g. with"
            " gffread) before using it."
        )

    log(f"Sampled {n_gene_rows} gene rows. Top attribute keys: {key_counter.most_common(8)}")

    found = [k for k in CANDIDATE_GENE_ID_KEYS if key_counter.get(k, 0) > 0]
    if not found:
        raise ValueError(
            f"None of the expected gene-name attribute keys {CANDIDATE_GENE_ID_KEYS}"
            f" were found. Actual keys: {sorted(key_counter)}"
        )
    log(f"Candidate gene-name attribute keys: {found}")

    if gene_id_attr is None:
        gene_id_attr = found[0]
        log(f"Auto-detected GTF_GENE_ID_ATTR: '{gene_id_attr}'")
    elif key_counter.get(gene_id_attr, 0) == 0:
        raise ValueError(
            f"The configured GTF_GENE_ID_ATTR='{gene_id_attr}' does not exist"
            f" in this GTF. Candidates: {found} (GENCODE uses 'gene_name',"
            " NCBI/RefSeq uses 'gene')"
        )
    else:
        log(f"GTF_GENE_ID_ATTR='{gene_id_attr}' exists in the GTF")

    if "gene_name" not in key_counter and "gene" in key_counter:
        warn(
            "No GENCODE-standard 'gene_name' attribute, but 'gene' is present."
            " This is the typical pattern for NCBI/RefSeq-format GTFs"
            " (e.g. the dog Dog10K_Boxer_Tasha reference)."
        )

    has_chr_prefix = any(c.startswith("chr") for c in chromosomes)
    sample_chroms = sorted(chromosomes)[:5]
    if has_chr_prefix:
        log(f"Chromosome naming: 'chr'-prefixed (e.g. {sample_chroms})")
    else:
        warn(
            f"Chromosome names are not 'chr'-prefixed (e.g. {sample_chroms})."
            " This looks like RefSeq accession format, which infercnvpy"
            " cannot recognize at all as-is. This script works around it via"
            " normalize_chromosome_names() by prepending 'chr', but if you"
            " want to exclude X/Y, pass an NCBI assembly report via"
            " --chromosome-map to convert to chr1../chrX/chrY."
        )

    gene_pos = parse_gtf_gene_positions(str(path), gene_id_type=gene_id_attr)
    gene_names = set(gene_pos.index.astype(str))
    log(f"Extracted {len(gene_names)} gene names from the GTF (attribute '{gene_id_attr}')")

    mito_chroms = detect_mito_chromosomes(gene_pos, extra=extra_mito_chromosomes)
    mito_genes = set(
        gene_pos.index[gene_pos["chromosome"].astype(str).isin(mito_chroms)].astype(str)
    )

    return {
        "gene_id_attr": gene_id_attr,
        "gene_names": gene_names,
        "gene_positions": gene_pos,
        "chromosomes": chromosomes,
        "has_chr_prefix": has_chr_prefix,
        "mito_chromosomes": mito_chroms,
        "mito_genes": mito_genes,
    }


def detect_mito_chromosomes(
    gene_pos: pd.DataFrame,
    extra: tuple[str, ...] | list[str] | None = None,
    max_length: int = MITO_MAX_LENGTH,
    min_genes: int = MITO_MIN_GENES,
) -> list[str]:
    """Identify the mitochondrial genome's sequence ID(s) from GTF gene coordinates.

    Relying on gene-name patterns (an "MT-" prefix, ND1/COX1, etc.) is
    fragile since naming varies across species/annotations. The GTF has
    sequence IDs and coordinates, which is a more direct basis for the
    decision. Detection is two-stage:

    1. Match against known IDs (MITO_CHROMOSOMES plus anything added via
       --mito-chromosome). For dog CanFam6 / Dog10K_Boxer_Tasha this is
       NC_002008.4.
    2. Inference from sequence length. Mammalian mtDNA is roughly 16.6 kb,
       so a sequence is called mtDNA if all of its genes fit within
       max_length (default 25 kb) and it has at least min_genes genes
       (mtDNA carries 13 protein-coding genes plus rRNA/tRNA, so its gene
       density distinguishes it from an unplaced scaffold of similar length).
    """
    known = {str(c) for c in MITO_CHROMOSOMES}
    if extra:
        known |= {str(c) for c in extra}

    found: list[str] = []
    reasons: dict[str, str] = {}
    for chrom, sub in gene_pos.groupby(gene_pos["chromosome"].astype(str)):
        name = str(chrom)
        if name in known or name[3:].upper() in {"M", "MT"} and name.lower().startswith("chr"):
            found.append(name)
            reasons[name] = "matches a known mtDNA ID"
            continue
        if len(sub) >= min_genes and int(sub["end"].max()) <= max_length:
            found.append(name)
            reasons[name] = (
                f"inferred from sequence length ({len(sub)} genes within {int(sub['end'].max()):,} bp)"
            )

    if found:
        for name in found:
            log(f"Identified as the mitochondrial genome: {name} ({reasons[name]}, "
                f"{int((gene_pos['chromosome'].astype(str) == name).sum())} genes)")
    else:
        warn(
            "Could not identify a mitochondrial genome sequence from the GTF."
            f" No match against known IDs {sorted(known)}, and no sequence"
            f" fits within {max_length:,} bp either. Specify one explicitly"
            " with --mito-chromosome (NC_002008.4 for dog CanFam6). Without"
            " it, detection falls back to gene-name matching."
        )
    return sorted(found)


def check_gene_overlap(
    var_names, gtf_gene_names: set[str], fail_below: float = 0.2, warn_below: float = 0.8
) -> float:
    """Report the overlap rate between var_names and the GTF gene-name set."""
    var_set = set(map(str, var_names))
    if not var_set or not gtf_gene_names:
        warn("Skipping overlap check (var_names or the GTF gene set is empty)")
        return 0.0
    overlap = var_set & gtf_gene_names
    rate = len(overlap) / len(var_set)
    msg = f"{len(overlap)} of {len(var_set)} var_names match the GTF ({rate:.1%})"
    if rate < fail_below:
        raise ValueError(
            msg + f" -- overlap is below {fail_below:.0%}. Most genes would be"
            " dropped from CNV estimation. Check GTF_GENE_ID_ATTR, and the"
            " format of var_names (symbol / LOC ID / Ensembl-style ID)."
        )
    if rate < warn_below:
        warn(msg + " -- overlap is somewhat low. Double-check gtf_gene_id and the naming convention.")
    else:
        log(msg)
    return rate


def load_chromosome_map(path: str | None) -> dict[str, str]:
    """Load a chromosome-name mapping table.

    Accepted formats:
    1. An NCBI assembly report (`*_assembly_report.txt`). Maps
       RefSeq-Accn -> 'chr'+Assigned-Molecule (assembled-molecule rows
       only). Mitochondria is mapped to 'chrM'.
    2. A two-column tab- or comma-separated text file
       (source-name<TAB>target-name).

    This is the practical fix for cases where a species' X/Y RefSeq
    accessions aren't otherwise identifiable: supplying the assembly
    report recovers chr1..chrN, chrX, chrY so that
    EXCLUDE_CHROMOSOMES=("chrX","chrY") works as intended.
    """
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Chromosome map file does not exist: {path}")

    mapping: dict[str, str] = {}
    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    is_assembly_report = any(l.startswith("# Sequence-Name") for l in lines)

    if is_assembly_report:
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 7:
                continue
            seq_role, assigned_molecule = fields[1].strip(), fields[2].strip()
            refseq_accn = fields[6].strip()
            if seq_role != "assembled-molecule" or refseq_accn in ("", "na"):
                continue
            name = "chrM" if assigned_molecule.upper() in ("MT", "M") else f"chr{assigned_molecule}"
            mapping[refseq_accn] = name
            if len(fields) > 4 and fields[4].strip() not in ("", "na"):
                mapping.setdefault(fields[4].strip(), name)  # map GenBank-Accn to the same name too
        log(f"Loaded {len(mapping)} chromosome mappings from the assembly report: {path}")
    else:
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.replace(",", "\t").split("\t")
            if len(fields) < 2:
                continue
            mapping[fields[0].strip()] = fields[1].strip()
        log(f"Loaded {len(mapping)} chromosome mappings from the map file: {path}")
    return mapping


def normalize_chromosome_names(
    var: pd.DataFrame,
    chromosome_map: dict[str, str] | None = None,
    mito_chromosomes: list[str] | None = None,
) -> pd.DataFrame:
    """Normalize var['chromosome'] to the 'chr'-prefixed form infercnvpy expects.

    Important implementation detail:
    infercnvpy only processes chromosomes satisfying `x.startswith("chr")`
    internally (`_running_mean_by_chromosome`). Left in RefSeq-accession
    form (e.g. 'NC_006583.4'), zero chromosomes would match and CNV
    estimation fails with the confusing
    `ValueError: not enough values to unpack (expected 2, got 0)`.
    So normalization here is mandatory, not optional -- without it,
    CNV estimation itself fails outright (it's not just that
    exclude_chromosomes silently becomes a no-op).

    - If a chromosome_map (e.g. from an assembly report) is given, it is
      applied first.
    - Any remaining non-'chr' names get 'chr' prepended (per-chromosome
      smoothing is independent, so biological chromosome order doesn't
      affect correctness, but it does affect the natural sort order used
      when plotting).
    """
    chromosome_map = dict(chromosome_map or {})
    # Map mtDNA to "chrM". infercnvpy excludes chrM internally, so this
    # alone keeps mitochondrial genes out of the CNV computation (mtDNA
    # copy number is unrelated to nuclear-genome copy number and its
    # expression is highly variable, so it should always be excluded).
    for name in mito_chromosomes or []:
        chromosome_map.setdefault(str(name), "chrM")
    chrom = var["chromosome"].astype(object)
    mapped = chrom.map(lambda c: chromosome_map.get(c, c) if isinstance(c, str) else c)

    n_mapped = int(sum(1 for a, b in zip(chrom, mapped) if isinstance(a, str) and a != b))
    if n_mapped:
        log(f"Applied the chromosome map: {n_mapped} genes")

    def _prefix(c):
        if not isinstance(c, str):
            return c
        return c if c.startswith("chr") else f"chr{c}"

    prefixed = mapped.map(_prefix)
    n_prefixed = int(sum(1 for a, b in zip(mapped, prefixed) if isinstance(a, str) and a != b))
    if n_prefixed:
        warn(
            f"Prepended 'chr' to the chromosome name of {n_prefixed} genes."
            " infercnvpy only processes chromosomes starting with 'chr', so"
            " CNV estimation would otherwise fail on RefSeq-accession-style"
            " names. To exclude X/Y, pass an NCBI assembly report via"
            " --chromosome-map to convert them to chrX/chrY."
        )
    var = var.copy()
    var["chromosome"] = prefixed
    return var


def check_exclude_chromosomes(
    exclude: tuple[str, ...], chromosomes: set[str]
) -> tuple[str, ...]:
    """Verify that EXCLUDE_CHROMOSOMES matches the chromosome names in the actual data."""
    if not exclude:
        log("EXCLUDE_CHROMOSOMES is empty (no sex chromosomes will be excluded)")
        return exclude
    matched = [c for c in exclude if c in chromosomes]
    missing = [c for c in exclude if c not in chromosomes]
    if missing:
        warn(
            f"EXCLUDE_CHROMOSOMES entries {missing} do not exist in the GTF"
            " chromosome names. infercnvpy silently ignores unmatched names"
            " rather than erroring, so these are effectively not excluded."
            f" Example chromosome names in the data: {sorted(chromosomes)[:5]}"
        )
    if matched:
        log(f"Chromosomes matched for exclusion: {matched}")
    return tuple(matched)


def validate_dimensionality(
    n_obs: int, n_vars: int, n_hvg: int, n_pcs: int, cells_per_metacell: int
) -> tuple[int, int, int]:
    """Adjust HVG count / PCA dimensionality / metacell size to fit the
    number of cells and genes."""
    log("=== Dimensionality/scale check ===")
    if n_obs > MAX_CELLS_WARN:
        warn(
            f"Cell count {n_obs:,} exceeds {MAX_CELLS_WARN:,}. SEACells/MC2 are"
            " reported to run out of memory or take a long time at this scale."
            " Consider splitting samples or a more scalable method such as CAMP."
        )

    new_n_pcs = int(min(n_pcs, max(2, min(n_obs, n_vars) - 1)))
    if new_n_pcs != n_pcs:
        warn(f"Reduced N_PCS from {n_pcs} to {new_n_pcs} (too large for {n_obs} cells)")
    new_n_hvg = int(min(n_hvg, n_vars))
    if new_n_hvg <= new_n_pcs:
        new_n_hvg = int(min(n_vars, new_n_pcs + 10))
    if new_n_hvg != n_hvg:
        warn(f"Adjusted N_HVG from {n_hvg} to {new_n_hvg} ({n_vars} genes available)")

    new_cpm = cells_per_metacell
    if n_obs < cells_per_metacell * 3:
        new_cpm = max(5, n_obs // 5)
        warn(
            f"CELLS_PER_METACELL={cells_per_metacell} is too large for"
            f" {n_obs} cells. Reduced to {new_cpm} (at least 3 metacells are required)."
        )
    log(f"Using: N_HVG={new_n_hvg}, N_PCS={new_n_pcs}, CELLS_PER_METACELL={new_cpm}")
    return new_n_hvg, new_n_pcs, new_cpm


# ---------------------------------------------------------------------------
# Step 1. Loading, post-QC processing, dimensionality reduction
# ---------------------------------------------------------------------------

def flag_mito_genes(
    adata: ad.AnnData,
    mito_genes: set[str] | None = None,
    mito_chromosomes: list[str] | None = None,
) -> int:
    """Flag mitochondrial genes in adata.var['mt'] and return the count.

    Prefers the GTF-derived gene set (mito_genes) -- genes located on the
    mtDNA sequence ID(s) identified by detect_mito_chromosomes(), e.g.
    NC_002008.4 for dog CanFam6. Because this is based on sequence ID, it
    is unaffected by gene-naming conventions (presence of an "MT-" prefix,
    LOC-style IDs, etc.).

    Falls back to gene-name pattern matching (an "MT-" prefix, or known
    symbols like ND1/COX1) only when mito_genes is not supplied (i.e. the
    GTF is unavailable).
    """
    if mito_genes:
        mask = np.isin(np.asarray(adata.var_names, dtype=object), sorted(mito_genes))
        source = (
            f"GTF sequence ID {mito_chromosomes}" if mito_chromosomes else "GTF-derived gene set"
        )
    else:
        upper = np.asarray(adata.var_names.str.upper(), dtype=object)
        mask = np.zeros(adata.n_vars, dtype=bool)
        for prefix in MITO_PREFIXES:
            mask |= np.array([str(name).startswith(prefix) for name in upper], dtype=bool)
        mask |= np.isin(upper, sorted(MITO_SYMBOLS))
        source = "gene-name pattern (fallback)"

    adata.var["mt"] = mask
    n_mito = int(mask.sum())
    if n_mito == 0:
        warn(
            f"Could not identify any mitochondrial genes (method: {source})."
            " pctMT-based QC is effectively disabled. Pass the mtDNA"
            " sequence ID via --mito-chromosome (NC_002008.4 for dog CanFam6)."
        )
    else:
        log(f"Identified {n_mito} mitochondrial genes (method: {source})")
        if mito_genes:
            missing = len(mito_genes) - n_mito
            if missing > 0:
                log(
                    f"  {missing} of {len(mito_genes)} mtDNA genes in the GTF are"
                    " not present in var_names (e.g. excluded from the CellRanger reference)"
                )
    return n_mito


def _warn_if_stale_cache(path: str | Path, cache_dir: str | Path = "cache") -> None:
    """Guard against scanpy's read_10x_mtx(cache=True) silently reading a stale cache.

    scanpy keys its cache purely by path, so re-running CellRanger at the
    same path leaves the old cache in place even though the matrix content
    changed. That swap would otherwise go unnoticed and invalidate the
    whole analysis, so warn whenever the cache is older than the matrix.
    """
    src = Path(path)
    cache = Path(cache_dir)
    if not cache.exists():
        return
    key = str(src.resolve()).lstrip("/").replace("/", "-")
    hits = [c for c in cache.glob("*.h5ad") if c.stem.startswith(key)]
    if not hits:
        return
    mtx = [m for m in (src / "matrix.mtx.gz", src / "matrix.mtx") if m.exists()]
    if not mtx:
        return
    newest_src = max(m.stat().st_mtime for m in mtx)
    for c in hits:
        if c.stat().st_mtime < newest_src:
            warn(
                f"scanpy cache {c} is older than the matrix file."
                " If CellRanger was re-run at the same path, the stale"
                f" matrix would be loaded. If in doubt, delete {cache}/ and re-run."
            )


def load_and_preprocess(
    paths: str | list[str],
    sample_ids: list[str] | None = None,
    n_hvg: int = N_HVG,
    n_pcs: int = N_PCS,
    mito_genes: set[str] | None = None,
    mito_chromosomes: list[str] | None = None,
    exclude_mito_from_hvg: bool = True,
    exclude_ribo_from_hvg: bool = False,
) -> ad.AnnData:
    """Load one or more CellRanger outputs and run normalization, HVG selection, and PCA.

    - Multiple samples are concatenated, and obs['sample_id'] is set.
    - obs['sample_id'] is set explicitly even for a single sample.
    - Raw counts are kept in adata.layers['counts'] (required for metacell
      aggregation, CNV, and DE).
    - Checks that `.X` is not already normalized.
    """
    if isinstance(paths, (str, Path)):
        paths = [str(paths)]
    if sample_ids is not None and len(sample_ids) != len(paths):
        raise ValueError("Length of sample_ids does not match the number of CellRanger output directories")

    adatas: list[ad.AnnData] = []
    used_ids: list[str] = []
    for i, path in enumerate(paths):
        _warn_if_stale_cache(path)
        sub = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
        sub.var_names_make_unique()
        if sample_ids is not None:
            sid = sample_ids[i]
        else:
            sid = _infer_sample_id(path, i)
        sub.obs["sample_id"] = sid
        used_ids.append(sid)
        log(f"Loaded: {path} -> sample_id='{sid}' ({sub.n_obs} cells x {sub.n_vars} genes)")
        adatas.append(sub)

    if len(adatas) == 1:
        adata = adatas[0]
    else:
        n_vars = {a.n_vars for a in adatas}
        if len(n_vars) > 1:
            warn(f"Gene count differs across samples {n_vars}. Concatenating on the common genes (inner join).")
        adata = ad.concat(
            adatas, join="inner", label=None, keys=used_ids, index_unique="-", merge="first"
        )
        log(f"After concat: {adata.n_obs} cells x {adata.n_vars} genes / sample_id: {used_ids}")

    adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")

    if not looks_like_raw_counts(adata.X):
        raise ValueError(
            "Loaded .X does not look like non-negative integers. Confirm"
            " that this points at CellRanger's filtered_feature_bc_matrix"
            " (raw counts). Applying normalize_total+log1p again to an"
            " already-normalized/log-transformed matrix would give incorrect values."
        )

    # Keep raw counts in a separate layer (used for aggregation, CNV estimation, and DE)
    adata.layers["counts"] = adata.X.copy()

    flag_mito_genes(adata, mito_genes=mito_genes, mito_chromosomes=mito_chromosomes)

    n_hvg, n_pcs, _ = validate_dimensionality(
        adata.n_obs, adata.n_vars, n_hvg, n_pcs, CELLS_PER_METACELL
    )

    # Standard normalization (used for visualization, PCA, and the SEACells kernel)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    hvg_kwargs = {"n_top_genes": n_hvg}
    if adata.obs["sample_id"].nunique() > 1:
        hvg_kwargs["batch_key"] = "sample_id"
    sc.pp.highly_variable_genes(adata, **hvg_kwargs)

    # --- Exclude mtDNA (and optionally ribosomal) genes from the HVG set used for PCA ---
    #
    # Why this matters: PCA is run on log1p values without sc.pp.scale, so
    # high-variance genes dominate the principal components. Even though
    # mtDNA genes are a small fraction of total counts, their expression is
    # high and highly variable across cells, so they tend to end up in the
    # HVG set. Leaving them in can noticeably reshape the nearest-neighbor
    # graph and correlate strongly with a principal component. Since the
    # neighbor structure underlies QC stratification, doublet-detection
    # batching, and the SEACells kernel, that effect propagates through the
    # whole downstream pipeline.
    #
    # Mitochondrial content reflects a cell-state axis (stress, death) that
    # is usually not a useful axis for separating cell types. This is worse
    # for rescue-aligned samples, where neighbors can end up driven by the
    # purely technical property of "how many reads survived rescue".
    excluded: list[str] = []
    hv_mask = np.asarray(adata.var["highly_variable"], dtype=bool)
    if exclude_mito_from_hvg and "mt" in adata.var:
        is_mt = np.asarray(adata.var["mt"], dtype=bool)
        drop = hv_mask & is_mt
        if drop.any():
            excluded += list(adata.var_names[drop])
            hv_mask &= ~is_mt
    if exclude_ribo_from_hvg:
        upper = np.asarray(adata.var_names.str.upper(), dtype=object)
        is_ribo = np.zeros(adata.n_vars, dtype=bool)
        for prefix in RIBO_PREFIXES:
            is_ribo |= np.array([str(n).startswith(prefix) for n in upper], dtype=bool)
        drop = hv_mask & is_ribo
        if drop.any():
            excluded += list(adata.var_names[drop])
            hv_mask &= ~is_ribo
    if excluded:
        log(
            f"Excluded {len(excluded)} genes from the HVG set (to keep the"
            f" embedding from being dominated by cell state): {sorted(excluded)[:12]}"
            f"{' ...' if len(excluded) > 12 else ''}"
        )
        log(f"  HVGs used for PCA: {int(hv_mask.sum())} genes")
    adata.var["hvg_for_pca"] = hv_mask

    try:  # scanpy >= 1.10 uses mask_var, earlier versions use use_highly_variable
        sc.tl.pca(adata, n_comps=n_pcs, mask_var="hvg_for_pca")
    except TypeError:  # pragma: no cover - older scanpy
        adata.var["highly_variable"] = hv_mask
        sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)

    check_embedding_mito_influence(adata)
    return adata


def check_embedding_mito_influence(
    adata: ad.AnnData, threshold: float = MITO_PC_CORR_ALERT
) -> float:
    """Measure the correlation between pctMT and the principal components to
    check whether the embedding is being dragged by mitochondrial content.

    Even after excluding mtDNA genes from the HVG set, some correlation can
    remain via nuclear-encoded stress-response genes etc. Recording the
    residual magnitude makes it possible to later check what is driving the
    neighbor structure. Also reports the partial correlation with total UMI
    (depth) regressed out.
    """
    if "X_pca" not in adata.obsm or "pct_counts_mt" not in adata.obs:
        return float("nan")
    P = np.asarray(adata.obsm["X_pca"])
    y = np.asarray(adata.obs["pct_counts_mt"], dtype=float)
    if not np.isfinite(y).any() or np.nanstd(y) == 0:
        return float("nan")
    depth = np.log1p(np.asarray(adata.obs["total_counts"], dtype=float))
    design = np.c_[np.ones_like(depth), depth]

    def _partial(x: np.ndarray) -> float:
        rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
        ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
        if np.std(rx) == 0 or np.std(ry) == 0:
            return 0.0
        return float(np.corrcoef(rx, ry)[0, 1])

    simple = np.array([abs(float(np.corrcoef(P[:, i], y)[0, 1])) for i in range(P.shape[1])])
    worst = int(np.argmax(simple))
    part = abs(_partial(P[:, worst]))
    log(
        f"Mitochondrial contribution to the embedding: |r| between pctMT and PC{worst + 1} = {simple[worst]:.3f}"
        f" (partial correlation with depth regressed out: {part:.3f}) / components with |r|>{threshold}: "
        f"{int((simple > threshold).sum())}/{P.shape[1]}"
    )
    if simple[worst] > threshold:
        warn(
            f"A principal component is strongly correlated with mitochondrial"
            f" content (|r|={simple[worst]:.3f}). Because the neighbor structure"
            " underlies QC stratification, doublet detection, and the SEACells"
            " kernel, this correlation propagates through the whole downstream"
            " pipeline. Consider adding --exclude-ribosomal-from-hvg, or"
            " regressing out pctMT during preprocessing."
        )
    return float(simple[worst])


def _infer_sample_id(path: str, index: int) -> str:
    """Infer a sample ID from a CellRanger output path.

    Assumes a structure like .../<sample>/outs/filtered_feature_bc_matrix;
    returns sample{index+1} if nothing usable is found.
    """
    parts = [p for p in Path(path).resolve().parts if p not in ("/", "")]
    ignore = {"filtered_feature_bc_matrix", "raw_feature_bc_matrix", "outs", "count"}
    for part in reversed(parts):
        if part.lower() not in ignore:
            return part
    return f"sample{index + 1}"


def ensure_sample_id(adata: ad.AnnData, default: str = "sample1") -> None:
    """Ensure obs['sample_id'] exists."""
    if "sample_id" not in adata.obs:
        warn(
            f"obs['sample_id'] does not exist. Setting it explicitly to"
            f" '{default}' since QC stratification and DE covariates depend on it."
        )
        adata.obs["sample_id"] = default
    adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
    log(f"sample_id: {list(adata.obs['sample_id'].cat.categories)}")


# ---------------------------------------------------------------------------
# Step 1b. Parametric QC filter (MAD-based adaptive thresholds)
# ---------------------------------------------------------------------------

def mad_outlier_mask(
    values: np.ndarray,
    nmads: float = 3.0,
    log: bool = True,
    direction: str = "both",  # "both" | "higher" | "lower"
) -> np.ndarray:
    """MAD-based outlier detection equivalent to scuttle::isOutlier. True = outlier (exclusion candidate).

    Uses the distribution's median and spread to set an adaptive threshold
    rather than a fixed cutoff (e.g. nCount_RNA > 20000), so it is less
    arbitrary across sequencing depths and samples. Note this detects
    outliers relative to the given population -- if the population itself
    contains a large distinct subgroup (e.g. a malignant clone), this must
    be applied stratified (per sample, per coarse cluster).
    """
    x = np.log1p(values) if log else np.asarray(values, dtype=float)
    x = np.asarray(x, dtype=float)
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826  # scaled to match SD under a normal-distribution assumption

    if mad == 0 or not np.isfinite(mad):
        # Degenerate distribution (all identical values): report no outliers
        return np.zeros_like(x, dtype=bool)

    lower = med - nmads * mad
    upper = med + nmads * mad

    if direction == "higher":
        return x > upper
    if direction == "lower":
        return x < lower
    return (x < lower) | (x > upper)


def parametric_qc_filter(
    adata: ad.AnnData,
    sample_key: str | None = "sample_id",
    coarse_cluster_key: str | None = None,
    nmads: float = 5.0,
    pct_mt_nmads: float = 3.0,
    min_group_size: int = 20,
    use_pctmt: bool = True,
) -> pd.DataFrame:
    """Apply MAD-based adaptive filtering on nCount/nFeature/pctMT.

    Stratification is the core idea here:
    - Median/MAD are computed independently per sample via sample_key
      (to absorb differences in sequencing depth).
    - If coarse_cluster_key is given, they are further computed
      independently per coarse cluster. This prevents a real biological
      shift (e.g. "the malignant cluster has uniformly higher RNA content")
      from being flagged as a technical outlier.

    High-count cells are evaluated more leniently (nmads) and are not
    included in `qc_fail`. Only low count, low feature count, and high
    pctMT are used for primary exclusion. Groups with fewer than
    min_group_size cells are skipped since their statistics would be unstable.

    Returns: a DataFrame with per-cell outlier=True/False columns.
    """
    qc_vars = ["mt"] if "mt" in adata.var else []
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True, percent_top=None)

    if "pct_counts_mt" not in adata.obs:
        if "mt" not in adata.var:
            flag_mito_genes(adata)  # fallback when GTF information is unavailable
        if adata.var["mt"].any():
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

    group_cols = [
        c for c in [sample_key, coarse_cluster_key] if c is not None and c in adata.obs
    ]
    if sample_key is not None and sample_key not in adata.obs:
        warn(
            f"obs['{sample_key}'] is missing, so QC stratification degrades"
            " to treating everything as one group, which may not be intended."
        )
    log(f"QC stratification keys: {group_cols if group_cols else '(no stratification)'}")

    result = pd.DataFrame(index=adata.obs_names)
    result["outlier_high_count"] = False
    result["outlier_low_count"] = False
    result["outlier_low_feature"] = False
    result["outlier_high_pctmt"] = False

    if group_cols:
        groups = adata.obs.groupby(group_cols, observed=True).indices
    else:
        groups = {"__all__": np.arange(adata.n_obs)}

    n_skipped = 0
    for key, idx in groups.items():
        idx = np.asarray(idx)
        if idx.size < min_group_size:
            n_skipped += 1
            continue
        sub = adata.obs.iloc[idx]

        # High-count side is evaluated leniently (nmads) to respect a real biological shift in a malignant cluster.
        result.iloc[idx, result.columns.get_loc("outlier_high_count")] = mad_outlier_mask(
            sub["total_counts"].values, nmads=nmads, direction="higher"
        )
        # Low-count side stays strict (~3 MAD) as usual; this flags empty/damaged cells.
        result.iloc[idx, result.columns.get_loc("outlier_low_count")] = mad_outlier_mask(
            sub["total_counts"].values, nmads=3.0, direction="lower"
        )
        result.iloc[idx, result.columns.get_loc("outlier_low_feature")] = mad_outlier_mask(
            sub["n_genes_by_counts"].values, nmads=3.0, direction="lower"
        )
        if "pct_counts_mt" in sub:
            result.iloc[idx, result.columns.get_loc("outlier_high_pctmt")] = mad_outlier_mask(
                sub["pct_counts_mt"].values, nmads=pct_mt_nmads, log=False, direction="higher"
            )
    if n_skipped:
        warn(f"Skipped QC evaluation for {n_skipped} group(s) with fewer than {min_group_size} cells")

    # High-count outliers can indicate malignancy, so they are not used
    # directly as an exclusion flag here; instead they are passed on as
    # input to downstream doublet calling and CNV-consistency checks.
    if use_pctmt:
        result["qc_fail"] = (
            result["outlier_low_count"]
            | result["outlier_low_feature"]
            | result["outlier_high_pctmt"]
        )
    else:
        warn(
            "pctMT-based exclusion is disabled (--no-pctmt-filter)."
            " This is appropriate for samples where the mtDNA composition is"
            " abnormal, since pctMT would not reflect true mitochondrial"
            " content there. The high_pctmt flag is still recorded for reference."
        )
        result["qc_fail"] = result["outlier_low_count"] | result["outlier_low_feature"]
    log(
        "QC: low_count={low}, low_feature={feat}, high_pctMT={mt}, "
        "high_count(not excluded)={high} / to be excluded={fail}".format(
            low=int(result["outlier_low_count"].sum()),
            feat=int(result["outlier_low_feature"].sum()),
            mt=int(result["outlier_high_pctmt"].sum()),
            high=int(result["outlier_high_count"].sum()),
            fail=int(result["qc_fail"].sum()),
        )
    )
    return result


# ---------------------------------------------------------------------------
# Step 1c. Cluster-aware doublet detection + CNV-consistency rescue
# ---------------------------------------------------------------------------

def detect_doublets_cluster_aware(
    adata: ad.AnnData,
    coarse_cluster_key: str = "coarse_cluster",
    expected_doublet_rate: float = 0.05,
    min_batch_size: int = 30,
    max_doublet_rate: float = 0.25,
) -> pd.DataFrame:
    """Run Scrublet per coarse cluster and return doublet_score / predicted_doublet.

    Using the coarse cluster as scanpy.pp.scrublet's batch_key keeps
    within-cluster (homotypic) pairs from dominating the simulated-doublet
    training set, concentrating detection power on between-cluster
    (heterotypic) doublets -- an approximation of the same idea behind
    Bioconductor's scDblFinder cluster-based mode.

    coarse_cluster_key should ideally be independent of the classification
    of interest (malignant/normal) -- e.g. a standard Leiden clustering result.

    Clusters with fewer than min_batch_size cells are merged into a single
    batch ('_small') before running Scrublet, since it is unstable on small batches.
    """
    if coarse_cluster_key not in adata.obs:
        raise KeyError(
            f"obs['{coarse_cluster_key}'] does not exist. Run coarse clustering"
            " before doublet detection."
        )

    # Build a minimal AnnData containing only raw counts to pass to Scrublet
    # (passing a normalized/log-transformed .X would give inaccurate doublet scores)
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    work = ad.AnnData(
        X=counts.copy(),
        obs=pd.DataFrame(index=adata.obs_names.copy()),
        var=pd.DataFrame(index=adata.var_names.copy()),
    )

    batch = adata.obs[coarse_cluster_key].astype(str)
    sizes = batch.value_counts()
    small = sizes[sizes < min_batch_size].index
    if len(small) > 0:
        warn(
            f"Merging coarse clusters with fewer than {min_batch_size} cells {list(small)} into the '_small' batch"
        )
        batch = batch.where(~batch.isin(small), "_small")
    work.obs["_doublet_batch"] = batch.values

    try:
        sc.pp.scrublet(
            work,
            batch_key="_doublet_batch",
            expected_doublet_rate=expected_doublet_rate,
            random_state=0,
        )
    except Exception as exc:
        warn(
            f"Cluster-stratified Scrublet failed ({exc}). Re-running without"
            " stratification (interpret results with caution, since"
            " homotypic pairs may be over-detected)."
        )
        work = ad.AnnData(
            X=counts.copy(),
            obs=pd.DataFrame(index=adata.obs_names.copy()),
            var=pd.DataFrame(index=adata.var_names.copy()),
        )
        sc.pp.scrublet(work, expected_doublet_rate=expected_doublet_rate, random_state=0)

    calls = work.obs[["doublet_score", "predicted_doublet"]].copy()
    calls["predicted_doublet"] = calls["predicted_doublet"].astype(bool)
    calls = calls.loc[adata.obs_names]  # align row order with the input
    n_called = int(calls["predicted_doublet"].sum())
    observed_rate = n_called / max(len(calls), 1)
    log(
        f"Doublet detection: predicted_doublet={n_called} / {len(calls)} cells"
        f" ({observed_rate:.1%}, expected {expected_doublet_rate:.1%})"
    )

    # --- Guard against automatic-threshold blowups ---
    # Scrublet's automatic threshold selection (call_doublets) can fail when
    # the simulated-doublet score distribution turns out unimodal, calling
    # most cells doublets (Scrublet itself recommends manually checking the
    # threshold). Running per-cluster makes this more likely for smaller
    # batches. 10x's expected doublet rate is roughly 1-10% based on loaded
    # cell count, so a call rate well above that is treated as threshold
    # failure, and we fall back to a conservative rule: call only the top
    # expected_doublet_rate fraction of cells by score as doublets.
    if observed_rate > max_doublet_rate:
        k = max(1, int(round(expected_doublet_rate * len(calls))))
        scores = np.sort(calls["doublet_score"].values)[::-1]
        threshold = float(scores[k - 1])
        warn(
            f"Doublet call rate was {observed_rate:.1%}, above the"
            f" {max_doublet_rate:.0%} cap. This likely means Scrublet's"
            f" automatic threshold selection failed (10x's expected doublet"
            f" rate is typically 1-10%), so falling back to calling only the"
            f" top {expected_doublet_rate:.0%} of cells by score as doublets"
            f" (score >= {threshold:.4g}, {k} cells). Check the doublet_score"
            " distribution and adjust --expected-doublet-rate /"
            " --max-doublet-rate if needed."
        )
        calls["predicted_doublet"] = (calls["doublet_score"] >= threshold).values
        calls.attrs["doublet_threshold_fallback"] = threshold
        log(f"After fallback: predicted_doublet={int(calls['predicted_doublet'].sum())}")
    return calls


def load_mito_reference_profile(path: str | None, out_dir: Path | None = None) -> dict | None:
    """Resolve the mtDNA composition reference profile to use.

    None uses the built-in reference (measured on a normally-aligned dog
    sample, Case1). 'auto' looks for mito_gene_profile.csv in out_dir.

    If the file cannot be found or parsed, this **warns and falls back to
    the built-in profile** rather than stopping the run: the reference is
    only used for rank-correlation rule 3 of the composition check, and the
    built-in profile supports the same check, so a bad path for this
    optional feature shouldn't abort a run that can take tens of minutes.
    Note the built-in profile is exactly the measured Case1 values, so
    results are unchanged whether the Case1-derived mito_gene_profile.csv
    is passed explicitly or this argument is left at the default.
    """
    if path is None:
        return None
    if path == "auto":
        cand = (out_dir / "mito_gene_profile.csv") if out_dir else None
        if cand and cand.exists():
            path = str(cand)
        else:
            log(
                "--mito-reference-profile auto: no reference file found,"
                " using the built-in reference values"
            )
            return None

    p = Path(path)
    if not p.exists():
        warn(
            f"The file given to --mito-reference-profile does not exist: {path}"
        )
        warn(
            "  Continuing with the built-in reference values (measured"
            " composition from a normally-aligned dog sample, Case1). This"
            " option can be omitted to run the same check. To use measured"
            " values from a different sample, point this at the"
            " mito_gene_profile.csv written to that run's --out-dir"
            " (this file is an output, not an input)."
        )
        return None
    try:
        table = pd.read_csv(p, index_col=0)
        col = "share_pct" if "share_pct" in table.columns else table.columns[-1]
        profile = {str(k): float(v) for k, v in table[col].items()}
        overlap = set(profile) & set(MITO_REFERENCE_PROFILE)
        if len(overlap) < 8:
            warn(
                f"Gene names in {path} overlap with fewer than 8 mtDNA genes"
                f" (overlap: {len(overlap)}). Continuing with the built-in reference values."
            )
            return None
        log(f"Loaded mtDNA reference profile: {path} (column '{col}', {len(profile)} genes)")
        return profile
    except Exception as exc:
        warn(f"Failed to load {path} ({exc}). Continuing with the built-in reference values.")
        return None


def check_mito_composition(
    adata: ad.AnnData,
    reference_profile: dict[str, float] | None = None,
    out_path: Path | None = None,
) -> dict:
    """Check whether the composition of the 13 mtDNA genes is biologically plausible.

    Why this check exists: for a sample whose fastq was corrupted and had
    to be rescue-aligned from a BAM that used a reference excluding
    mitochondria, mt-derived reads were picked up by NUMTs in the nuclear
    genome, and only a subset were recovered as mt genes. The resulting
    composition was badly skewed -- ND6 accounted for 47% of mt counts
    (normally ~0.3%, near the bottom), while the usually-dominant COX3,
    ATP6, and COX2 were nearly absent.

    Notably, this sample's **overall mt fraction was in the normal range**
    (0.157% vs. 0.127% for a normal sample) -- looking at pctMT alone would
    not catch this; only the per-gene breakdown reveals it. And because
    pctMT is evaluated with a MAD-based relative threshold, the small
    residual mt counts that did survive get flagged en masse as
    "possibly dead cells".

    Three rules are used (all based on composition share, not absolute counts):
      1. ND6 ranks within the top MITO_ND6_RANK_ALERT genes by share
      2. More than MITO_MAX_NEAR_ZERO genes have a share below MITO_NEAR_ZERO_PCT%
      3. Rank correlation with the reference profile falls below MITO_MIN_SPEARMAN
    """
    import scipy.sparse as sp

    if "mt" not in adata.var or not bool(np.asarray(adata.var["mt"], dtype=bool).any()):
        log("No mtDNA genes identified; skipping the composition check")
        return {"status": "skipped"}

    log("=== mtDNA gene composition check ===")
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
    is_mt = np.asarray(adata.var["mt"], dtype=bool)
    per_gene = pd.Series(
        np.asarray(counts[:, is_mt].sum(axis=0)).ravel().astype(float),
        index=adata.var_names[is_mt].astype(str),
    )
    total_mt = float(per_gene.sum())
    total_all = float(counts.sum())
    if total_mt <= 0:
        warn("mtDNA counts are 0; cannot run the composition check.")
        return {"status": "no_counts"}

    share = (per_gene / total_mt * 100).sort_values(ascending=False)
    log(f"mtDNA share of total counts: {total_mt / total_all:.4%}")
    log("Per-gene composition (%):")
    print(share.round(2).to_string())

    ref = dict(reference_profile or MITO_REFERENCE_PROFILE)
    problems: list[str] = []

    # Rule 1: ND6 rank
    ranks = {g: i + 1 for i, g in enumerate(share.index)}
    nd6_rank = ranks.get("ND6")
    if nd6_rank is not None and nd6_rank <= MITO_ND6_RANK_ALERT:
        problems.append(
            f"ND6 ranks #{nd6_rank} by composition share ({share.get('ND6', 0):.1f}%)."
            " ND6 is transcribed from the light strand, so it is normally near"
            " the bottom (reference: 0.3%); ranking near the top has no"
            " biological explanation"
        )

    # Rule 2: number of near-undetected genes
    near_zero = share.index[share < MITO_NEAR_ZERO_PCT].tolist()
    if len(near_zero) > MITO_MAX_NEAR_ZERO:
        problems.append(
            f"{len(near_zero)} genes have a composition share below"
            f" {MITO_NEAR_ZERO_PCT}% ({near_zero}). Because mtDNA is"
            " transcribed as a single polycistronic transcript, individual"
            " genes disappearing like this is not expected"
        )

    # Rule 3: rank correlation with the reference profile
    rho = float("nan")
    common = [g for g in share.index if g in ref]
    if len(common) >= 8:
        try:
            from scipy.stats import spearmanr

            rho = float(spearmanr([share[g] for g in common], [ref[g] for g in common])[0])
            log(f"Rank correlation with the reference profile: Spearman rho = {rho:.3f}")
            if np.isfinite(rho) and rho < MITO_MIN_SPEARMAN:
                problems.append(
                    f"Rank correlation with the reference profile is low"
                    f" ({rho:.3f}, threshold {MITO_MIN_SPEARMAN}). Samples"
                    " measuring the same underlying biology should correlate strongly"
                )
        except Exception as exc:  # pragma: no cover
            log(f"Skipped rank-correlation computation: {exc}")

    if out_path is not None:
        pd.DataFrame(
            {"counts": per_gene.astype(int), "share_pct": (per_gene / total_mt * 100).round(3)}
        ).sort_values("share_pct", ascending=False).to_csv(out_path)
        log(f"Saved mtDNA composition: {out_path.name}"
            " (can be passed to --mito-reference-profile in a future run)")

    status = "suspect" if problems else "ok"
    if problems:
        warn("mtDNA composition looks biologically implausible:")
        for i, msg in enumerate(problems, 1):
            warn(f"  ({i}) {msg}")
        warn(
            "Possible causes: alignment against a reference excluding"
            " mitochondria, sequence rescue from a BAM, mismapping to NUMTs,"
            " or a mismatch between the reference GTF and the BAM. In this"
            " state, pctMT does not reflect true mitochondrial content, so"
            " disable pctMT-based QC with --no-pctmt-filter, and do not"
            " interpret the mitochondrial/nuclear balance anomaly check (Step 3b)."
        )
    else:
        log("Composition looks plausible (consistent with the known mtDNA transcription pattern)")
    return {
        "status": status,
        "share_pct": share.to_dict(),
        "spearman_to_reference": rho,
        "problems": problems,
        "mt_fraction": total_mt / total_all,
    }


# ---------------------------------------------------------------------------
# Step 4a. CNV estimation branch (infercnvpy)
# ---------------------------------------------------------------------------

def parse_gtf_gene_positions(gtf_path: str, gene_id_type: str = "gene") -> pd.DataFrame:
    """Lightweight, self-contained GTF parser that extracts per-gene chromosome/start/end.

    infercnvpy.io.genomic_position_from_gtf()'s internal attribute parser
    can fail to extract attribute keys for non-GENCODE formats (NCBI/RefSeq,
    CellRanger's filtered GTF, etc.), so this is provided as a fallback.
    Only 'gene' rows are parsed, with the attribute column (field 9) parsed
    as ' key "value";' pairs (only the standard GTF format is supported;
    GFF3 is not).
    """
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    records = []

    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            chrom, start, end, attr_str = fields[0], fields[3], fields[4], fields[8]
            attrs = {}
            for kv in attr_str.strip().split(";"):
                kv = kv.strip()
                if not kv or " " not in kv:
                    continue
                k, v = kv.split(" ", 1)
                attrs[k] = v.strip('"')

            gene_key = attrs.get(gene_id_type)
            if gene_key is None:
                continue
            records.append((gene_key, chrom, int(start), int(end)))

    if not records:
        raise ValueError(
            f"Could not extract any gene rows with a '{gene_id_type}' attribute"
            " from the GTF. Confirm the file is valid GTF and check the"
            " attribute key name (GENCODE: gene_name / NCBI-RefSeq: gene)."
        )

    df = pd.DataFrame(records, columns=["gene", "chromosome", "start", "end"])
    # If the same gene appears in multiple rows (e.g. different versions), keep the first occurrence
    df = df.drop_duplicates(subset="gene", keep="first").set_index("gene")
    return df


def _write_h5ad(adata: ad.AnnData, path: Path) -> None:
    """Save an h5ad file. Mixed object dtypes in var/obs can cause this to
    fail; on failure, retry after stringifying those columns, and if that
    still fails, warn and continue rather than raising."""
    try:
        adata.write_h5ad(path)
        log(f"Saved: {path}")
        return
    except Exception as exc:
        warn(f"Failed to save {path.name} ({exc}). Retrying after stringifying columns.")
    try:
        copy = adata.copy()
        for frame in (copy.obs, copy.var):
            for col in frame.columns:
                if frame[col].dtype == object:
                    frame[col] = frame[col].astype(str)
        copy.write_h5ad(path)
        log(f"Saved: {path}")
    except Exception as exc:  # pragma: no cover
        warn(f"Giving up saving {path.name}: {exc}")


def _leiden_supports_igraph() -> bool:
    """Check whether scanpy's sc.tl.leiden accepts flavor='igraph'
    (used to avoid a FutureWarning on scanpy >= 1.10)."""
    import inspect

    try:
        return "flavor" in inspect.signature(sc.tl.leiden).parameters
    except (TypeError, ValueError):  # pragma: no cover
        return False



# ---------------------------------------------------------------------------
# Externalized marker panels (v3.3)
# ---------------------------------------------------------------------------
# Hardcoding the markers used for cell-type calling would require editing
# code every time the species changes. Instead they live in
# `markers/markers_<species>.csv` and can be swapped via `--markers`.
# Columns:
#   cell_type, gene, role, use_as_normal_reference, use_for_labels, validation
# role                     primary (default) / fallback (used only if primary markers are absent)
# use_as_normal_reference  whether this cell type may be used as the CNV normal
#                          reference (Fibroblast / Epithelial can themselves be
#                          malignant, so FALSE for those)
# use_for_labels           whether this type may be used for metacell labeling
#                          (validated types only)
# validation               depth_matched_case1_dog / orthology_conversion
#                          (the latter triggers a warning if used for labeling,
#                          since it is unvalidated in this project)

BUILTIN_MARKER_SPECIES = ("dog", "human", "mouse")


class MarkerTable:
    """Marker table. A thin wrapper exposing different views of the panel as properties."""

    def __init__(self, df: pd.DataFrame, source: str):
        self.df = df
        self.source = source

    def _sets(self, mask) -> dict[str, list[str]]:
        sub = self.df[mask]
        return {ct: list(dict.fromkeys(g["gene"].tolist()))
                for ct, g in sub.groupby("cell_type", sort=False)}

    @property
    def marker_sets(self) -> dict[str, list[str]]:
        return self._sets(self.df["role"] == "primary")

    @property
    def fallback_sets(self) -> dict[str, list[str]]:
        return self._sets(self.df["role"] == "fallback")

    @property
    def normal_celltypes(self) -> list[str]:
        m = self.df["use_as_normal_reference"]
        return list(dict.fromkeys(self.df.loc[m, "cell_type"].tolist()))

    @property
    def label_panels(self) -> dict[str, list[str]]:
        """Panels approved for labeling (validated)."""
        return self._sets(self.df["use_for_labels"])

    @property
    def reported_panels(self) -> dict[str, list[str]]:
        """Panels that are scored but not used for labeling."""
        return self._sets(~self.df["use_for_labels"] & (self.df["role"] == "primary"))


def resolve_marker_path(spec: str | None, search_dir: str | Path | None = None) -> Path:
    """Resolve a species name like 'dog', or a CSV path, to an actual file."""
    base = Path(search_dir) if search_dir else Path(__file__).resolve().parent
    if spec is None:
        spec = "dog"
    p = Path(spec)
    if p.suffix.lower() == ".csv":
        if not p.exists():
            raise FileNotFoundError(f"Marker CSV not found: {p}")
        return p
    name = str(spec).strip().lower()
    cand = base / "markers" / f"markers_{name}.csv"
    if cand.exists():
        return cand
    raise FileNotFoundError(
        f"Marker table '{spec}' not found. Built-in options are"
        f" {', '.join(BUILTIN_MARKER_SPECIES)}. A CSV path may also be given"
        f" (looked for: {cand})")


def load_marker_table(spec: str | None = "dog",
                      search_dir: str | Path | None = None) -> MarkerTable:
    """Load a marker table, checking for missing columns and unvalidated panels used for labeling."""
    path = resolve_marker_path(spec, search_dir)
    df = pd.read_csv(path)
    need = {"cell_type", "gene", "role", "use_as_normal_reference",
            "use_for_labels", "validation"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"{path.name} is missing required columns: {sorted(missing)}")
    for col in ("use_as_normal_reference", "use_for_labels"):
        df[col] = (df[col].astype(str).str.strip().str.lower()
                   .isin({"true", "1", "yes", "y"}))
    df["cell_type"] = df["cell_type"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()
    df["role"] = df["role"].astype(str).str.strip().str.lower()
    bad = set(df["role"]) - {"primary", "fallback"}
    if bad:
        raise ValueError(f"role must be 'primary' or 'fallback': {sorted(bad)}")

    tbl = MarkerTable(df, str(path))
    n_lab = len(tbl.label_panels)
    log(f"Marker table: {path.name} / {df['cell_type'].nunique()} types /"
        f" {len(df)} gene rows / {n_lab} types approved for labeling /"
        f" {len(tbl.normal_celltypes)} types usable as normal reference")
    unval = df[df["use_for_labels"] & (df["validation"] != "depth_matched_case1_dog")]
    if len(unval):
        warn("Some panels are marked for labeling but are unvalidated in this project: "
             + ", ".join(sorted(unval["cell_type"].unique()))
             + ". Panels built purely by naming-convention conversion (human, mouse)"
             " are in this state. Confirm specificity before relying on their labels")
    if n_lab == 0:
        warn("No cell types have use_for_labels set to TRUE;"
             " metacell labeling cannot proceed")
    return tbl
