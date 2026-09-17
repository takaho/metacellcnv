#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""metacellcnv_scanpy.py -- Seurat-free unified preprocessing module.

Runs cell selection -> normalization -> HVG -> PCA -> UMAP -> clustering
entirely in scanpy (previously done in R/Seurat), so that downstream steps
(the CNV/DE pipeline, metacell splitting, R-side visualization) all read
the same cell set, gene set, and embedding.

Depends only on `scrna_common.py`; does not import SEACells, infercnvpy,
or pyDESeq2, so it can be used standalone for preprocessing.

Key differences from the legacy pipeline:
  * Mitochondrial genes are identified from the GTF sequence ID (default
    NC_002008.4), not from a `MT-` name prefix -- needed because that
    prefix is a primate convention, and non-primate mtDNA gene symbols
    (e.g. dog COX1, ND1, ...) don't carry it. Uses the same CLI flag names
    as the main pipeline (--gtf, --gtf-gene-id, --mito-chromosome,
    --mito-reference-profile, --mito-nmads, --keep-mito-in-hvg,
    --exclude-ribosomal-from-hvg, --no-pctmt-filter), so the same command
    line can be reused.
  * All outputs go into one directory with a `prep_manifest.json` that
    records every parameter and each file's sha256, so which analyses
    used which cell set can be checked mechanically after the fact.
  * `clusters.tsv.gz` from this module is the canonical cluster
    assignment for metacell splitting; containment of an existing
    metacell assignment within these clusters can be measured with
    `metacell_cluster_containment()` / `--containment`.

Numerical relationship to Seurat: this module does not try to reproduce
Seurat's output. Defaults match the existing pipeline; pass
`--seurat-compat` to move toward Seurat's defaults instead. Differences:

| Aspect | Default here (= existing pipeline) | Seurat default (`--seurat-compat`) |
|---|---|---|
| Normalization scale | `normalize_total()` = median counts | `scale.factor=1e4` |
| HVG | scanpy `seurat` flavor | `vst` = scanpy `seurat_v3` (needs `scikit-misc`) |
| Scaling before PCA | none (PCA on log1p values) | `ScaleData`-equivalent (z-score, clip at 10) |
| UMAP | `n_neighbors=15, min_dist=0.5` | `n_neighbors=30, min_dist=0.3` |

Clustering always uses Leiden here (Seurat uses Louvain/SLM), so cluster
counts will not match exactly; `--cluster-algo louvain` gets closer, but
SLM has no scanpy equivalent.

Pipeline steps:
  1. Load + keep raw counts + identify mtDNA + normalize + HVG + PCA (all cells)
  2. Provisional clustering (for QC stratification)
  3. mtDNA composition check
  4. MAD-based QC (stratified by sample x provisional cluster)
  5. Doublet detection (optional)
  6. Recompute normalization/HVG/PCA on the final cell set (`--no-refit` to disable)
  7. Neighbor graph -> clustering -> UMAP (3D by default)
  8. Export

Step 6 has no equivalent in the existing pipeline, which keeps using the
HVGs and PCs chosen from all cells before QC, so excluded cells still
influence gene selection and the principal components. Recomputing the
embedding from only the final cell set is more self-consistent, so it is
the default; pass `--no-refit` to match the existing pipeline's numbers
exactly.

See README.md for further background.

Usage
-----
    python3 metacellcnv_scanpy.py --cellranger-dir <sample>/outs/filtered_feature_bc_matrix \
        --sample-id <sample> --gtf genes.gtf --mito-chromosome <mito_contig> \
        --out-dir prep/<sample>

    # Measure containment of an existing metacell assignment only
    python3 metacellcnv_scanpy.py --containment results/<sample>/cell_to_metacell.csv \
        --clusters prep/<sample>/clusters.tsv.gz
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

__version__ = "1.1"

# ---------------------------------------------------------------------------
# Uses shared helpers from scrna_common.py. Does not depend on the CNV/DE
# analysis pipeline, so this module can be used standalone for
# preprocessing only.
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

try:
    import scrna_common as common
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "Could not import scrna_common.py."
        f" Place it in the same directory: {exc}"
    )

import scanpy as sc  # noqa: E402  (fine to import after common)

log = common.log
warn = common.warn

# Seurat's defaults (used when --seurat-compat is passed)
SEURAT_SCALE_FACTOR = 1e4
SEURAT_UMAP_NEIGHBORS = 30
SEURAT_UMAP_MIN_DIST = 0.3
SEURAT_SCALE_MAX = 10.0
SEURAT_MIN_CELLS_PER_GENE = 3


# ===========================================================================
# Small utilities
# ===========================================================================
def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    """Hash used to verify outputs; lets later steps detect cell-set mismatches."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def write_tsv_gz(df: pd.DataFrame, path: Path, index: bool = True,
                 index_label: str | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as fh:
        df.to_csv(fh, sep="\t", index=index, index_label=index_label)
    return path


def write_lines_gz(values, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as fh:
        for v in values:
            fh.write(f"{v}\n")
    return path


def write_mtx_gz(matrix, path: Path, orient: str = "cells-genes",
                 integer: bool = False) -> Path:
    """Save a sparse matrix in MatrixMarket format (.mtx.gz).

    orient='cells-genes' writes n_cells x n_genes; 'genes-cells' transposes
    it to match the orientation used by 10x / Seurat's ReadMtx. The
    orientation is also recorded in the manifest.
    """
    if orient not in ("cells-genes", "genes-cells"):
        raise ValueError(f"Invalid orient: {orient}")
    m = sp.csr_matrix(matrix)
    if orient == "genes-cells":
        m = m.T
    m = m.tocoo()
    if integer:
        data = np.rint(m.data)
        if np.abs(data - m.data).max(initial=0.0) > 1e-8:
            raise ValueError("integer=True but values are not integers")
        m = sp.coo_matrix((data.astype(np.int64), (m.row, m.col)), shape=m.shape)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as fh:
        sio.mmwrite(fh, m, field="integer" if integer else "real",
                    comment=f" orientation={orient} (rows x cols)")
    return path


def _timed(label: str, fn, *args, **kwargs):
    t0 = time.time()
    out = fn(*args, **kwargs)
    log(f"{label}: {time.time() - t0:.1f}s")
    return out


# ===========================================================================
# Clustering and UMAP
# ===========================================================================
def build_neighbors(adata, n_neighbors: int = 15, n_pcs: int | None = None,
                    use_rep: str = "X_pca", seed: int = 0) -> None:
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                    use_rep=use_rep, random_state=seed)


def cluster_cells(adata, resolution: float = 0.5, key_added: str = "cluster",
                  algo: str = "leiden", seed: int = 0) -> pd.Series:
    """Cluster cells with Leiden (default) or Louvain.

    Uses the same calling convention as `coarse_cluster` in the main
    pipeline (scanpy >=1.10 recommends flavor='igraph'; older versions
    don't have that argument at all).
    """
    if algo == "leiden":
        kwargs = {"resolution": resolution, "key_added": key_added,
                  "random_state": seed}
        if common._leiden_supports_igraph():
            kwargs.update({"flavor": "igraph", "n_iterations": 2, "directed": False})
        try:
            sc.tl.leiden(adata, **kwargs)
        except (ImportError, TypeError, ValueError) as exc:
            warn(f"Recommended leiden settings unavailable, retrying with defaults: {exc}")
            sc.tl.leiden(adata, resolution=resolution, key_added=key_added,
                         random_state=seed)
    elif algo == "louvain":
        # Seurat's FindClusters offers Louvain/SLM; scanpy has no SLM, so Louvain is as close as it gets.
        try:
            sc.tl.louvain(adata, resolution=resolution, key_added=key_added,
                          random_state=seed)
        except ImportError as exc:
            raise SystemExit(
                "--cluster-algo louvain requires the louvain package"
                f" (pip install louvain). Use --cluster-algo leiden for Leiden instead: {exc}"
            )
    else:
        raise ValueError(f"Invalid --cluster-algo: {algo}")
    n = adata.obs[key_added].nunique()
    log(f"Clusters ({algo}, resolution={resolution}): {n}")
    return adata.obs[key_added]


def embed_umap(adata, n_components: int = 3, n_neighbors: int | None = None,
               min_dist: float = 0.5, seed: int = 0,
               key_added: str = "X_umap") -> np.ndarray:
    """Compute a UMAP embedding, 3D by default (per spec).

    If n_neighbors is given, rebuilds the neighbor graph for UMAP.
    Otherwise reuses the graph built for clustering, so clusters and UMAP
    come from the same graph and their boundaries stay consistent.
    """
    if n_neighbors is not None:
        build_neighbors(adata, n_neighbors=n_neighbors, seed=seed)
    sc.tl.umap(adata, n_components=n_components, min_dist=min_dist,
               random_state=seed)
    if key_added != "X_umap":
        adata.obsm[key_added] = adata.obsm["X_umap"]
    return np.asarray(adata.obsm["X_umap"])


# ===========================================================================
# Recompute the embedding (on the final cell set)
# ===========================================================================
def refit_embedding(adata, n_hvg: int, n_pcs: int, *,
                    norm_target_sum: float | None = None,
                    hvg_flavor: str = "seurat",
                    scale_before_pca: bool = False,
                    scale_max_value: float | None = None,
                    exclude_mito_from_hvg: bool = True,
                    exclude_ribo_from_hvg: bool = False,
                    batch_key: str | None = "sample_id") -> None:
    """Recompute normalization, HVG, and PCA from `layers['counts']`.

    Call this after QC and doublet exclusion so the exported embedding is
    determined only by the final cell set. mtDNA / ribosomal genes are
    excluded from HVG for the same reason as in the main pipeline (to keep
    the embedding from being driven by cell state rather than cell identity).
    """
    if "counts" not in adata.layers:
        raise ValueError("layers['counts'] is missing. Raw counts must be kept.")
    adata.X = adata.layers["counts"].copy()
    # Clear any log1p record left in uns from the first pass, or scanpy will
    # wrongly warn that the data is already log-transformed (X has been
    # reset to raw counts here).
    adata.uns.pop("log1p", None)

    if norm_target_sum is not None:
        sc.pp.normalize_total(adata, target_sum=norm_target_sum)
    else:
        sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    hvg_kwargs: dict = {"n_top_genes": n_hvg, "flavor": hvg_flavor}
    if batch_key and batch_key in adata.obs and adata.obs[batch_key].nunique() > 1:
        hvg_kwargs["batch_key"] = batch_key
    if hvg_flavor == "seurat_v3":
        # seurat_v3 requires raw counts
        adata.X = adata.layers["counts"].copy()
        try:
            sc.pp.highly_variable_genes(adata, **hvg_kwargs)
        except ImportError as exc:
            raise SystemExit(
                "HVG flavor='seurat_v3' requires scikit-misc"
                f" (pip install scikit-misc): {exc}"
            )
        adata.X = adata.layers["lognorm"].copy()
    else:
        sc.pp.highly_variable_genes(adata, **hvg_kwargs)

    hv = np.asarray(adata.var["highly_variable"], dtype=bool)
    excluded: list[str] = []
    if exclude_mito_from_hvg and "mt" in adata.var:
        is_mt = np.asarray(adata.var["mt"], dtype=bool)
        if (hv & is_mt).any():
            excluded += list(adata.var_names[hv & is_mt])
            hv &= ~is_mt
    if exclude_ribo_from_hvg:
        upper = adata.var_names.str.upper()
        is_ribo = np.zeros(adata.n_vars, dtype=bool)
        for prefix in common.RIBO_PREFIXES:
            is_ribo |= np.asarray(upper.str.startswith(prefix), dtype=bool)
        if (hv & is_ribo).any():
            excluded += list(adata.var_names[hv & is_ribo])
            hv &= ~is_ribo
    if excluded:
        log(f"Excluded {len(excluded)} genes from HVG: {sorted(excluded)[:12]}"
            f"{' ...' if len(excluded) > 12 else ''}")
    adata.var["hvg_for_pca"] = hv
    log(f"HVGs used for PCA: {int(hv.sum())}")

    if scale_before_pca:
        # Equivalent to Seurat's ScaleData. Off by default (main pipeline does PCA on log1p values).
        sub = adata[:, hv].copy()
        sc.pp.scale(sub, max_value=scale_max_value)
        sc.tl.pca(sub, n_comps=n_pcs)
        adata.obsm["X_pca"] = sub.obsm["X_pca"]
        adata.uns["pca"] = sub.uns.get("pca", {})
        del sub
    else:
        try:
            sc.tl.pca(adata, n_comps=n_pcs, mask_var="hvg_for_pca")
        except TypeError:  # pragma: no cover - old scanpy
            adata.var["highly_variable"] = hv
            sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)
    common.check_embedding_mito_influence(adata)


# ===========================================================================
# Metacell / cluster containment
# ===========================================================================
def metacell_cluster_containment(assignments: pd.Series | pd.DataFrame,
                                 clusters: pd.Series) -> dict:
    """Measure how well each metacell is contained within a single cluster.

    Quantifies the requirement that metacells be interpretable as subsets
    of conventional clusters. Boundary cells always exist, so 1.0 is not
    expected.

    containment(m) = max_c |{cells in m} intersect {cells in cluster c}| / |{cells in m}|

    Returns:
      per_metacell : per-metacell containment / dominant cluster / cell count / clusters spanned
      summary      : weighted-mean containment, fraction fully contained, spanning distribution
    """
    if isinstance(assignments, pd.DataFrame):
        col = None
        # "metacell" is the new name; "SEACell" is kept for compatibility with
        # results from before the rename; "metacell_id"/"seacell" are aliases
        # this function has accepted on its own.
        for cand in ("metacell", "SEACell", "metacell_id", "seacell"):
            if cand in assignments.columns:
                col = cand
                break
        if col is None:
            raise ValueError(
                f"No metacell column found (expected one of: SEACell/metacell/metacell_id)."
                f" Actual columns: {list(assignments.columns)}"
            )
        assignments = assignments[col]
    assign = pd.Series(assignments).astype(str)
    clu = pd.Series(clusters).astype(str)

    common = assign.index.intersection(clu.index)
    if len(common) == 0:
        raise ValueError("No overlapping barcodes between the metacell assignment and clusters")
    if len(common) < len(assign):
        warn(f"Barcode overlap {len(common)}/{len(assign)}."
             " Check prep_manifest.json to confirm these are outputs from the same cell set.")
    assign, clu = assign.loc[common], clu.loc[common]

    ct = pd.crosstab(assign, clu)
    sizes = ct.sum(axis=1)
    top = ct.max(axis=1)
    per = pd.DataFrame({
        "n_cells": sizes,
        "containment": top / sizes,
        "dominant_cluster": ct.idxmax(axis=1),
        "n_clusters_spanned": (ct > 0).sum(axis=1),
    }).sort_values("containment")
    weighted = float((per["containment"] * per["n_cells"]).sum() / per["n_cells"].sum())
    summary = {
        "n_metacells": int(len(per)),
        "n_cells": int(per["n_cells"].sum()),
        "n_clusters": int(clu.nunique()),
        "containment_weighted_mean": weighted,
        "containment_median": float(per["containment"].median()),
        "frac_fully_contained": float((per["containment"] >= 1.0).mean()),
        "frac_containment_ge_0.9": float((per["containment"] >= 0.9).mean()),
        "frac_containment_ge_0.8": float((per["containment"] >= 0.8).mean()),
        "spanned_clusters_distribution": {
            str(k): int(v) for k, v in
            per["n_clusters_spanned"].value_counts().sort_index().items()
        },
    }
    return {"per_metacell": per, "summary": summary}


def report_containment(result: dict) -> None:
    s = result["summary"]
    log(f"Metacell/cluster containment ({s['n_metacells']} metacells / "
        f"{s['n_cells']:,} cells / {s['n_clusters']} clusters)")
    log(f"  weighted mean containment : {s['containment_weighted_mean']:.3f}")
    log(f"  median                    : {s['containment_median']:.3f}")
    log(f"  fully contained fraction  : {s['frac_fully_contained']:.1%}")
    log(f"  >=0.9 / >=0.8             : {s['frac_containment_ge_0.9']:.1%}"
        f" / {s['frac_containment_ge_0.8']:.1%}")
    log(f"  clusters-spanned distribution : {s['spanned_clusters_distribution']}")
    worst = result["per_metacell"].head(5)
    if len(worst):
        log("  lowest-containment metacells:")
        for mc, row in worst.iterrows():
            log(f"    {mc}: {row['containment']:.2f}"
                f" (n={int(row['n_cells'])}, spans {int(row['n_clusters_spanned'])} clusters)")


# ===========================================================================
# MT-derived mRNA dead-cell judgment (adaptive threshold)
# ===========================================================================
#
# Why a fixed threshold (e.g. pctMT 5%) is not used here:
#
# (a) The baseline mitochondrial fraction varies by orders of magnitude
#     across species, reference genome, and alignment. Non-primate mtDNA
#     gene symbols (e.g. dog COX1, ND1) also don't carry the `MT-` prefix
#     that a name-based fixed filter (like Seurat's percent.mt regex)
#     expects, so such a filter can silently flag almost nothing.
#
# (b) A fixed threshold is depth-biased: at low total counts, sampling
#     noise alone can push a healthy cell's pctMT well above threshold,
#     while the same cell at higher depth would pass.
#
# (c) High pctMT does not necessarily mean a dead cell. The expected
#     signature of cell death (membrane rupture) is falling nuclear-
#     derived counts alongside a roughly unchanged mtDNA count; a cell
#     that is simply mitochondria-rich instead shows mtDNA counts rising
#     while nuclear counts stay flat. pctMT rises in both cases, so the
#     ratio alone cannot tell them apart.
#
# (d) Mitochondrial fraction is strongly overdispersed across cells. A
#     plain binomial null flags far too many cells as significant; the
#     null needs to account for that overdispersion (beta-binomial), or
#     a large fraction of otherwise normal cells looks like an outlier.
#
# So this module only marks cells dead when ALL FOUR gates below pass,
# and by default only diagnoses and recommends -- it removes nothing.
#
#   Gate 1 measurability    mtDNA genes identified, with enough counts to judge
#   Gate 2 separability     outliers actually exist under the (overdispersion-aware) null
#   Gate 3 mechanism fit    flagged cells look like they lost cytoplasmic mRNA, not just gained mtDNA
#   Gate 4 non-confounding  the flagged set isn't enriched for malignant / doublet / one cluster

# Gate defaults (all overridable from the CLI)
MITO_MIN_GENES = 5            # minimum number of identified mtDNA genes
MITO_MIN_MEDIAN_COUNTS = 3    # minimum median mtDNA count (below this, noise dominates)
MITO_MAX_ZERO_FRAC = 0.20     # max fraction of cells with zero mtDNA counts
MITO_FDR = 0.01               # FDR for the beta-binomial test
MITO_MIN_FLAGGED = 20         # minimum number of flagged cells (absolute)
MITO_MIN_FLAGGED_FRAC = 1e-3  # minimum fraction of cells flagged
MITO_MAX_FLAGGED_FRAC = 0.30  # if too many cells are flagged, suspect the null model
MITO_MAX_NUC_RATIO = 0.85     # max ratio of flagged cells' nuclear counts to the stratum median (cytoplasmic loss)
MITO_MAX_MT_RATIO = 3.0       # max ratio of flagged cells' mtDNA counts to the stratum median
                              # (membrane rupture leaves mtDNA roughly unchanged, not elevated)
MITO_MIN_RESID_DROP_MAD = 0.5 # min drop in complexity residual (MAD units) for flagged cells
MITO_CONFOUND_OR = 2.0        # odds ratio considered confounding
MITO_CONFOUND_P = 0.01        # significance level for the confounding Fisher test
MITO_MIN_STRATUM = 50         # strata smaller than this are pooled together
MITO_FIT_QUANTILE = 0.90      # upper quantile trimmed when fitting the null distribution (robustness)


def complexity_residual(total_counts, n_genes, deg: int = 2) -> np.ndarray:
    """Residual of detected-gene count after removing the part explained by depth.

    Dead cells lose cytoplasmic mRNA diversity, so they show fewer
    detected genes than other cells with the same total count. The
    relationship with total counts is strongly curved, so a degree-2 fit
    is done in log-log space before taking the residual.
    """
    x = np.log10(np.maximum(np.asarray(total_counts, dtype=float), 1.0))
    y = np.log10(np.maximum(np.asarray(n_genes, dtype=float), 1.0))
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < deg + 2:
        return np.zeros_like(y)
    coef = np.polyfit(x[ok], y[ok], deg)
    return y - np.polyval(coef, x)


def _betabinom_mom(k: np.ndarray, n: np.ndarray,
                   fit_quantile: float | None = 0.90) -> tuple[float, float] | None:
    """Method-of-moments estimate of beta-binomial (alpha, beta); None if no overdispersion.

    Cell-to-cell variation in mitochondrial fraction combines sampling
    noise (binomial) with genuine biological variation (overdispersion).
    Without absorbing the latter into the null distribution, most cells
    end up "significant".

    Fitting on the **upper-trimmed** distribution (fit_quantile) matters:
    because the moment estimator relies on the variance, a small number of
    dead cells can inflate the variance enough to widen the null and hide
    themselves (a failure mode of the moment method). Trimming the upper
    tail before fitting avoids this, at the cost of making the resulting
    FDR control approximate -- gates 3 (mechanism fit) and 4
    (non-confounding) cover that gap.
    """
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    ok = n > 0
    k, n = k[ok], n[ok]
    if len(n) < 10 or n.sum() <= 0:
        return None
    if fit_quantile is not None and 0 < fit_quantile < 1 and len(n) >= 20:
        frac = k / np.maximum(n, 1.0)
        keep = frac <= np.quantile(frac, fit_quantile)
        if keep.sum() >= 10:
            k, n = k[keep], n[keep]
    m = k.sum() / n.sum()
    if not (0 < m < 1):
        return None
    w = n / n.sum()
    v = float(np.average((k / n - m) ** 2, weights=w))
    n_h = float(len(n) / np.sum(1.0 / np.maximum(n, 1.0)))   # harmonic mean
    excess = v * n_h / max(m * (1 - m), 1e-12) - 1.0
    if excess <= 1e-9:
        return None                                  # explained by the binomial alone
    rho = max((n_h - 1.0) / excess - 1.0, 1e-3)      # rho = alpha + beta
    return m * rho, (1 - m) * rho


def _bh(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    q = np.empty(m)
    ranked = p[order] * m / np.arange(1, m + 1)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.clip(q, 0.0, 1.0)


def mito_dead_cell_report(
    total_counts, mt_counts, n_genes, *,
    stratum=None, index=None, n_mito_genes: int | None = None,
    labels: dict | None = None,
    fdr: float = MITO_FDR,
    min_genes: int = MITO_MIN_GENES,
    min_median_counts: float = MITO_MIN_MEDIAN_COUNTS,
    max_zero_frac: float = MITO_MAX_ZERO_FRAC,
    min_flagged: int = MITO_MIN_FLAGGED,
    min_flagged_frac: float = MITO_MIN_FLAGGED_FRAC,
    max_flagged_frac: float = MITO_MAX_FLAGGED_FRAC,
    max_nuc_ratio: float = MITO_MAX_NUC_RATIO,
    max_mt_ratio: float = MITO_MAX_MT_RATIO,
    min_resid_drop_mad: float = MITO_MIN_RESID_DROP_MAD,
    confound_or: float = MITO_CONFOUND_OR,
    confound_p: float = MITO_CONFOUND_P,
    min_stratum: int = MITO_MIN_STRATUM,
    fit_quantile: float = MITO_FIT_QUANTILE,
    composition_ok: bool | None = None,
) -> dict:
    """Evaluate MT-derived dead-cell status via four gates plus an adaptive threshold.

    Removes nothing. Returns the judgment and a recommendation only;
    adaptive_mito_filter() performs the actual removal, according to `mode`.

    `labels` carries boolean vectors used for the confounding check (e.g.
    {"malignant": bool array, "doublet": bool array}).

    Limitation of gate 4: when the flagged set correlates strongly with a
    malignancy label, this test cannot distinguish "MT flagging is
    preferentially catching tumor cells (should not be removed)" from
    "tumor cells are genuinely dying more often (removal is fine)". Since
    it cannot tell these apart, the default is to stop; proceed with
    mode="force" only if you have other evidence to decide.

    Returns a dict:
      per_cell   : DataFrame with pct_mt / p_value / q_value / flagged / complexity_resid
      gates      : gate name -> {"passed": bool, "detail": str, ...}
      recommend  : "remove" | "do_not_remove"
      reason     : human-readable justification (emitted to the log as-is)
      stats      : key summary statistics
    """
    from scipy import stats as _st

    tot = np.asarray(total_counts, dtype=float)
    mt = np.asarray(mt_counts, dtype=float)
    ng = np.asarray(n_genes, dtype=float)
    n_cells = len(tot)
    if not (len(mt) == len(ng) == n_cells):
        raise ValueError("total_counts / mt_counts / n_genes must be the same length")
    if index is None:
        index = pd.RangeIndex(n_cells)
    idx = pd.Index(index)

    pct_mt = 100.0 * mt / np.maximum(tot, 1.0)
    nuc = tot - mt
    resid = complexity_residual(tot, ng)

    if stratum is None:
        strat = pd.Series(["all"] * n_cells, index=idx)
    else:
        strat = pd.Series(np.asarray(stratum, dtype=object), index=idx).astype(str)
        counts = strat.value_counts()
        small = set(counts.index[counts < min_stratum])
        if small:
            strat = strat.where(~strat.isin(small), "__pooled__")

    # --- adaptive threshold: upper-tail probability under a per-stratum beta-binomial null ---
    pvals = np.ones(n_cells)
    fits: dict[str, dict] = {}
    for gname, pos in strat.groupby(strat).groups.items():
        sel = idx.get_indexer(pd.Index(pos))
        ab = _betabinom_mom(mt[sel], tot[sel], fit_quantile=fit_quantile)
        if ab is None:
            # Strata with no detectable overdispersion fall back to a binomial null
            # (baseline rate taken from the lower half only, so that a stratum that
            # is mostly dead cells doesn't inflate its own baseline)
            thr = np.quantile(pct_mt[sel], 0.5)
            base = sel[pct_mt[sel] <= thr]
            p0 = max(mt[base].sum() / max(tot[base].sum(), 1.0), 1e-6)
            pvals[sel] = _st.binom.sf(mt[sel] - 1, tot[sel], p0)
            fits[gname] = {"model": "binom", "p0": float(p0), "n": int(len(sel))}
        else:
            a, b = ab
            pvals[sel] = _st.betabinom.sf(mt[sel] - 1, tot[sel], a, b)
            fits[gname] = {"model": "betabinom", "alpha": float(a), "beta": float(b),
                           "rho": float(a + b), "mean_pct": float(100 * a / (a + b)),
                           "n": int(len(sel))}
    qvals = _bh(pvals)
    flagged = qvals < fdr

    # For reference: candidate count under the null fit without trimming. A large
    # difference from the robust (trimmed) version is evidence that the candidates
    # themselves were widening the null.
    pv_full = np.ones(n_cells)
    for gname, pos in strat.groupby(strat).groups.items():
        sel = idx.get_indexer(pd.Index(pos))
        ab = _betabinom_mom(mt[sel], tot[sel], fit_quantile=None)
        if ab is None:
            pv_full[sel] = pvals[sel]
        else:
            pv_full[sel] = _st.betabinom.sf(mt[sel] - 1, tot[sel], *ab)
    n_flag_full = int((_bh(pv_full) < fdr).sum())

    per_cell = pd.DataFrame(
        {"total_counts": tot, "mt_counts": mt, "n_genes": ng,
         "pct_mt": pct_mt, "complexity_resid": resid,
         "stratum": strat.values, "p_value": pvals, "q_value": qvals,
         "mito_flagged": flagged},
        index=idx)

    gates: dict[str, dict] = {}

    # --- gate 1: measurability ---
    med_mt = float(np.median(mt))
    zero_frac = float((mt <= 0).mean())
    g1_reasons = []
    if n_mito_genes is not None and n_mito_genes < min_genes:
        g1_reasons.append(f"Only {n_mito_genes} mtDNA genes identified, below the minimum of {min_genes}")
    if med_mt < min_median_counts:
        g1_reasons.append(
            f"Median mtDNA count is {med_mt:.0f}, below the minimum of {min_median_counts:.0f}"
            " (per-cell signal is dominated by sampling noise)")
    if zero_frac > max_zero_frac:
        g1_reasons.append(f"{zero_frac:.1%} of cells have zero mtDNA counts, above the {max_zero_frac:.0%} limit")
    if composition_ok is False:
        g1_reasons.append("mtDNA gene composition check failed"
                          " (pctMT does not reflect true mitochondrial content)")
    gates["1_測定可能性"] = {
        "passed": not g1_reasons,
        "detail": " / ".join(g1_reasons) if g1_reasons
                  else f"{n_mito_genes if n_mito_genes is not None else '?'} mtDNA genes, "
                       f"median count {med_mt:.0f}, zero-count fraction {zero_frac:.1%}",
        "median_mt_counts": med_mt, "zero_frac": zero_frac,
        "n_mito_genes": n_mito_genes,
    }

    # --- gate 2: separability ---
    n_flag = int(flagged.sum())
    frac_flag = n_flag / max(n_cells, 1)
    g2_reasons = []
    if n_flag < min_flagged:
        g2_reasons.append(f"FDR<{fdr}: only {n_flag} cells flagged, below the minimum of {min_flagged}"
                          " (essentially no outliers once overdispersion is accounted for)")
    if frac_flag < min_flagged_frac:
        g2_reasons.append(f"Flagged fraction {frac_flag:.3%} is below the minimum of {min_flagged_frac:.1%}")
    if frac_flag > max_flagged_frac:
        g2_reasons.append(f"Flagged fraction {frac_flag:.1%} exceeds the {max_flagged_frac:.0%} limit"
                          " (either the null model is misspecified, or the whole sample is compromised)")
    gates["2_分離可能性"] = {
        "passed": not g2_reasons,
        "detail": " / ".join(g2_reasons) if g2_reasons
                  else f"FDR<{fdr}: {n_flag} cells flagged ({frac_flag:.2%})",
        "n_flagged": n_flag, "frac_flagged": frac_flag,
    }

    # --- gate 3: mechanism fit ---
    #
    # The signature of cell death (membrane rupture) is: nuclear-derived
    # counts fall, mtDNA counts stay roughly constant, and detected-gene
    # complexity falls. A cell that's simply mitochondria-rich instead
    # shows mtDNA counts rising while nuclear counts don't drop. pctMT
    # rises in both cases, so the ratio alone can't distinguish them.
    #
    # The check is run on the flagged set itself rather than as a
    # whole-dataset Spearman correlation -- with only a few percent of
    # cells flagged, a whole-dataset correlation gets diluted by the
    # healthy majority and loses power. It is still computed and kept
    # below as a reference value.
    rho = float("nan")
    if np.std(pct_mt) > 0 and np.std(resid) > 0:
        rho = float(_st.spearmanr(pct_mt, resid).statistic)
    nuc_ratio = mt_ratio = resid_drop_mad = float("nan")
    if n_flag > 0:
        ref: dict[str, tuple[float, float, float, float]] = {}
        for gname, pos in strat.groupby(strat).groups.items():
            sel = idx.get_indexer(pd.Index(pos))
            keep = sel[~flagged[sel]]
            if len(keep) < 5:
                keep = sel
            r_mad = float(_st.median_abs_deviation(resid[keep], scale="normal"))
            ref[gname] = (float(np.median(nuc[keep])),
                          float(np.median(resid[keep])),
                          float(max(np.median(mt[keep]), 0.5)),
                          r_mad if r_mad > 1e-9 else float(np.std(resid[keep]) or 1.0))
        fsel = np.where(flagged)[0]
        gsel = strat.values
        nuc_ratio = float(np.median([nuc[i] / max(ref[gsel[i]][0], 1.0) for i in fsel]))
        mt_ratio = float(np.median([mt[i] / ref[gsel[i]][2] for i in fsel]))
        resid_drop_mad = float(np.median(
            [(resid[i] - ref[gsel[i]][1]) / ref[gsel[i]][3] for i in fsel]))
    g3_reasons = []
    if n_flag == 0:
        g3_reasons.append("No flagged cells to evaluate the mechanism against")
    else:
        if not np.isnan(nuc_ratio) and nuc_ratio > max_nuc_ratio:
            g3_reasons.append(
                f"Flagged cells' nuclear-derived counts are {nuc_ratio:.2f}x the stratum median, "
                f"above the {max_nuc_ratio:.2f}x limit -- no sign of cytoplasmic mRNA loss")
        if not np.isnan(mt_ratio) and mt_ratio > max_mt_ratio:
            g3_reasons.append(
                f"Flagged cells' mtDNA counts are {mt_ratio:.1f}x the stratum median, "
                f"above the {max_mt_ratio:.1f}x limit. Membrane rupture leaves mtDNA roughly "
                "unchanged rather than elevated, so this looks like mitochondria-rich cells, not dead ones")
        if not np.isnan(resid_drop_mad) and resid_drop_mad > -min_resid_drop_mad:
            g3_reasons.append(
                f"Flagged cells' complexity-residual drop is {resid_drop_mad:+.2f} MAD, "
                f"short of {-min_resid_drop_mad:+.2f} MAD -- transcript diversity has not fallen")
    gates["3_機構整合性"] = {
        "passed": not g3_reasons,
        "detail": " / ".join(g3_reasons) if g3_reasons
                  else (f"Flagged cells: nuclear {nuc_ratio:.2f}x, mtDNA {mt_ratio:.1f}x, "
                        f"complexity {resid_drop_mad:+.2f} MAD = consistent with membrane rupture"),
        "nuc_ratio": nuc_ratio, "mt_ratio": mt_ratio,
        "resid_drop_mad": resid_drop_mad,
        "complexity_rho_all_cells": rho,
    }

    # --- gate 4: non-confounding ---
    conf: dict[str, dict] = {}
    g4_reasons = []
    for name, vec in (labels or {}).items():
        v = pd.Series(np.asarray(vec), index=idx).reindex(idx)
        mask = v.notna().values
        if mask.sum() < 20 or n_flag == 0:
            continue
        b = v.fillna(False).astype(bool).values
        a11 = int((flagged & b & mask).sum()); a12 = int((flagged & ~b & mask).sum())
        a21 = int((~flagged & b & mask).sum()); a22 = int((~flagged & ~b & mask).sum())
        if min(a11 + a12, a21 + a22) == 0:
            continue
        orr, pv = _st.fisher_exact([[a11, a12], [a21, a22]], alternative="two-sided")
        rate_in = a11 / max(a11 + a12, 1)
        rate_out = a21 / max(a21 + a22, 1)
        conf[name] = {"odds_ratio": float(orr), "p_value": float(pv),
                      "rate_in_flagged": rate_in, "rate_in_rest": rate_out}
        if pv < confound_p and (orr >= confound_or or (orr > 0 and orr <= 1 / confound_or)):
            g4_reasons.append(
                f"Flagged set is confounded with '{name}' (OR={orr:.2f}, p={pv:.2g}; "
                f"{rate_in:.1%} within flagged vs {rate_out:.1%} in the rest). "
                "Filtering on MT would preferentially discard that group")
    gates["4_非交絡"] = {
        "passed": not g4_reasons,
        "detail": " / ".join(g4_reasons) if g4_reasons
                  else ("no evidence of confounding" if conf else "not evaluated -- no labels available to check"),
        "tests": conf,
    }

    all_passed = all(g["passed"] for g in gates.values())
    lines = [f"MT dead-cell call: recommendation = "
             f"{'remove' if all_passed else 'do not remove'}"]
    for name, g in gates.items():
        lines.append(f"  [{'pass' if g['passed'] else 'fail'}] {name}: {g['detail']}")
    if not all_passed:
        lines.append("  -> Given the failed gate(s) above, removal based on MT lacks support. "
                     "Rely on complexity/total-count QC instead.")
    return {
        "per_cell": per_cell,
        "gates": gates,
        "recommend": "remove" if all_passed else "do_not_remove",
        "reason": "\n".join(lines),
        "fits": fits,
        "stats": {
            "n_cells": n_cells, "fdr": fdr,
            "pct_mt_median": float(np.median(pct_mt)),
            "pct_mt_q99": float(np.percentile(pct_mt, 99)),
            "median_mt_counts": med_mt, "zero_frac": zero_frac,
            "n_flagged": n_flag, "frac_flagged": frac_flag,
            "complexity_rho_all_cells": rho, "nuc_ratio": nuc_ratio,
            "mt_ratio": mt_ratio, "resid_drop_mad": resid_drop_mad,
            "n_over_1pct": int((pct_mt > 1).sum()),
            "n_over_5pct": int((pct_mt > 5).sum()),
            "n_over_10pct": int((pct_mt > 10).sum()),
            "n_flagged_untrimmed_null": n_flag_full,
            "fit_quantile": fit_quantile,
        },
    }


def adaptive_mito_filter(adata, mode: str = "report", *,
                         stratum_key: str | None = "cluster_provisional",
                         label_keys: tuple[str, ...] = ("malignant", "predicted_doublet"),
                         fixed_pct: float | None = None,
                         composition_ok: bool | None = None,
                         n_mito_genes: int | None = None,
                         out_dir: Path | str | None = None,
                         **kwargs) -> dict:
    """Run the MT-based dead-cell judgment on an AnnData and write results to obs.

    mode:
      "off"    no diagnosis at all (adds nothing)
      "report" default. Diagnoses and logs a recommendation but **removes nothing**
      "auto"   sets mito_dead=True only when all 4 gates pass
      "force"  sets mito_dead=True from the adaptive judgment, ignoring the gates (with a warning)
      "fixed"  sets mito_dead=True using a Seurat-style fixed threshold (--mito-filter-pct)

    No mode actually drops cells; the caller filters on obs['mito_dead'].
    """
    if mode == "off":
        log("MT dead-cell judgment: --mito-filter off, skipping diagnosis")
        return {"mode": "off", "recommend": None, "n_flagged": 0}
    if mode not in ("report", "auto", "force", "fixed"):
        raise ValueError(f"Invalid --mito-filter: {mode}")

    counts = adata.layers.get("counts", adata.X)
    is_mt = (np.asarray(adata.var["mt"], dtype=bool) if "mt" in adata.var
             else np.zeros(adata.n_vars, dtype=bool))
    if n_mito_genes is None:
        n_mito_genes = int(is_mt.sum())
    tot = np.asarray(counts.sum(axis=1)).ravel()
    mtc = (np.asarray(counts[:, is_mt].sum(axis=1)).ravel() if is_mt.any()
           else np.zeros(adata.n_obs))
    ng = (np.asarray(adata.obs["n_genes_by_counts"].values, dtype=float)
          if "n_genes_by_counts" in adata.obs
          else np.asarray((counts > 0).sum(axis=1)).ravel().astype(float))

    strat = (adata.obs[stratum_key].values
             if stratum_key and stratum_key in adata.obs else None)
    labels = {}
    for key in label_keys:
        if key in adata.obs:
            col = adata.obs[key]
            if col.dtype == object or str(col.dtype) == "category":
                labels[key] = col.astype(str).isin(("malignant", "True", "true")).values
            else:
                labels[key] = np.asarray(col.values, dtype=bool)

    rep = mito_dead_cell_report(
        tot, mtc, ng, stratum=strat, index=adata.obs_names,
        n_mito_genes=n_mito_genes, labels=labels or None,
        composition_ok=composition_ok, **kwargs)

    for line in rep["reason"].split("\n"):
        (log if rep["recommend"] == "remove" or line.startswith("MT ") else warn)(line)
    st = rep["stats"]
    log(f"  observed: pctMT median {st['pct_mt_median']:.3f}% / 99th pct "
        f"{st['pct_mt_q99']:.2f}% / median mtDNA count {st['median_mt_counts']:.0f}")
    log(f"  at fixed thresholds: pctMT>1% {st['n_over_1pct']} cells / >5% {st['n_over_5pct']} cells"
        f" / >10% {st['n_over_10pct']} cells")

    per = rep["per_cell"]
    adata.obs["pct_counts_mt_adaptive"] = per["pct_mt"].values
    adata.obs["mito_complexity_resid"] = per["complexity_resid"].values
    adata.obs["mito_qvalue"] = per["q_value"].values
    adata.obs["mito_flagged"] = per["mito_flagged"].values

    if mode == "fixed":
        thr = MITO_FIXED_PCT if fixed_pct is None else fixed_pct
        dead = per["pct_mt"].values > thr
        note = f"fixed threshold pctMT>{thr}%"
        warn(f"MT: {note} marks {int(dead.sum())} cells for removal"
             " (Seurat-compatible; applies independently of the diagnosis above)")
    elif mode == "auto":
        if rep["recommend"] == "remove":
            dead = per["mito_flagged"].values
            note = f"adaptive judgment FDR<{rep['stats']['fdr']} (all 4 gates passed)"
            log(f"MT: {note} marks {int(dead.sum())} cells for removal")
        else:
            dead = np.zeros(adata.n_obs, dtype=bool)
            note = "adaptive judgment (no removal -- gate(s) failed)"
            warn("MT: --mito-filter auto but the gates did not pass, so"
                 " no cells are removed. Use --mito-filter force to override")
    elif mode == "force":
        dead = per["mito_flagged"].values
        note = f"adaptive judgment FDR<{rep['stats']['fdr']} (gates ignored)"
        warn(f"MT: --mito-filter force ignores the gates and marks {int(dead.sum())}"
             " cells for removal. Review the failure reasons above carefully")
    else:
        dead = np.zeros(adata.n_obs, dtype=bool)
        note = "diagnosis only (no removal)"
    adata.obs["mito_dead"] = dead

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_tsv_gz(per.assign(mito_dead=dead), out_dir / "mito_filter_percell.tsv.gz",
                     index_label="barcode")
        summary = {"mode": mode, "note": note,
                   "recommend": rep["recommend"], "reason": rep["reason"],
                   "n_marked_dead": int(dead.sum()),
                   "gates": {k: {kk: vv for kk, vv in v.items()}
                             for k, v in rep["gates"].items()},
                   "stats": st, "fits": rep["fits"]}
        (out_dir / "mito_filter_report.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2, default=str),
            encoding="utf-8")
        log(f"  Wrote MT judgment details: {out_dir / 'mito_filter_report.json'}")

    return {"mode": mode, "note": note, "recommend": rep["recommend"],
            "reason": rep["reason"], "n_flagged": int(per["mito_flagged"].sum()),
            "n_marked_dead": int(dead.sum()), "gates": rep["gates"],
            "stats": st, "per_cell": per}


MITO_FIXED_PCT = 5.0   # default for --mito-filter fixed (Seurat's conventional value; weakly justified)


# ===========================================================================
# Export
# ===========================================================================
def export_prep(adata, out_dir: Path, *, orient: str = "cells-genes",
                cluster_key: str = "cluster", gene_pos: pd.DataFrame | None = None,
                params: dict | None = None, write_h5ad: bool = False) -> dict:
    """Write all output files plus the manifest."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    files: dict[str, Path] = {}

    files["barcodes"] = write_lines_gz(adata.obs_names, out_dir / "barcodes.tsv.gz")

    var = pd.DataFrame(index=adata.var_names.astype(str))
    var.index.name = "gene"
    for col in ("mt", "highly_variable", "hvg_for_pca"):
        if col in adata.var:
            var[col] = np.asarray(adata.var[col], dtype=bool)
    if gene_pos is not None:
        for col in ("chromosome", "start", "end"):
            if col in gene_pos.columns:
                var[col] = gene_pos[col].reindex(var.index).values
    files["features"] = write_tsv_gz(var, out_dir / "features.tsv.gz",
                                    index=True, index_label="gene")

    files["matrix_raw"] = write_mtx_gz(adata.layers["counts"],
                                       out_dir / "matrix_raw.mtx.gz",
                                       orient=orient, integer=True)
    lognorm = adata.layers.get("lognorm", adata.X)
    files["matrix_lognorm"] = write_mtx_gz(lognorm,
                                           out_dir / "matrix_lognorm.mtx.gz",
                                           orient=orient, integer=False)

    pca = pd.DataFrame(np.asarray(adata.obsm["X_pca"]), index=adata.obs_names,
                       columns=[f"PC{i + 1}" for i in range(adata.obsm["X_pca"].shape[1])])
    files["pca"] = write_tsv_gz(pca, out_dir / "pca.tsv.gz", index_label="barcode")

    if "X_umap" in adata.obsm:
        um = np.asarray(adata.obsm["X_umap"])
        umap = pd.DataFrame(um, index=adata.obs_names,
                            columns=[f"UMAP{i + 1}" for i in range(um.shape[1])])
        files[f"umap{um.shape[1]}d"] = write_tsv_gz(
            umap, out_dir / f"umap{um.shape[1]}d.tsv.gz", index_label="barcode")

    clu_cols = [c for c in (cluster_key, "cluster_provisional", "cell_type",
                            "sample_id") if c in adata.obs]
    clu = adata.obs[clu_cols].astype(str)
    files["clusters"] = write_tsv_gz(clu, out_dir / "clusters.tsv.gz",
                                     index_label="barcode")

    qc_cols = [c for c in adata.obs.columns
               if c in ("sample_id", "n_counts", "n_genes", "total_counts",
                        "n_genes_by_counts", "pct_counts_mt", "pct_mt",
                        "doublet_score", "predicted_doublet", "final_doublet",
                        "qc_fail")]
    if qc_cols:
        files["qc_metrics"] = write_tsv_gz(adata.obs[qc_cols],
                                           out_dir / "qc_metrics.tsv.gz",
                                           index_label="barcode")

    if write_h5ad:
        h5 = out_dir / "prep.h5ad"
        common._write_h5ad(adata, h5)
        files["h5ad"] = h5

    manifest = {
        "metacellcnv_scanpy_version": __version__,
        "prep_scanpy_version": __version__,   # legacy key (kept for compatibility up to v3.3)
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "matrix_orientation": orient,
        "cluster_key": cluster_key,
        "n_clusters": int(adata.obs[cluster_key].nunique()) if cluster_key in adata.obs else None,
        "params": params or {},
        "versions": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "scanpy": sc.__version__,
            "anndata": __import__("anndata").__version__,
            "scrna_common": getattr(common, "__version__", "unknown"),
        },
        "files": {},
    }
    for key, path in files.items():
        manifest["files"][key] = {
            "path": path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    mpath = out_dir / "prep_manifest.json"
    mpath.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    log(f"Wrote outputs: {out_dir} ({len(files)} files + manifest)")
    for key, path in files.items():
        log(f"  {key:16s} {path.name:24s} {path.stat().st_size / 1e6:8.2f} MB")
    return manifest


# ===========================================================================
# Main processing
# ===========================================================================
def run_prep(args: argparse.Namespace):
    t_start = time.time()
    out_dir = Path(args.out_dir or "prep")
    out_dir.mkdir(parents=True, exist_ok=True)

    paths = args.cellranger_dir or list(getattr(common, "CELLRANGER_DIRS", []) or [])
    if not paths:
        raise SystemExit("Please provide --cellranger-dir")
    sample_ids = args.sample_id if args.sample_id else None

    # --- mtDNA identification: by GTF sequence ID, not by gene name (core of dog support) ---
    gene_pos = None
    mito_genes = None
    mito_chroms = list(args.mito_chromosome or [])
    if args.gtf and Path(args.gtf).exists():
        info = _timed("GTF parsing", common.inspect_gtf, args.gtf,
                      gene_id_attr=args.gtf_gene_id,
                      extra_mito_chromosomes=mito_chroms or None)
        mito_genes = info.get("mito_genes")
        mito_chroms = list(info.get("mito_chromosomes") or mito_chroms)
        log(f"mtDNA sequences: {mito_chroms} / mtDNA genes: {len(mito_genes or [])}")
        try:
            gene_pos = common.parse_gtf_gene_positions(args.gtf, args.gtf_gene_id)
            if args.chromosome_map:
                cmap = common.load_chromosome_map(args.chromosome_map)
                if cmap:
                    gene_pos["chromosome"] = (
                        gene_pos["chromosome"].astype(str).map(lambda c: cmap.get(c, c))
                    )
        except Exception as exc:
            warn(f"Failed to get gene coordinates (features will have no coordinates): {exc}")
    else:
        warn(f"GTF not found ({args.gtf}). mtDNA will be identified from known IDs/prefixes only.")

    # --- 1. load, normalize, HVG, PCA (all cells) ---
    adata = _timed(
        "Load and preprocess", common.load_and_preprocess, paths, sample_ids,
        n_hvg=args.n_hvg, n_pcs=args.n_pcs, mito_genes=mito_genes,
        mito_chromosomes=mito_chroms,
        exclude_mito_from_hvg=not args.keep_mito_in_hvg,
        exclude_ribo_from_hvg=args.exclude_ribosomal_from_hvg,
    )
    common.ensure_sample_id(adata)

    # Minimum-cells-per-gene filter (equivalent to Seurat's CreateSeuratObject(min.cells=...))
    if args.min_cells_per_gene > 0:
        n_before = adata.n_vars
        detected = np.asarray((adata.layers["counts"] > 0).sum(axis=0)).ravel()
        keep = detected >= args.min_cells_per_gene
        adata = adata[:, keep].copy()
        log(f"Gene filter (detected in >= {args.min_cells_per_gene} cells): "
            f"{n_before} -> {adata.n_vars}")

    # --- 2. provisional clustering (for QC stratification) ---
    build_neighbors(adata, n_neighbors=args.n_neighbors, seed=args.seed)
    cluster_cells(adata, resolution=args.cluster_resolution,
                  key_added="cluster_provisional", algo=args.cluster_algo,
                  seed=args.seed)

    # --- 3. mtDNA composition check (feeds gate 1 of the MT judgment) ---
    composition_ok: bool | None = None
    if mito_genes:
        try:
            ref = common.load_mito_reference_profile(args.mito_reference_profile, out_dir)
            comp = common.check_mito_composition(
                adata, reference_profile=ref,
                out_path=out_dir / "mito_gene_profile.csv")
            # check_mito_composition() returns {"status", "problems", ...}.
            # When status is "skipped"/"no_counts" it can't be judged (stays None).
            if isinstance(comp, dict) and "problems" in comp:
                composition_ok = not comp["problems"]
                log(f"  composition check: {'passed' if composition_ok else 'failed'}"
                    f" (status={comp.get('status')})")
        except Exception as exc:
            warn(f"Skipped mtDNA composition check: {exc}")

    # --- 4a. MT-derived dead-cell judgment (diagnosis only, no removal, by default) ---
    # Run after provisional clustering, since the null is fit per stratum.
    mito_res = adaptive_mito_filter(
        adata, mode=args.mito_filter,
        stratum_key="cluster_provisional",
        fixed_pct=args.mito_filter_pct,
        composition_ok=composition_ok,
        n_mito_genes=int(len(mito_genes or [])) or None,
        fdr=args.mito_filter_fdr,
        fit_quantile=args.mito_fit_quantile,
        out_dir=out_dir,
    )

    # --- 4b. MAD-based QC (stratified by sample x provisional cluster) ---
    qc = common.parametric_qc_filter(
        adata, sample_key="sample_id", coarse_cluster_key="cluster_provisional",
        pct_mt_nmads=args.mito_nmads, use_pctmt=not args.no_pctmt_filter,
    )
    n_before = adata.n_obs
    fail = np.asarray(qc["qc_fail"].values, dtype=bool)
    if args.mito_filter != "off" and "mito_dead" in adata.obs:
        dead = np.asarray(adata.obs["mito_dead"].values, dtype=bool)
        if dead.any():
            log(f"Cells marked for removal by MT judgment: {int(dead.sum())}"
                f" (of which {int((dead & fail).sum())} also fail MAD QC)")
        fail = fail | dead
    adata.obs["qc_fail"] = fail
    adata = adata[~fail].copy()
    log(f"Cells after QC: {adata.n_obs} (excluded {n_before - adata.n_obs})")

    # --- 5. doublet detection ---
    if not args.no_doublet:
        try:
            calls = common.detect_doublets_cluster_aware(
                adata, coarse_cluster_key="cluster_provisional",
                expected_doublet_rate=args.expected_doublet_rate,
                max_doublet_rate=args.max_doublet_rate,
            )
            adata.obs["doublet_score"] = calls["doublet_score"].values
            adata.obs["predicted_doublet"] = np.asarray(
                calls["predicted_doublet"].values, dtype=bool)
            n_before = adata.n_obs
            adata = adata[~adata.obs["predicted_doublet"].values].copy()
            log(f"Cells after doublet exclusion: {adata.n_obs} (excluded {n_before - adata.n_obs})")
        except Exception as exc:
            warn(f"Doublet detection failed, skipping: {exc}")
    else:
        log("Skipping doublet detection (--no-doublet)")

    # --- 6. recompute the embedding on the final cell set ---
    if args.no_refit:
        warn("--no-refit: keeping the HVGs and PCA chosen before QC"
             " (matches the existing pipeline's numbers)")
        if "lognorm" not in adata.layers:
            adata.layers["lognorm"] = adata.X.copy()
    else:
        _timed("Recompute on final cell set (normalize/HVG/PCA)", refit_embedding,
               adata, args.n_hvg, args.n_pcs,
               norm_target_sum=args.norm_target_sum,
               hvg_flavor=args.hvg_flavor,
               scale_before_pca=args.scale_before_pca,
               scale_max_value=args.scale_max_value,
               exclude_mito_from_hvg=not args.keep_mito_in_hvg,
               exclude_ribo_from_hvg=args.exclude_ribosomal_from_hvg)

    # --- 7. final clustering and UMAP ---
    build_neighbors(adata, n_neighbors=args.n_neighbors, seed=args.seed)
    cluster_cells(adata, resolution=args.cluster_resolution,
                  key_added="cluster", algo=args.cluster_algo, seed=args.seed)
    _timed(f"UMAP {args.umap_components}D", embed_umap, adata,
           n_components=args.umap_components,
           n_neighbors=args.umap_neighbors, min_dist=args.umap_min_dist,
           seed=args.seed)
    if args.umap_neighbors is not None:
        # If the graph was rebuilt for UMAP, restore the clustering graph
        build_neighbors(adata, n_neighbors=args.n_neighbors, seed=args.seed)

    # --- 8. export ---
    params = {k: (str(v) if isinstance(v, Path) else v)
              for k, v in vars(args).items()}
    params["mito_chromosomes"] = mito_chroms
    params["n_mito_genes"] = int(len(mito_genes or []))
    params["mito_filter_result"] = {
        "mode": mito_res.get("mode"),
        "recommend": mito_res.get("recommend"),
        "n_flagged": mito_res.get("n_flagged"),
        "n_marked_dead": mito_res.get("n_marked_dead"),
    }
    manifest = export_prep(adata, out_dir, orient=args.orient,
                           cluster_key="cluster", gene_pos=gene_pos,
                           params=params, write_h5ad=args.write_h5ad)

    log(f"Done: {time.time() - t_start:.1f}s")
    return adata, manifest


def run_containment(args: argparse.Namespace) -> dict:
    """Measure containment of an existing metacell assignment against clusters only."""
    apath = Path(args.containment)
    assign = pd.read_csv(apath, index_col=0)
    if args.clusters is None:
        raise SystemExit("--containment also requires --clusters")
    cpath = Path(args.clusters)
    clu_df = pd.read_csv(cpath, sep="\t" if ".tsv" in cpath.name else ",", index_col=0)
    col = args.cluster_column or ("cluster" if "cluster" in clu_df.columns
                                  else clu_df.columns[0])
    res = metacell_cluster_containment(assign, clu_df[col])
    report_containment(res)
    if args.out_dir:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_tsv_gz(res["per_metacell"], out_dir / "metacell_containment.tsv.gz",
                     index_label="metacell")
        (out_dir / "metacell_containment.json").write_text(
            json.dumps(res["summary"], ensure_ascii=False, indent=2), encoding="utf-8")
        log(f"Wrote containment results: {out_dir}")
    return res


# ===========================================================================
# CLI
# ===========================================================================
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Seurat-free unified preprocessing (scanpy only, metacellcnv_scanpy v%s)" % __version__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # --- arguments with the same names as the main pipeline (so command lines can be reused) ---
    p.add_argument("--cellranger-dir", action="append", default=None,
                   help="CellRanger output directory (repeatable)")
    p.add_argument("--sample-id", action="append", default=None,
                   help="sample_id corresponding to each --cellranger-dir")
    p.add_argument("--gtf", default=getattr(common, "GTF_PATH", "genes.gtf"),
                   help="Path to the gene-coordinate annotation GTF")
    p.add_argument("--gtf-gene-id", default=getattr(common, "GTF_GENE_ID_ATTR", "gene"),
                   help="GTF gene-name attribute key ('auto' to autodetect)")
    p.add_argument("--chromosome-map", default=getattr(common, "CHROMOSOME_MAP_PATH", None),
                   help="NCBI assembly report, or a 2-column chromosome-name mapping table")
    p.add_argument("--mito-chromosome", action="append", default=None,
                   help="Sequence ID of the mitochondrial genome (repeatable). "
                        "NC_002008.4 for dog CanFam6 / Dog10K_Boxer_Tasha")
    p.add_argument("--mito-reference-profile", default=None,
                   help="[Optional input] Reference CSV of mtDNA gene composition")
    p.add_argument("--mito-nmads", type=float, default=3.5,
                   help="Number of MADs used to flag pctMT outliers")
    p.add_argument("--no-pctmt-filter", action="store_true",
                   help="Disable cell exclusion based on pctMT")
    p.add_argument("--keep-mito-in-hvg", action="store_true",
                   help="Keep mtDNA genes in HVG/PCA (excluded by default)")
    p.add_argument("--exclude-ribosomal-from-hvg", action="store_true",
                   help="Also exclude cytoplasmic ribosomal genes (RPL*/RPS*) from HVG")
    p.add_argument("--n-hvg", type=int, default=getattr(common, "N_HVG", 1500))
    p.add_argument("--n-pcs", type=int, default=getattr(common, "N_PCS", 50))
    p.add_argument("--expected-doublet-rate", type=float, default=0.05)
    p.add_argument("--max-doublet-rate", type=float, default=0.25)
    p.add_argument("--out-dir", default=None,
                   help="Output directory for results (default 'prep'). With --containment, "
                        "results are written only if this is given")

    # --- options specific to this module ---
    p.add_argument("--n-neighbors", type=int, default=15,
                   help="k for the neighbor graph (used for clustering)")
    p.add_argument("--cluster-resolution", type=float, default=0.5,
                   help="Clustering resolution")
    p.add_argument("--cluster-algo", choices=["leiden", "louvain"], default="leiden",
                   help="Clustering method (Seurat's default SLM has no scanpy equivalent)")
    p.add_argument("--umap-components", type=int, default=3,
                   help="Number of UMAP dimensions")
    p.add_argument("--umap-neighbors", type=int, default=None,
                   help="k for rebuilding the neighbor graph for UMAP. If unset, reuses the clustering graph")
    p.add_argument("--umap-min-dist", type=float, default=0.5)
    p.add_argument("--norm-target-sum", type=float, default=None,
                   help="Normalization target sum. Unset uses the scanpy default (median counts). "
                        "Use 10000 to match Seurat's LogNormalize")
    p.add_argument("--hvg-flavor", choices=["seurat", "seurat_v3", "cell_ranger"],
                   default="seurat",
                   help="HVG method. seurat_v3 is the equivalent of Seurat's vst (needs scikit-misc)")
    p.add_argument("--scale-before-pca", action="store_true",
                   help="Z-score before PCA (equivalent to Seurat's ScaleData). "
                        "Off by default, matching the existing pipeline")
    p.add_argument("--scale-max-value", type=float, default=SEURAT_SCALE_MAX,
                   help="Clip value used with --scale-before-pca")
    p.add_argument("--mito-filter",
                   choices=["off", "report", "auto", "force", "fixed"],
                   default="report",
                   help="MT-derived dead-cell removal. "
                        "report=default; diagnoses and recommends but removes nothing / "
                        "auto=removes only if all 4 gates (measurability, separability, "
                        "mechanism fit, non-confounding) pass / force=removes, ignoring gates / "
                        "fixed=removes using a Seurat-style fixed threshold / off=no diagnosis")
    p.add_argument("--mito-filter-pct", type=float, default=MITO_FIXED_PCT,
                   help="pctMT threshold for --mito-filter fixed. Not recommended, since the "
                        "baseline varies by orders of magnitude across species/reference/alignment")
    p.add_argument("--mito-filter-fdr", type=float, default=MITO_FDR,
                   help="FDR for the adaptive judgment (beta-binomial null)")
    p.add_argument("--mito-fit-quantile", type=float, default=MITO_FIT_QUANTILE,
                   help="Upper quantile trimmed when fitting the null. 1.0 fits on all cells, "
                        "but a small number of dead cells can widen their own null and "
                        "hide themselves (breaking the moment estimator)")
    p.add_argument("--min-cells-per-gene", type=int, default=0,
                   help="Keep only genes detected in at least this many cells (0 to disable)")
    p.add_argument("--orient", choices=["cells-genes", "genes-cells"],
                   default="cells-genes",
                   help="Orientation of the mtx file. genes-cells matches 10x / Seurat's ReadMtx")
    p.add_argument("--no-refit", action="store_true",
                   help="Don't recompute normalization/HVG/PCA after QC (matches the existing pipeline)")
    p.add_argument("--no-doublet", action="store_true",
                   help="Skip doublet detection")
    p.add_argument("--write-h5ad", action="store_true",
                   help="Also write prep.h5ad (for downstream code that uses AnnData directly)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--seurat-compat", action="store_true",
                   help="Move toward Seurat defaults: --norm-target-sum 10000, "
                        "--hvg-flavor seurat_v3, --scale-before-pca, "
                        "--umap-neighbors 30, --umap-min-dist 0.3, "
                        "--min-cells-per-gene 3. Clustering stays Leiden "
                        "(Seurat's default SLM has no scanpy equivalent, and Louvain won't match either)")

    # --- mode that only measures containment ---
    p.add_argument("--containment", default=None,
                   help="Measure cluster containment for an existing metacell assignment CSV"
                        " (cell_to_metacell.csv) and exit")
    p.add_argument("--clusters", default=None,
                   help="Cluster file used with --containment (clusters.tsv.gz)")
    p.add_argument("--cluster-column", default=None,
                   help="Column name to use within --clusters (default 'cluster')")

    args = p.parse_args(argv)
    if args.seurat_compat:
        log("--seurat-compat: moving toward Seurat defaults (numbers will differ from the existing pipeline)")
        args.norm_target_sum = SEURAT_SCALE_FACTOR
        args.hvg_flavor = "seurat_v3"
        args.scale_before_pca = True
        args.umap_neighbors = SEURAT_UMAP_NEIGHBORS
        args.umap_min_dist = SEURAT_UMAP_MIN_DIST
        args.min_cells_per_gene = SEURAT_MIN_CELLS_PER_GENE
        warn("Clustering stays Leiden. Seurat's default SLM has no scanpy equivalent, "
             "and Louvain (--cluster-algo louvain) won't match cluster counts either.")
    return args


def main(argv: list[str] | None = None) -> int:
    warnings.filterwarnings("ignore", category=FutureWarning)
    # Performance warning from scanpy/scipy internals when assigning into sparse
    # matrices. Not something we can fix here and it doesn't affect results, so
    # suppress it to keep the log readable.
    try:
        from scipy.sparse import SparseEfficiencyWarning
        warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
    except Exception:
        pass
    args = parse_args(argv)
    if args.containment:
        run_containment(args)
        return 0
    run_prep(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
