#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""metacell_annotation.py -- Classify metacells as "confidently typed", "mixed", or "unclassified".

Why this module exists
=======================
Legacy metacell labels (e.g. `cell_type` values like `Myeloid_3`) suffered
from three problems:

(1) **Panel scores were driven by genes with no discriminative power.**
    Detection rates for several "typing" genes were essentially the same in
    tumor and non-tumor cells, so moving the panel score didn't actually
    distinguish cell type.

(2) **`Admixture` claimed to measure something it didn't.** When mixture was
    actually measured at the cell level, most metacells were effectively
    pure and only a small fraction were truly mixed -- the old `Admixture`
    label conflated "measured as mixed" with "insufficient evidence to call
    a type," giving two different situations the same name.

(3) **Thresholds were ad hoc** ("more than 2x", "minority side >= 10%",
    "eyeball the valley").

This module fixes all three: discriminative power is checked with a
statistical test, mixture is measured (not assumed), and thresholds are
derived from the data. When a call can't be made, classification is refused
and the reason is recorded. See README.md for validation results.

Pipeline
========
    (1) Per-cell karyotype projection -> threshold from a 2-component
        mixture model (refused if the distribution isn't bimodal)
    (2) Per-metacell tumor-cell fraction -> "mixed" is measured via a
        binomial test
    (3) Non-mixed metacells are typed from aggregated panel scores; the
        threshold is set via FDR over a null distribution built by
        shuffling cell->metacell assignment
    (4) Metacells that can't be typed are labeled Unclassified, by reason

Output labels
=============
    typed            : Tumor:<type> / <type>            a type could be called
    measured-mixture : Mixed:tumor+normal                mixture was measured
    unclassified      : Unclassified:low-depth            too few callable cells
                        Unclassified:no-dominant-type    cells present but no type stands out
                        Unclassified:multiple-types      two or more types significant at once
                        Tumor:Unclassified               tumor, but phenotype undetermined

Design choice: metacells are not split
=======================================
Per-cell karyotype projection is accurate enough that splitting cells into
tumor/normal before metacell construction is technically possible. This
module does not do that, because:

  - Most metacells are already CNV-pure; splitting would only fix a small
    minority of them.
  - The main cause of "unclassified" is insufficient depth, not mixture --
    splitting a metacell doesn't add information per cell.
  - Splitting removes candidates from the neighbor graph, shrinking exactly
    the resource that's already scarce: the number of similar cells grouped
    into one metacell.

So mixture is not designed away -- it is measured and reported.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

__version__ = "1.0"

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import scrna_common as common  # noqa: E402

log, warn = common.log, common.warn


# ===========================================================================
# Cell-type panels
# ===========================================================================
# Only genes whose specificity was validated against a depth-matched control
# (candidate cells vs. a control matched on number of detected genes) in the
# reference species are used for classification. The Epithelial panel was
# instead validated via its association with karyotype (tumor side vs.
# normal side).
def panels_from_table(spec: str | None = "dog") -> tuple[dict, dict]:
    """Build (label panels, reference panels) from the marker CSV.

    The two dicts below used to be hardcoded directly (species-specific).
    They were externalized so the species can be swapped via `--markers`,
    but the hardcoded dicts remain here as defaults.
    """
    tbl = common.load_marker_table(spec)
    return tbl.label_panels, tbl.reported_panels


VALIDATED_PANELS: dict[str, list[str]] = {
    "Macrophage": ["CD68", "AIF1", "TYROBP", "CD14", "FCER1G", "LYZ",
                   "C1QA", "C1QB", "C1QC", "CTSS"],
    "Endothelial": ["PECAM1", "CDH5", "VWF", "KDR", "CLDN5", "EGFL7",
                    "TEK", "ERG", "FLT1", "ESAM", "RAMP2", "PLVAP"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "DCN",
                   "LUM", "POSTN", "SPARC", "PDGFRB", "THY1", "FAP"],
    "Epithelial": ["KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "EPCAM",
                   "CLDN4", "SFN"],
}
# Panels that are aggregated for reference but not used for classification --
# specificity has not been validated for these: canonical marker coverage is
# sparse in the reference species, and per-cell detection counts are too low.
REPORTED_PANELS: dict[str, list[str]] = {
    "TNK": ["CD3D", "CD3E", "CD3G", "CD2", "LCK", "ITK", "THEMIS", "SKAP1",
            "GZMA", "NKG7", "KLRD1", "IL7R"],
    "BPlasma": ["CD19", "PAX5", "BANK1", "BLK", "CD22", "FCRL1",
                "JCHAIN", "XBP1", "DERL3", "SDC1"],
    "SmoothMuscle": ["ACTA2", "TAGLN", "MYH11", "RGS5", "NOTCH3", "CNN1", "DES"],
}

# Boundary between "low depth" and "no dominant type" among unclassified
# metacells: below this many callable cells there isn't enough basis to
# discuss type at all.
MIN_CALLED_CELLS = 3
# Number of positive genes required for a per-cell panel call. This can
# instead be set via FDR against a Chung-Lu null (the `fdr` argument of
# call_cell_types); the default is this fixed, previously validated value.
MIN_PANEL_POSITIVE = 3
# FDR used for typing and for mixture calls.
PANEL_FDR = 0.01
MIXTURE_FDR = 0.05
# Floor on the "per-cell karyotype call error rate" used in the mixture
# null. It is estimated from the mixture model's overlap, but is floored
# here so the estimate can't become unrealistically optimistic.
MIN_MISCALL_RATE = 0.02


# ===========================================================================
# Panel availability and per-cell classification
# ===========================================================================
def panel_gene_availability(var_names, panels: dict[str, list[str]]) -> pd.DataFrame:
    """Tabulate, for each panel, whether each marker gene is present in the reference.

    In non-model reference species, canonical markers can be annotated under
    a locus ID instead of a gene symbol, or be missing entirely. Checking
    this before classification matters -- otherwise a panel's score can end
    up driven purely by whichever non-specific genes happen to remain.
    """
    have = {str(g).upper() for g in var_names}
    rows = []
    for name, genes in panels.items():
        present = [g for g in genes if g.upper() in have]
        rows.append({"panel": name, "n_genes": len(genes),
                     "n_present": len(present),
                     "missing": ",".join(g for g in genes if g.upper() not in have)})
    return pd.DataFrame(rows).set_index("panel")


def _resolve(adata, genes: list[str]) -> list[str]:
    idx = {str(g).upper(): g for g in adata.var_names}
    return [idx[g.upper()] for g in genes if g.upper() in idx]


def _detected(adata, genes: list[str]) -> np.ndarray:
    """Return, per cell, the number of the given genes that are detected."""
    gg = _resolve(adata, genes)
    if not gg:
        return np.zeros(adata.n_obs, dtype=int)
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    sub = X[:, [adata.var_names.get_loc(g) for g in gg]]
    return np.asarray((sub > 0).sum(1)).ravel().astype(int)


def panel_positive_counts(adata, panels: dict[str, list[str]]) -> pd.DataFrame:
    """Cell x panel matrix of positive-gene counts."""
    return pd.DataFrame({p: _detected(adata, g) for p, g in panels.items()},
                        index=adata.obs_names)


def _n_genes(adata) -> np.ndarray:
    if "n_genes_by_counts" in adata.obs:
        return np.asarray(adata.obs["n_genes_by_counts"].values, dtype=float)
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    return np.asarray((X > 0).sum(1)).ravel().astype(float)


def chung_lu_expected(adata, genes: list[str]) -> np.ndarray:
    """Expected number of panel genes detected in cell i, under a Chung-Lu null.

    P(gene g detected in cell i) ~= r_i * c_g / T
      r_i = number of genes detected in cell i, c_g = number of cells gene g
      is detected in, T = sum of r_i
    This preserves per-cell depth (r_i) and per-gene detectability (c_g), so
    the null captures the effect of "a deep cell tends to be positive for
    almost anything."
    """
    gg = _resolve(adata, genes)
    if not gg:
        return np.zeros(adata.n_obs)
    r = _n_genes(adata)
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    cols = [adata.var_names.get_loc(g) for g in gg]
    c = np.asarray((X[:, cols] > 0).sum(0)).ravel().astype(float)
    return np.clip(np.outer(r, c) / max(r.sum(), 1.0), 0, 1).sum(1)


def call_cell_types(adata, panels: dict[str, list[str]] | None = None, *,
                    min_positive: int = MIN_PANEL_POSITIVE,
                    fdr: float | None = None) -> pd.DataFrame:
    """Call at most one type per cell (positive on 2+ panels -> doublet candidate).

    If `fdr` is given, the per-panel positive-gene-count threshold is chosen
    from an FDR against a Chung-Lu null, per panel. If not given, the fixed
    `min_positive` value is used.

    Important limitation: this per-cell call is not highly reproducible.
    At typical detected-gene depths, which 2-3 genes of a 10-gene panel
    happen to be detected is close to chance. **Use this for population-
    level aggregation, not as a label to trust for an individual cell.**
    Metacell classification relies mainly on the aggregated panel score.
    """
    panels = panels or VALIDATED_PANELS
    K = panel_positive_counts(adata, panels)
    thr = {}
    for p, genes in panels.items():
        if fdr is None:
            thr[p] = min_positive
            continue
        exp = chung_lu_expected(adata, genes)
        pv = stats.poisson.sf(K[p].values - 1, np.maximum(exp, 1e-9))
        q = _bh(pv)
        sig = K[p].values[q < fdr]
        thr[p] = int(sig.min()) if len(sig) else min_positive
    pos = np.column_stack([K[p].values >= thr[p] for p in panels])
    npos = pos.sum(1)
    names = np.array(list(panels))
    called = np.where(npos == 1, names[np.argmax(pos, axis=1)],
                      np.where(npos >= 2, "doublet", "none"))
    out = K.copy()
    out["n_panels_positive"] = npos
    out["called"] = called
    out.attrs["thresholds"] = thr
    return out


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    m = len(p)
    if m == 0:
        return p
    o = np.argsort(p)
    q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.clip(q, 0.0, 1.0)


# ===========================================================================
# Metacell-level panel scores and their null distribution
# ===========================================================================
def metacell_panel_scores(adata, metacell_key: str,
                          panels: dict[str, list[str]] | None = None,
                          return_ratios: bool = False):
    """Per metacell, the log2 median of obs/exp across each panel's genes.

    Using the **median** rather than the sum is the key point. A sum is
    dominated by a single broadly-detected gene (a panel can end up
    "significant" purely on the strength of one gene like that, even when
    most of the panel is not actually elevated). A median only rises once a
    majority of the panel's genes move.

    Also returns the number of genes whose ratio exceeds 2x, since accepting
    a type requires several genes to be elevated together.
    """
    panels = panels or VALIDATED_PANELS
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    mc = pd.Series(np.asarray(adata.obs[metacell_key].values, dtype=object),
                   index=adata.obs_names).astype(str)
    r = _n_genes(adata)
    T = r.sum()
    codes, uniq = pd.factorize(mc.values)
    n_mc = len(uniq)
    sum_r = np.bincount(codes, weights=r, minlength=n_mc)
    out = {}; ratios = {}
    for p, genes in panels.items():
        gg = _resolve(adata, genes)
        if len(gg) < 3:
            out[p] = np.full(n_mc, np.nan)
            out[p + "_nsig"] = np.zeros(n_mc, dtype=int)
            ratios[p] = np.full((n_mc, 0), np.nan)
            continue
        cols = [adata.var_names.get_loc(g) for g in gg]
        D = (X[:, cols] > 0)
        Dd = np.asarray(D.todense()) if hasattr(D, "todense") else np.asarray(D)
        c = Dd.sum(0).astype(float)                       # per-gene detected-cell count
        obs = np.zeros((n_mc, len(gg)))
        for j in range(len(gg)):
            obs[:, j] = np.bincount(codes, weights=Dd[:, j].astype(float),
                                    minlength=n_mc)
        exp = np.outer(sum_r, c) / max(T, 1.0)
        ratio = np.log2((obs + 0.5) / (exp + 0.5))
        out[p] = np.median(ratio, axis=1)
        out[p + "_nsig"] = (ratio > 1.0).sum(1)     # reference: number of genes over 2x
        ratios[p] = ratio
    S = pd.DataFrame(out, index=pd.Index(uniq, name=metacell_key))
    if return_ratios:
        return S, {p: pd.DataFrame(v, index=S.index,
                                   columns=[g for g in _resolve(adata, panels[p])])
                   for p, v in ratios.items()}
    return S


def calibrate_panel_thresholds(adata, metacell_key: str,
                               panels: dict[str, list[str]] | None = None, *,
                               n_perm: int = 30, fdr: float = PANEL_FDR,
                               depth_bins: int = 5,
                               seed: int = 0) -> dict:
    """Set panel-score thresholds from a null built by shuffling cell->metacell assignment.

    Shuffling is done **within depth bins**. Shuffling across depth bins
    would erase the real structure of "deep cells cluster into the same
    metacell," making the null too permissive.

    Returns, per panel, the threshold (the 1-fdr quantile of the null) and a
    summary of the null distribution. Because it is derived from the data
    rather than a fixed "2x" cutoff, it automatically adapts to each
    sample's depth distribution and panel gene availability.
    """
    panels = panels or VALIDATED_PANELS
    rng = np.random.default_rng(seed)
    ng = _n_genes(adata)
    q = np.quantile(ng, np.linspace(0, 1, depth_bins + 1)[1:-1])
    strat = np.digitize(ng, q)
    orig = np.asarray(adata.obs[metacell_key].values, dtype=object).astype(str)
    null = {p: [] for p in panels}
    tmp_key = "__perm_mc__"
    for _ in range(n_perm):
        perm = orig.copy()
        for s in np.unique(strat):
            m = strat == s
            perm[m] = rng.permutation(perm[m])
        adata.obs[tmp_key] = perm
        S = metacell_panel_scores(adata, tmp_key, panels)
        for p in panels:
            v = S[p].values
            null[p].append(v[np.isfinite(v)])
    if tmp_key in adata.obs:
        del adata.obs[tmp_key]
    res = {}
    for p in panels:
        v = np.concatenate(null[p]) if null[p] else np.array([0.0])
        res[p] = {"threshold": float(np.quantile(v, 1.0 - fdr)),
                  "null_median": float(np.median(v)),
                  "null_sd": float(np.std(v)),
                  "n_null": int(len(v))}
    return res


# ===========================================================================
# Karyotype and per-cell tumor calling
# ===========================================================================
def chromosome_karyotype(adata, groups, group_a, group_b, *,
                         chrom_key: str = "chromosome",
                         min_counts: int = 200) -> pd.Series:
    """Build a "karyotype profile": per-chromosome mean of per-gene log2 ratio between two groups.

    Computed directly from raw counts rather than from infercnvpy output.
    infercnvpy expresses values as a log-ratio to the reference group, which
    forces the reference group's own value to 0 by construction -- masking
    any tumor content that may actually be present within the reference
    group itself.
    """
    if chrom_key not in adata.var:
        raise ValueError(f"var['{chrom_key}'] is missing."
                         " Add genomic coordinates first (e.g. via scrna_common.add_genomic_positions)")
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    g = pd.Series(np.asarray(groups), index=adata.obs_names).astype(str)
    ia = (g == str(group_a)).values
    ib = (g == str(group_b)).values
    if ia.sum() < 10 or ib.sum() < 10:
        raise ValueError(f"Group sizes are too small: {group_a}={ia.sum()} / {group_b}={ib.sum()}")
    sa = np.asarray(X[ia].sum(0)).ravel().astype(float)
    sb = np.asarray(X[ib].sum(0)).ravel().astype(float)
    ch = adata.var[chrom_key].astype(str).values
    ok = (ch != "nan") & (ch != "NA") & ((sa + sb) >= min_counts)
    if ok.sum() < 200:
        raise ValueError(f"Genes with coordinates and >= {min_counts} counts: only {ok.sum()} found")
    lr = np.log2(((sa[ok] + 1) / (sa.sum() + ok.sum())) /
                 ((sb[ok] + 1) / (sb.sum() + ok.sum())))
    prof = pd.Series(lr).groupby(ch[ok]).mean()
    return prof - prof.mean()


def cell_karyotype_projection(adata, karyotype: pd.Series, *,
                              chrom_key: str = "chromosome") -> pd.Series:
    """Per-cell projection of chromosome composition onto the karyotype direction.

    Each cell's chromosome composition is expressed as a log2 ratio to the
    cohort mean, then projected onto the karyotype vector. A substantial
    share of the projection's variance is explained by expression program
    and sequencing depth rather than CNV, which is why **this single 1-D
    projection is used instead of adding one dimension per chromosome**.
    """
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    ch = adata.var[chrom_key].astype(str).values
    keep = pd.Index(karyotype.index.astype(str))
    W = np.zeros((adata.n_obs, len(keep)))
    for j, c in enumerate(keep):
        cols = np.where(ch == c)[0]
        if len(cols) == 0:
            continue
        W[:, j] = np.asarray(X[:, cols].sum(1)).ravel()
    tot = np.maximum(W.sum(1, keepdims=True), 1.0)
    frac = W / tot
    gl = W.sum(0) / max(W.sum(), 1.0)
    L = np.log2((frac + 1e-9) / (gl + 1e-9))
    L = L - L.mean(1, keepdims=True)
    k = karyotype.reindex(keep).values.astype(float)
    k = k - k.mean()
    return pd.Series((L @ k) / max(k @ k, 1e-12), index=adata.obs_names)


def fit_bimodal_threshold(x, *, min_weight: float = 0.10,
                          min_separation: float = 2.0, seed: int = 0) -> dict:
    """Estimate a threshold from a 2-component Gaussian mixture; refuse if not bimodal.

    This replaces eyeballing the valley between two modes. If the returned
    dict's `usable` is False, it means a tumor/normal threshold cannot be
    determined from this data, and the caller should give up on mixture
    calling.

    `min_separation` is the minimum allowed distance between component
    means, in pooled SDs.
    """
    from sklearn.mixture import GaussianMixture
    v = np.asarray(x, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 100:
        return {"usable": False, "reason": f"only {len(v)} samples, too few"}
    gm = GaussianMixture(2, n_init=5, random_state=seed).fit(v.reshape(-1, 1))
    mu = gm.means_.ravel()
    sd = np.sqrt(gm.covariances_.ravel())
    w = gm.weights_.ravel()
    o = np.argsort(mu)
    mu, sd, w = mu[o], sd[o], w[o]
    pooled = np.sqrt((sd[0] ** 2 + sd[1] ** 2) / 2)
    sep = (mu[1] - mu[0]) / max(pooled, 1e-9)
    reasons = []
    if w.min() < min_weight:
        reasons.append(f"smaller component weight is {w.min():.3f}, below {min_weight}")
    if sep < min_separation:
        reasons.append(f"component separation is {sep:.2f} SD, below {min_separation} SD")
    grid = np.linspace(mu[0], mu[1], 2001)
    post = gm.predict_proba(grid.reshape(-1, 1))[:, o[1]]
    thr = float(grid[np.argmin(np.abs(post - 0.5))])
    # miscall rate: probability each component falls on the wrong side of the threshold
    err_lo = float(1 - stats.norm.cdf(thr, mu[0], max(sd[0], 1e-9)))
    err_hi = float(stats.norm.cdf(thr, mu[1], max(sd[1], 1e-9)))
    return {"usable": not reasons,
            "reason": " / ".join(reasons) if reasons else
                      f"components separated by {sep:.2f} SD (weights {w[0]:.2f}/{w[1]:.2f})",
            "threshold": thr, "separation": float(sep),
            "means": mu.tolist(), "sds": sd.tolist(), "weights": w.tolist(),
            "miscall_low": err_lo, "miscall_high": err_hi,
            "miscall_rate": float(max(err_lo, err_hi))}


# ===========================================================================
# Metacell annotation
# ===========================================================================
def _mixture_pvalues(k: np.ndarray, n: np.ndarray, eps: float) -> tuple[np.ndarray, np.ndarray]:
    """Binomial-test p-values against two nulls: "all normal" and "all tumor".

    k = number of cells called tumor in the metacell, n = number of cells,
    eps = per-cell call error rate. Even a purely normal metacell has k>0
    due to miscalls, so the null uses p=eps rather than p=0. A metacell is
    called "mixed" only when both nulls are rejected. This handles metacell
    size correctly, unlike a fixed "minority side >= 10%" rule (the same
    fraction is weaker evidence in a 20-cell metacell than in a 150-cell
    one).
    """
    p_not_normal = stats.binom.sf(k - 1, n, eps)          # probability of k or more
    p_not_tumor = stats.binom.cdf(k, n, 1.0 - eps)        # probability of k or fewer
    return p_not_normal, p_not_tumor


def annotate_metacells(adata, metacell_key: str, *,
                       karyotype: pd.Series | None = None,
                       tumor_groups: tuple[str, str] | None = None,
                       group_key: str | None = None,
                       chrom_key: str = "chromosome",
                       panels: dict[str, list[str]] | None = None,
                       reported_panels: dict[str, list[str]] | None = None,
                       markers: str | None = None,
                       panel_fdr: float = PANEL_FDR,
                       mixture_fdr: float = MIXTURE_FDR,
                       min_called_cells: int = MIN_CALLED_CELLS,
                       n_perm: int = 30, seed: int = 0,
                       out_dir: Path | str | None = None) -> pd.DataFrame:
    """Classify metacells into typed / measured-mixture / unclassified.

    If `karyotype` is not given, it is computed from `group_key` and
    `tumor_groups` (e.g. group_key="putative_malignant",
    tumor_groups=("malignant", "normal")). If the karyotype cannot be
    computed, mixture calling is skipped and only typing is performed.
    """
    if markers is not None and panels is None and reported_panels is None:
        panels, reported = panels_from_table(markers)
    else:
        panels = panels or VALIDATED_PANELS
        reported = reported_panels if reported_panels is not None else REPORTED_PANELS

    avail = panel_gene_availability(adata.var_names, {**panels, **reported})
    log("Panel gene availability:")
    for p, r in avail.iterrows():
        tag = "used for typing" if p in panels else "reference only"
        miss = f" missing: {r.missing}" if r.missing else ""
        log(f"  {p:<14} {r.n_present}/{r.n_genes} present ({tag}){miss}")
    usable = [p for p in panels if avail.loc[p, "n_present"] >= 3]
    dropped = [p for p in panels if p not in usable]
    if dropped:
        warn(f"Dropping from classification (fewer than 3 genes present in reference): {dropped}")
    panels = {p: panels[p] for p in usable}
    if not panels:
        raise ValueError("No panels remain usable for classification")

    # --- Per-cell typing (used for population aggregation and doublet tracking) ---
    cells = call_cell_types(adata, panels, fdr=None)
    mc = pd.Series(np.asarray(adata.obs[metacell_key].values, dtype=object),
                   index=adata.obs_names).astype(str)
    ng = _n_genes(adata)

    # --- Metacell-level panel scores, and thresholds from the shuffle null ---
    S, RATIOS = metacell_panel_scores(adata, metacell_key, {**panels, **reported},
                                      return_ratios=True)
    thr = calibrate_panel_thresholds(adata, metacell_key, panels,
                                     n_perm=n_perm, fdr=panel_fdr, seed=seed)
    log(f"Panel score thresholds (cell->metacell shuffle x{n_perm}, FDR {panel_fdr}):")
    for p in panels:
        t = thr[p]
        log(f"  {p:<14} threshold {t['threshold']:+.3f}"
            f" (null median {t['null_median']:+.3f} +/- {t['null_sd']:.3f})")

    # --- Karyotype and per-cell tumor calling ---
    idx_pre = S.index
    proj = None
    fit = {"usable": False, "reason": "no karyotype provided"}
    if karyotype is None and group_key is not None and tumor_groups is not None:
        try:
            karyotype = chromosome_karyotype(adata, adata.obs[group_key].values,
                                             tumor_groups[0], tumor_groups[1],
                                             chrom_key=chrom_key)
            log(f"Computed karyotype profile ({len(karyotype)} chromosomes, "
                f"amplitude {karyotype.max()-karyotype.min():.3f})")
        except Exception as exc:
            warn(f"Failed to compute karyotype: {exc}")
    mc_proj = None
    if karyotype is not None:
        try:
            proj = cell_karyotype_projection(adata, karyotype, chrom_key=chrom_key)
            # Bimodality is assessed on the **metacell mean**, not the
            # per-cell projection: per-cell values have high variance, and
            # the mixture model can latch onto tails that aren't actually
            # tumor/normal. Averaging over the cells in a metacell makes the
            # two modes come out cleanly.
            mc_proj = proj.groupby(mc.values).mean().reindex(idx_pre)
            fit = fit_bimodal_threshold(mc_proj.dropna().values, seed=seed)
            (log if fit["usable"] else warn)(
                f"Bimodality of karyotype projection (metacell mean): {fit['reason']}")
            if fit["usable"]:
                # The per-cell threshold is the midpoint between component
                # means. The mean is unchanged by aggregation, so the
                # metacell-level estimate can be reused directly (only the
                # variance differs).
                mu_n, mu_t = fit["means"]
                cell_thr = 0.5 * (mu_n + mu_t)
                lo = proj.values[proj.values <= cell_thr]
                hi = proj.values[proj.values > cell_thr]
                sd_n = float(np.std(lo)) if len(lo) > 20 else 1.0
                sd_t = float(np.std(hi)) if len(hi) > 20 else 1.0
                e_n = float(1 - stats.norm.cdf(cell_thr, mu_n, max(sd_n, 1e-9)))
                e_t = float(stats.norm.cdf(cell_thr, mu_t, max(sd_t, 1e-9)))
                fit["cell_threshold"] = cell_thr
                fit["cell_miscall_normal"] = e_n
                fit["cell_miscall_tumor"] = e_t
                fit["miscall_rate"] = float(max(e_n, e_t))
                log(f"  metacell-mean components {mu_n:+.3f} / {mu_t:+.3f}"
                    f" -> per-cell threshold {cell_thr:+.3f}")
                log(f"  Estimated per-cell miscall rate: normal->tumor {e_n:.3f}"
                    f" / tumor->normal {e_t:.3f}")
        except Exception as exc:
            warn(f"Per-cell karyotype projection failed: {exc}")
            proj = None

    # --- Per-metacell aggregation ---
    idx = S.index
    g = cells.assign(_mc=mc.values, _ng=ng)
    grp = g.groupby("_mc")
    A = pd.DataFrame(index=idx)
    A["n_cells"] = grp.size().reindex(idx).fillna(0).astype(int)
    A["n_callable"] = grp._ng.apply(lambda s: int((s >= 2000).sum())).reindex(idx).fillna(0).astype(int)
    A["median_genes"] = grp._ng.median().reindex(idx)
    for p in panels:
        A["ncell_" + p] = grp.called.apply(lambda s, p=p: int((s == p).sum())).reindex(idx).fillna(0).astype(int)
    A["ncell_doublet"] = grp.called.apply(lambda s: int((s == "doublet").sum())).reindex(idx).fillna(0).astype(int)
    A["n_called"] = A[["ncell_" + p for p in panels]].sum(1)
    for c in S.columns:
        A[("score_" if not c.endswith("_nsig") else "") + c] = S[c].values
    # Number of genes exceeding the calibrated threshold (not a fixed "2x"),
    # kept consistent with the score threshold -- otherwise a metacell could
    # pass the median test but fail on gene count against a mismatched cutoff.
    for p in panels:
        R = RATIOS.get(p)
        A["ngene_over_thr_" + p] = (
            (R.values > thr[p]["threshold"]).sum(1) if R is not None and R.shape[1]
            else 0)

    # --- Measuring mixture (binomial test) ---
    A["cnv_frac_tumor_cells"] = np.nan
    A["cnv_mixture_q"] = np.nan
    A["karyotype_proj"] = (mc_proj.reindex(idx).values if mc_proj is not None
                           else np.nan)
    if proj is not None and fit["usable"]:
        is_t = (proj.values > fit["cell_threshold"])
        kt = pd.Series(is_t, index=adata.obs_names).groupby(mc.values).sum().reindex(idx).fillna(0).values
        nn = A.n_cells.values.astype(float)
        A["cnv_frac_tumor_cells"] = kt / np.maximum(nn, 1)
        eps = max(fit["miscall_rate"], MIN_MISCALL_RATE)
        p1, p2 = _mixture_pvalues(kt, nn, eps)
        # Both nulls must be rejected, so take the larger p-value (intersection test)
        A["cnv_mixture_q"] = _bh(np.maximum(p1, p2))
        log(f"Mixture calling: using error rate {eps:.3f} in the null, FDR {mixture_fdr} ->"
            f" {int((A.cnv_mixture_q < mixture_fdr).sum())} metacells called mixed")

    # --- Labeling ---
    def decide(r):
        mixq = r.cnv_mixture_q
        f = r.cnv_frac_tumor_cells
        if np.isfinite(mixq) and mixq < mixture_fdr:
            return ("Mixed:tumor+normal", "measured-mixture",
                    "per-cell karyotype call: tumor %.0f%% / normal %.0f%% (binomial test q=%.2g)"
                    % (f * 100, (1 - f) * 100, mixq))
        # `side` is decided by comparing the metacell-mean projection to the
        # midpoint between component means. The cell fraction f is
        # sensitive to exactly where the threshold sits, so it is reported
        # but not used to decide.
        mp = r.karyotype_proj
        if np.isfinite(mp) and fit.get("usable"):
            side = "tumor" if mp > fit["cell_threshold"] else "normal"
            side_txt = "karyotype projection %+.2f (%s side)" % (mp, "tumor" if side == "tumor" else "normal")
        else:
            side = "normal"
            side_txt = "no karyotype information"
        # Conditions to accept a type: (a) the per-gene ratio median is at
        # or above the calibrated threshold; (b) at least 3 genes
        # individually exceed their threshold. (a) requires that a majority
        # of the panel is elevated; (b) rules out accepting a type on just
        # 1-2 genes.
        passing = [p for p in panels
                   if np.isfinite(r["score_" + p])
                   and r["score_" + p] >= thr[p]["threshold"]
                   and r["ngene_over_thr_" + p] >= 3]
        if len(passing) == 1:
            t = passing[0]
            lbl = ("Tumor:" + t) if side == "tumor" else t
            return (lbl, "typed",
                    "%s. aggregate score %s median %+.2f is at/above threshold %+.2f"
                    " (%d genes over threshold)"
                    % (side_txt, t, r["score_" + t], thr[t]["threshold"],
                       int(r["ngene_over_thr_" + t])))
        if len(passing) >= 2:
            lbl = ("Tumor:" if side == "tumor" else "") + "Unclassified:multiple-types"
            return (lbl, "unclassified",
                    "%d types simultaneously exceed threshold on aggregate score (%s); cannot narrow to one"
                    % (len(passing), "+".join(sorted(passing))))
        if side == "tumor" and np.isfinite(r.karyotype_proj):
            return ("Tumor:Unclassified", "unclassified",
                    "%s, but no type exceeds threshold; phenotype undetermined" % side_txt)
        if r.n_callable == 0 or r.n_called < min_called_cells:
            return ("Unclassified:low-depth", "unclassified",
                    "%d callable cells (%d with >=2000 detected genes); insufficient depth"
                    % (int(r.n_called), int(r.n_callable)))
        return ("Unclassified:no-dominant-type", "unclassified",
                "%d callable cells present, but no type exceeds threshold" % int(r.n_called))

    res = A.apply(decide, axis=1, result_type="expand")
    res.columns = ["label", "label_kind", "evidence"]
    A = pd.concat([res, A], axis=1)
    log("\nAnnotation results:")
    for k, v in A.label_kind.value_counts().items():
        log(f"  [{k}] {v} metacells")
        for l, c in A[A.label_kind == k].label.value_counts().items():
            log(f"      {c:>4}  {l}")

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        A.to_csv(out_dir / "metacell_annotation.csv")
        summary = {"version": __version__,
                   "panels_used": list(panels), "panels_reported": list(reported),
                   "panel_availability": avail.to_dict(orient="index"),
                   "panel_thresholds": thr,
                   "karyotype": (karyotype.to_dict() if karyotype is not None else None),
                   "projection_fit": fit,
                   "label_counts": A.label.value_counts().to_dict(),
                   "kind_counts": A.label_kind.value_counts().to_dict(),
                   "params": {"panel_fdr": panel_fdr, "mixture_fdr": mixture_fdr,
                              "min_called_cells": min_called_cells, "n_perm": n_perm}}
        (out_dir / "metacell_annotation_report.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
        log(f"Wrote: {out_dir/'metacell_annotation.csv'} and report.json")
    A.attrs["panel_thresholds"] = thr
    A.attrs["projection_fit"] = fit
    return A
