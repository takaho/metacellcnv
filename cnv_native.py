#!/usr/bin/env python3
"""cnv_native.py — Native CNV estimation, without an infercnvpy dependency.

infercnvpy is explicitly marked "experimental" by its own maintainers, and
in practice this caused concrete problems for this project: reference-group
values get rounded to exactly 0 (hiding tumor contamination in the
reference), the dynamic noise threshold depends on chunk size (breaking
reproducibility), and the windowed running-mean signal is easily confounded
by co-expressed gene clusters (operons, HLA, keratins, immunoglobulins) that
mimic a CNV-shaped signal without being one.

This module provides two estimators.

`infercnv_scores()` reproduces infercnvpy's formula (subtract reference mean
-> bounded -> clip -> pyramid running mean -> per-cell median centering ->
dynamic threshold), but computes the dynamic threshold's standard deviation
once over all cells rather than per 5,000-cell chunk, so results no longer
depend on chunk size or cell ordering. Kept for comparison with prior results.

Significance in both estimators is assessed from an empirical per-bin null
(genes randomly reassigned to bins, preserving bin size) rather than a
parametric normal approximation. The null distribution's tails are
substantially heavier than normal, driven by cell-to-cell variation in
expression programs rather than by sampling/count noise, so no distributional
assumption removes it. See README.md for a summary and validation results
and calibration tables.

`bin_composition_cnv()` (recommended) tests raw-count genomic bin composition
directly. It needs no reference group, so it does not suffer from the
reference-rounds-to-zero problem, and it includes three validated controls:
  (a) depth independence: uses each bin's share of a cell's total counts, so
      per-cell depth affects only the variance of that share (modeled with a
      beta-binomial).
  (b) an expression-program control: recomputes the same statistic after
      randomly reassigning genes to bins (bin size preserved), isolating
      amplitude coming from co-expression rather than genuine CNV.
  (c) a genomic-continuity test: checks whether permuting gene order within
      a chromosome changes the amplitude, to test whether sub-chromosomal
      structure is supported (vs. whole-chromosome aneuploidy only).

Both estimators can proceed through `cnv_embedding()` for PCA -> neighbors ->
Leiden clustering (a thin scanpy wrapper, since infercnvpy's own `cnv.tl.pca`
etc. are themselves thin scanpy wrappers).
"""

from __future__ import annotations

import re
from typing import Literal, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

try:
    from scrna_common import log, warn
except Exception:  # pragma: no cover
    def log(msg: str) -> None:
        print(f"[cnv] {msg}", flush=True)

    def warn(msg: str) -> None:
        print(f"[cnv][WARNING] {msg}", flush=True)


__version__ = "1.0"

DEFAULT_WINDOW = 100
DEFAULT_STEP = 10
DEFAULT_LFC_CLIP = 3.0
DEFAULT_DYNAMIC_THRESHOLD = 1.5
DEFAULT_EXCLUDE = ("chrX", "chrY", "chrM")
# Defaults for the bin-composition method
DEFAULT_BIN_GENES = 100
DEFAULT_MIN_BIN_COUNTS = 20
BIN_FDR = 0.05
# Number of iterations for the expression-program control and continuity test
DEFAULT_N_CONTROL = 20
# Default bin-splitting strategy and distribution model, both chosen by
# calibration on real data (see README)
DEFAULT_BIN_BY = "counts"
DEFAULT_DIST = "betabinom"
# Significance is assessed from the empirical null-shuffle distribution
# (a normal approximation misses the tails by 20-160x; see the module
# docstring for the derivation).
DEFAULT_P_METHOD = "empirical"
DEFAULT_N_NULL_SHUFFLE = 100


# ---------------------------------------------------------------------------
# Shared: genomic gene ordering
# ---------------------------------------------------------------------------

def _natural_key(name: str):
    return [int(t) if t.isdigit() else t.lower()
            for t in re.split(r"(\d+)", str(name))]


def ordered_genes(var: pd.DataFrame, chrom_key: str = "chromosome",
                  start_key: str = "start",
                  exclude: Sequence[str] = DEFAULT_EXCLUDE) -> dict[str, np.ndarray]:
    """Return, per chromosome, gene row indices ordered by position."""
    for k in (chrom_key, start_key):
        if k not in var:
            raise ValueError(f"var['{k}'] not found. Annotate gene coordinates first")
    ch = var[chrom_key].astype(str)
    ok = ch.notna() & (ch != "nan") & (ch != "NA") & ~ch.isin([str(e) for e in exclude])
    if not ok.any():
        raise ValueError("No genes with coordinates remain after exclusion")
    out: dict[str, np.ndarray] = {}
    pos = np.arange(len(var))
    for c in sorted(ch[ok].unique(), key=_natural_key):
        m = ok.values & (ch.values == c)
        idx = pos[m]
        order = np.argsort(var[start_key].values[m].astype(float), kind="stable")
        out[c] = idx[order]
    return out


# ---------------------------------------------------------------------------
# Method 1: infercnv-equivalent (windowed running mean)
# ---------------------------------------------------------------------------

def running_mean_pyramid(X: np.ndarray, window: int = DEFAULT_WINDOW,
                         step: int = DEFAULT_STEP) -> np.ndarray:
    """Pyramid-weighted running mean, same definition as infercnvpy's `_running_mean`.

    Weights are 1,2,...,n/2,...,2,1. Uses 'valid' convolution, so the output
    width is n_genes - window + 1, then subsampled every `step`. If there are
    fewer genes than `window`, collapses to one uniform-weight value (same as
    the reference implementation).
    """
    n_genes = X.shape[1]
    if n_genes == 0:
        return np.zeros((X.shape[0], 0))
    if window < n_genes:
        r = np.arange(1, window + 1)
        pyramid = np.minimum(r, r[::-1]).astype(float)
        sm = np.apply_along_axis(
            lambda row: np.convolve(row, pyramid, mode="valid"), 1, X) / pyramid.sum()
        return sm[:, np.arange(0, sm.shape[1], step)]
    w = np.ones(n_genes, dtype=float)
    sm = np.apply_along_axis(lambda row: np.convolve(row, w, mode="valid"), 1, X) / w.sum()
    return sm


def _bounded_center(X: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Difference from reference. With 2+ reference groups, uses the bounded
    method (values within the reference range are treated as 0)."""
    if reference.ndim == 1:
        reference = reference[None, :]
    if reference.shape[0] == 1:
        return X - reference[0, :]
    lo, hi = reference.min(0), reference.max(0)
    out = np.zeros(X.shape, dtype=float)
    above, below = X > hi, X < lo
    out[above] = (X - hi)[above]
    out[below] = (X - lo)[below]
    return out


def infercnv_scores(adata, *, reference_key: str | None = None,
                    reference_cat: Sequence[str] | None = None,
                    reference: np.ndarray | None = None,
                    lfc_clip: float = DEFAULT_LFC_CLIP,
                    window_size: int = DEFAULT_WINDOW,
                    step: int = DEFAULT_STEP,
                    dynamic_threshold: float | None = DEFAULT_DYNAMIC_THRESHOLD,
                    exclude_chromosomes: Sequence[str] = DEFAULT_EXCLUDE,
                    chrom_key: str = "chromosome", start_key: str = "start",
                    layer: str | None = None, key_added: str = "cnv",
                    inplace: bool = True, block: int = 2000):
    """CNV scores equivalent to infercnvpy's `cnv.tl.infercnv`, but chunk-independent.

    Assumes X is **log-normalized** (same as infercnv).

    Differs from the reference implementation in exactly one respect: the
    dynamic threshold's standard deviation is computed once over all cells,
    rather than per 5,000-cell chunk — so which elements get zeroed out no
    longer depends on cell ordering or chunk size.
    """
    X = adata.layers[layer] if layer else adata.X
    if sp.issparse(X):
        X = X.tocsr()

    ref = _resolve_reference(adata, reference_key, reference_cat, reference, layer)
    order = ordered_genes(adata.var, chrom_key, start_key, exclude_chromosomes)

    chr_pos, parts = {}, []
    offset = 0
    for c, idx in order.items():
        chr_pos[c] = offset
        blocks = []
        for s in range(0, adata.n_obs, block):
            e = min(s + block, adata.n_obs)
            xb = X[s:e][:, idx]
            xb = np.asarray(xb.todense()) if sp.issparse(xb) else np.asarray(xb, dtype=float)
            xb = _bounded_center(xb, ref[:, idx])
            np.clip(xb, -lfc_clip, lfc_clip, out=xb)
            blocks.append(running_mean_pyramid(xb, window_size, step))
        part = np.vstack(blocks)
        parts.append(part)
        offset += part.shape[1]
    res = np.hstack(parts)
    log(f"infercnv-equivalent: {res.shape[0]:,} cells x {res.shape[1]:,} windows"
        f" / {len(order)} chromosomes / window {window_size} genes, step {step}")

    # Center by per-cell median
    res -= np.median(res, axis=1, keepdims=True)

    if dynamic_threshold is not None:
        thr = float(dynamic_threshold * np.std(res))   # computed once, over all cells
        n_before = int((res != 0).sum())
        res[np.abs(res) < thr] = 0.0
        log(f"dynamic threshold {thr:.4f} (SD computed once over all cells): "
            f"nonzero {n_before:,} -> {int((res != 0).sum()):,}")

    out = sp.csr_matrix(res)
    if inplace:
        adata.obsm[f"X_{key_added}"] = out
        adata.uns[key_added] = {"chr_pos": chr_pos}
        return None
    return chr_pos, out


def _resolve_reference(adata, reference_key, reference_cat, reference, layer):
    """Reference profile (category x gene). Equivalent to infercnvpy's `_get_reference`."""
    X = adata.layers[layer] if layer else adata.X
    if reference is not None:
        ref = np.asarray(reference, dtype=float)
        return ref[None, :] if ref.ndim == 1 else ref
    if reference_key is None or reference_cat is None:
        warn("No reference group specified; using the mean over all cells as"
             " reference. CNV will be underestimated if the reference contains tumor cells")
        m = X.mean(axis=0)
        return np.asarray(m).reshape(1, -1)
    obs = adata.obs[reference_key].astype(str)
    cats = [reference_cat] if isinstance(reference_cat, str) else list(reference_cat)
    missing = [c for c in cats if c not in set(obs)]
    if missing:
        raise ValueError(f"Reference categories not found in obs['{reference_key}']: {missing}")
    rows = []
    for c in cats:
        m = (obs.values == c)
        if m.sum() == 0:
            raise ValueError(f"No cells found for reference category {c}")
        if m.sum() < 3:
            warn(f"Reference category {c} has only {m.sum()} cells."
                 " The reference mean will be unstable (the bounded method's"
                 " lower/upper bound would be set by a single cell's noise)")
        rows.append(np.asarray(X[m].mean(axis=0)).ravel())
    if len(rows) == 1:
        warn("Only one reference category, so this reduces to a simple"
             " difference. Cell-type-specific gene clusters are more likely"
             " to be mistaken for CNV (the bounded method would help with 2+ categories)")
    return np.vstack(rows)


# ---------------------------------------------------------------------------
# Method 2: genomic bin count composition (recommended)
# ---------------------------------------------------------------------------

def _split_by_weight(idx: np.ndarray, w: np.ndarray, n_bins: int) -> list[np.ndarray]:
    """Split genome-ordered genes into intervals of equal cumulative weight.

    With uniform weight (w=1) this reduces to splitting genes evenly. Using
    each bin's total count as the weight instead equalizes expected detection
    power (E[k]) across bins.
    """
    if n_bins <= 1 or idx.size == 0:
        return [idx]
    cw = np.cumsum(np.maximum(w, 0.0))
    total = cw[-1]
    if total <= 0:
        return [a for a in np.array_split(idx, n_bins) if a.size]
    edges = np.searchsorted(cw, np.linspace(0, total, n_bins + 1)[1:-1], side="left")
    parts, prev = [], 0
    for e in list(edges) + [idx.size]:
        e = int(min(max(e, prev), idx.size))
        if e > prev:
            parts.append(idx[prev:e])
            prev = e
    return parts or [idx]


def genomic_bins(var: pd.DataFrame, bin_genes: int = DEFAULT_BIN_GENES,
                 chrom_key: str = "chromosome", start_key: str = "start",
                 exclude: Sequence[str] = DEFAULT_EXCLUDE,
                 whole_chromosome: bool = False,
                 gene_weight: np.ndarray | None = None,
                 bin_by: Literal["genes", "counts", "detected"] = "genes",
                 ) -> tuple[pd.Series, pd.DataFrame]:
    """Group genes into genomic-order bins.

    `bin_by` controls how bins are split, which directly affects detection
    power uniformity:

    - ``"genes"``: splits by gene count (same idea as infercnv's fixed
      window_size=100-gene window). **Regions of high vs. low gene density
      then get very different total counts per bin, so detection power is
      not uniform across bins.**
    - ``"counts"``: splits by total counts per bin, equalizing E[k_ib] and
      therefore the binomial/beta-binomial variance ≈ n p(1-p), giving more
      uniform detection power.
    - ``"detected"``: splits by number of detections (nonzero cell x gene
      entries) per bin, downweighting high-dropout genes.

    `bin_genes` only sets the *target number of bins*
    (n_bins = chromosome gene count / bin_genes), so switching `bin_by`
    keeps roughly the same bin count and comparisons stay fair.

    `whole_chromosome=True` collapses each chromosome to one bin (whole-
    chromosome aneuploidy only).
    """
    order = ordered_genes(var, chrom_key, start_key, exclude)
    if bin_by != "genes" and gene_weight is None:
        raise ValueError(f"bin_by='{bin_by}' requires gene_weight (a per-gene weight)")
    gene_bin = pd.Series(pd.NA, index=var.index, dtype="object")

    # Bin count allocation. Assigning "gene count / bin_genes" bins per
    # chromosome independently would give low-expression chromosomes the
    # same bin count as high-expression ones, so expected counts per bin
    # wouldn't match. When bin_by is weighted, allocate the total bin count
    # by each chromosome's share of the weight instead.
    total_bins = max(1, int(round(sum(len(i) for i in order.values()) / bin_genes)))
    n_bins_per: dict[str, int] = {}
    if bin_by == "genes" or whole_chromosome:
        for c, idx in order.items():
            n_bins_per[c] = max(1, int(round(len(idx) / bin_genes)))
    else:
        w = np.asarray(gene_weight, float)
        share = {c: max(w[idx].sum(), 0.0) for c, idx in order.items()}
        tot = sum(share.values()) or 1.0
        for c, idx in order.items():
            n_bins_per[c] = int(min(len(idx), max(1, round(total_bins * share[c] / tot))))

    rows = []
    for c, idx in order.items():
        if whole_chromosome:
            chunks = [idx]
        else:
            n_bins = n_bins_per[c]
            if bin_by == "genes":
                chunks = [a for a in np.array_split(idx, n_bins) if a.size]
            else:
                chunks = _split_by_weight(idx, np.asarray(gene_weight, float)[idx], n_bins)
        for b, ch in enumerate(chunks):
            if len(ch) == 0:
                continue
            bid = f"{c}:{b:03d}"
            gene_bin.iloc[ch] = bid
            row = {"bin": bid, "chromosome": c, "n_genes": len(ch),
                   "start": float(var[start_key].values[ch].astype(float).min()),
                   "end": float(var[start_key].values[ch].astype(float).max())}
            if gene_weight is not None:
                row["weight"] = float(np.asarray(gene_weight, float)[ch].sum())
            rows.append(row)
    info = pd.DataFrame(rows).set_index("bin")
    msg = (f"Genomic bins ({bin_by}): {len(info)} bins"
           f" / {info.chromosome.nunique()} chromosomes"
           f" / median genes/bin {info.n_genes.median():.0f}")
    if "weight" in info:
        cv = info["weight"].std() / max(info["weight"].mean(), 1e-12)
        msg += f" / bin weight CV {cv:.3f}"
    log(msg)
    return gene_bin, info


def _bin_counts(X, gene_bin: pd.Series, bins: pd.Index) -> np.ndarray:
    """Cell x bin raw count totals, computed via a single sparse matrix product."""
    codes = pd.Categorical(gene_bin.values, categories=list(bins)).codes
    keep = codes >= 0
    S = sp.csr_matrix((np.ones(keep.sum()),
                       (np.where(keep)[0], codes[keep])),
                      shape=(len(gene_bin), len(bins)))
    W = X @ S
    return np.asarray(W.todense() if sp.issparse(W) else W, dtype=float)


def _nb_alpha(k: np.ndarray, mu: np.ndarray, trim: float = 0.90) -> float:
    """Estimate the negative binomial dispersion parameter alpha by the method of moments.

    Uses Var(k) = mu (1 + alpha*mu), solved as Σ[(k-mu)² - mu] / Σ mu².

    Why offer a negative-binomial option: subtracting a mean expression
    value (infercnv's approach) is unstable for genes that drop out in most
    cells, since a single gene's mean log-normalized value is then dominated
    by the zero fraction. Modeling counts directly avoids that.

    This module's primary unit is a bin total k_ib (a sum over ~100 genes),
    so individual-gene dropout is averaged out within the sum. Both the
    beta-binomial (the marginal of a Dirichlet-multinomial) and the negative
    binomial are overdispersion models for bin totals, and neither is
    sensitive to zero-inflation the way mean-subtraction is. Which is better
    is decided empirically (see `null_calibration()`).
    """
    ok = np.isfinite(k) & np.isfinite(mu) & (mu > 0)
    if ok.sum() < 20:
        return 0.0
    ratio = np.where(ok, k / np.maximum(mu, 1e-9), np.nan)
    cut = np.nanquantile(ratio[ok], trim)
    use = ok & (ratio <= cut)
    if use.sum() < 20:
        use = ok
    num = float(np.sum((k[use] - mu[use]) ** 2 - mu[use]))
    den = float(np.sum(mu[use] ** 2))
    if den <= 0:
        return 0.0
    return float(max(num / den, 0.0))


def _betabinom_phi(k: np.ndarray, n: np.ndarray, p: float) -> float:
    """Estimate per-bin overdispersion phi by the method of moments.

    Uses Var(k/n) = p(1-p)/n * (1 + phi(n-1)); phi=0 is the binomial case.
    Trims the top 10% before fitting, to keep cells that carry a CNV from
    inflating the null they're being tested against (the same fix used for
    dead-cell/MT-fraction calls elsewhere in this pipeline).
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(n > 0, k / np.maximum(n, 1), np.nan)
    ok = np.isfinite(frac) & (n > 0)
    if ok.sum() < 20:
        return 0.0
    cut = np.nanquantile(frac[ok], 0.90)
    use = ok & (frac <= cut)
    if use.sum() < 20:
        use = ok
    obs_var = float(np.nanvar(frac[use]))
    exp_var = float(np.mean(p * (1 - p) / n[use]))
    if exp_var <= 0:
        return 0.0
    nbar = float(np.mean(n[use]))
    phi = (obs_var / exp_var - 1.0) / max(nbar - 1.0, 1.0)
    return float(max(phi, 0.0))


def _betabinom_phi_vec(W: np.ndarray, n: np.ndarray, p_b: np.ndarray,
                       trim: float = 0.90) -> np.ndarray:
    """Estimate per-bin overdispersion phi for **all bins at once** (vectorized `_betabinom_phi`).

    Running `_betabinom_phi` per bin in a Python loop measured 20 ms per
    shuffle (365 metacells x 352 bins), 57% of total CNV compute time.
    Vectorizing the same formula across bins is an order of magnitude
    faster, with identical results.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = W / np.maximum(n, 1)[:, None]
    frac = np.where(n[:, None] > 0, frac, np.nan)
    cut = np.nanquantile(frac, trim, axis=0)
    use = np.isfinite(frac) & (frac <= cut[None, :])
    enough = use.sum(0) >= 20
    use = np.where(enough[None, :], use, np.isfinite(frac))
    F = np.where(use, frac, np.nan)
    obs_var = np.nanvar(F, axis=0)
    inv_n = np.where(use, 1.0 / np.maximum(n, 1)[:, None], np.nan)
    exp_var = (p_b * (1 - p_b)) * np.nanmean(inv_n, axis=0)
    nbar = np.nanmean(np.where(use, n[:, None], np.nan), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        phi = (obs_var / exp_var - 1.0) / np.maximum(nbar - 1.0, 1.0)
    phi = np.where(np.isfinite(phi) & (exp_var > 0), phi, 0.0)
    return np.maximum(phi, 0.0)


def _nb_alpha_vec(W: np.ndarray, mu: np.ndarray, trim: float = 0.90) -> np.ndarray:
    """Estimate negative-binomial alpha for all bins at once (vectorized `_nb_alpha`)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = W / np.maximum(mu, 1e-9)
    ratio = np.where(mu > 0, ratio, np.nan)
    cut = np.nanquantile(ratio, trim, axis=0)
    use = np.isfinite(ratio) & (ratio <= cut[None, :])
    enough = use.sum(0) >= 20
    use = np.where(enough[None, :], use, np.isfinite(ratio))
    num = np.nansum(np.where(use, (W - mu) ** 2 - mu, np.nan), axis=0)
    den = np.nansum(np.where(use, mu ** 2, np.nan), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        a = num / den
    return np.maximum(np.where(np.isfinite(a) & (den > 0), a, 0.0), 0.0)


def _dispersion_z(W: np.ndarray, exp: np.ndarray, n_i: np.ndarray,
                  p_b: np.ndarray, dist: str,
                  fit_mask: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Estimate per-bin overdispersion and return z.

    - ``binom``     : Var = n p (1-p). Ignores overdispersion, so z comes
                      out far too large (measured null SD(z) = 6.3). Do not use.
    - ``betabinom`` : Var = n p (1-p) [1 + phi (n-1)]. The natural choice for
                      compositional data; corresponds to the marginal of a
                      Dirichlet-multinomial.
    - ``nb``        : Var = mu (1 + alpha*mu), the standard scRNA-seq count
                      model, with mu = n_i p_b here. Ignores the compositional
                      constraint (Σ_b k_ib = n_i), so it's conservative when
                      there are few bins.

    `fit_mask` selects **which cells are used to estimate overdispersion**.
    This matters a great deal: estimating overdispersion including cells
    that actually carry a CNV lets those cells inflate their own null and
    become undetectable (measured, with chr1 amplified 1.5x in 20% of
    metacells: median detection power drops to 0.014 when all cells are used
    vs. a normal-only reference — see README for the full table).

    Trimming the top 10% alone (`_betabinom_phi`'s `fit_quantile`) is not
    enough once the CNV-carrying fraction exceeds 10%; the same issue shows
    up in this pipeline's dead-cell/MT-fraction calls.
    """
    nb = len(p_b)
    m = np.ones(W.shape[0], dtype=bool) if fit_mask is None else np.asarray(fit_mask, bool)
    if m.sum() < 20:
        m = np.ones(W.shape[0], dtype=bool)
    if dist == "binom":
        disp = np.zeros(nb)
        var = exp * (1 - p_b)[None, :]
    elif dist == "betabinom":
        disp = _betabinom_phi_vec(W[m], n_i[m], p_b)
        var = exp * (1 - p_b)[None, :] * (1 + disp[None, :] * np.maximum(n_i - 1, 1)[:, None])
    elif dist == "nb":
        disp = _nb_alpha_vec(W[m], exp[m])
        var = exp * (1 + disp[None, :] * exp)
    else:
        raise ValueError("dist must be one of 'betabinom' / 'nb' / 'binom'")
    z = (W - exp) / np.sqrt(np.maximum(var, 1e-9))
    return z, disp


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if not ok.any():
        return q
    x = p[ok]
    o = np.argsort(x)
    m = x.size
    adj = np.minimum.accumulate((x[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    tmp = np.empty(m)
    tmp[o] = np.clip(adj, 0, 1)
    q[ok] = tmp
    return q


def _shuffle_gene_bin(gene_bin: pd.Series, info: pd.DataFrame,
                      rng: np.random.Generator) -> pd.Series:
    """Randomly reassign genes to bins, preserving bin sizes.

    Breaks only genomic adjacency; per-gene expression level, dropout rate,
    and per-cell depth are all preserved. This makes a null that still
    contains any signal coming from co-expression rather than true CNV.
    """
    assigned = np.where(gene_bin.notna().values)[0]
    perm = rng.permutation(assigned)
    out = pd.Series(pd.NA, index=gene_bin.index, dtype="object")
    s = 0
    for bid, sz in zip(info.index, info["n_genes"].values):
        out.iloc[perm[s:s + int(sz)]] = bid
        s += int(sz)
    return out


def bin_composition_cnv(adata, *, bin_genes: int = DEFAULT_BIN_GENES,
                        whole_chromosome: bool = False,
                        reference_mask: np.ndarray | None = None,
                        layer: str | None = "counts",
                        chrom_key: str = "chromosome", start_key: str = "start",
                        exclude_chromosomes: Sequence[str] = DEFAULT_EXCLUDE,
                        min_bin_counts: int = DEFAULT_MIN_BIN_COUNTS,
                        fdr: float = BIN_FDR,
                        n_control: int = DEFAULT_N_CONTROL,
                        n_null_shuffle: int = DEFAULT_N_NULL_SHUFFLE,
                        p_method: Literal["empirical", "normal"] = DEFAULT_P_METHOD,
                        bin_by: Literal["genes", "counts", "detected"] = DEFAULT_BIN_BY,
                        dist: Literal["betabinom", "nb", "binom"] = DEFAULT_DIST,
                        seed: int = 0, key_added: str = "cnv",
                        inplace: bool = True) -> dict:
    """Derive CNV log2 ratios from genomic bin count composition, with FDR.

    For cell i, bin b:
        k_ib = bin b's total raw count,  n_i = total count over coordinate-annotated genes
        p_b  = bin b's count share in the reference (default: all cells)
        log2 ratio = log2( (k_ib + 1) / (n_i p_b + 1) )
    Tested per bin with a beta-binomial (overdispersion phi_b estimated by
    the method of moments).

    A key difference from infercnv: the reference can be "all cells".
    infercnv rounds reference-group values to 0, hiding tumor contamination
    in the reference; here, using all cells as reference means the statistic
    is a deviation from the overall mean, so if tumor cells dominate, normal
    cells necessarily show up on the opposite side.
    """
    X = adata.layers[layer] if (layer and layer in adata.layers) else adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsr()

    gw = None
    if bin_by == "counts":
        gw = np.asarray(X.sum(0)).ravel().astype(float)
    elif bin_by == "detected":
        gw = np.asarray((X > 0).sum(0)).ravel().astype(float)
    gene_bin, info = genomic_bins(adata.var, bin_genes, chrom_key, start_key,
                                  exclude_chromosomes, whole_chromosome,
                                  gene_weight=gw, bin_by=bin_by)
    bins = info.index
    W = _bin_counts(X, gene_bin, bins)                  # cells x bins
    n_i = W.sum(1)
    if np.median(n_i) < min_bin_counts * len(bins) / 50:
        warn(f"Median coordinate-annotated count is only {np.median(n_i):.0f}."
             " Larger bins (increase bin_genes) would be more stable")

    if reference_mask is None:
        ref = np.ones(adata.n_obs, dtype=bool)
        warn("No reference cells specified. Bin expected fraction p_b and"
             " overdispersion phi_b will be estimated from all cells, so"
             " CNV-carrying cells inflate their own null and detection power"
             " drops sharply (measured: power 0.014 at 20% spike-in)."
             " Pass reference_mask if any cells are known to be normal")
    else:
        ref = np.asarray(reference_mask, dtype=bool)
    if ref.sum() < 10:
        raise ValueError(f"Only {ref.sum()} reference cells available")
    p_b = W[ref].sum(0) / max(W[ref].sum(), 1.0)

    exp = np.outer(n_i, p_b)
    lr = np.log2((W + 1.0) / (exp + 1.0))

    # --- Per-bin overdispersion model ---
    z, disp = _dispersion_z(W, exp, n_i, p_b, dist, fit_mask=ref)

    # Parametric p-values are not used: the null's SD(z) (gene-to-bin
    # shuffle) is not 1 for any of the distribution choices tested — the
    # excess variance comes from cell-to-cell expression-program
    # differences, not sampling noise, so no distributional assumption
    # removes it. z is instead calibrated by dividing by the null's SD
    # (see README.md for a summary).
    z_scale, resolution = 1.0, None
    if n_null_shuffle > 0:
        rng0 = np.random.default_rng(seed + 7919)
        null = np.empty((n_null_shuffle, adata.n_obs, len(bins)), dtype=np.float32)
        for i in range(n_null_shuffle):
            gb = _shuffle_gene_bin(gene_bin, info, rng0)
            Wn = _bin_counts(X, gb, bins)
            nn = Wn.sum(1)
            pn = Wn[ref].sum(0) / max(Wn[ref].sum(), 1.0)
            zn, _ = _dispersion_z(Wn, np.outer(nn, pn), nn, pn, dist, fit_mask=ref)
            null[i] = zn
        N = np.abs(null.reshape(-1, len(bins)).astype(float))
        z_scale = float(max(np.nanstd(N), 1e-9))
        log(f"null (genes randomly reassigned to bins), {n_null_shuffle} shuffles:"
            f" SD(|z|) {z_scale:.3f} / {N.shape[0]:,} draws/bin")

    if p_method == "empirical" and n_null_shuffle > 0:
        # p is built per bin from that bin's own null ECDF. See the module
        # docstring for why a normal approximation is not used.
        pv = np.empty_like(z, dtype=float)
        for j in range(len(bins)):
            nul = np.sort(N[:, j])
            r = np.searchsorted(nul, np.abs(z[:, j]), side="left")
            pv[:, j] = ((nul.size - r) + 1.0) / (nul.size + 1.0)
        resolution = 1.0 / (N.shape[0] + 1.0)
        log(f"Empirical p resolution floor {resolution:.2e}"
            f" (cannot produce smaller p; increase n_null_shuffle if needed)")
    else:
        if p_method == "empirical":
            warn("n_null_shuffle=0, so empirical p cannot be built. Falling back to normal approximation")
        z = z / z_scale
        pv = 2 * stats.norm.sf(np.abs(z))

    q = _bh(pv.ravel()).reshape(pv.shape)
    n_sig = int((q < fdr).sum())
    if resolution is not None and n_sig:
        # BH compares the k-th smallest p to k*alpha/m. If the smallest
        # required p is below the empirical resolution, that region can't
        # be resolved.
        need = fdr * max(n_sig, 1) / q.size
        if need < resolution:
            warn(f"BH requires p ({need:.2e}) below the empirical resolution"
                 f" ({resolution:.2e}). Increase n_null_shuffle to at least"
                 f" {int(np.ceil(1.0 / need / adata.n_obs)) + 1}")
    info["dispersion"] = disp
    info["frac_cells_significant"] = (q < fdr).mean(0)
    info["z_sd"] = z.std(0)

    res: dict = {"log2_ratio": lr, "z": z, "q": q, "bins": info,
                 "p_bin": p_b, "n_counts": n_i, "dist": dist, "bin_by": bin_by,
                 "gene_bin": gene_bin, "z_scale": z_scale,
                 "p_value": pv, "p_method": p_method, "resolution": resolution,
                 "amplitude": float(np.mean(np.abs(lr)))}

    # --- (b) Expression-program control: genes randomized, bin size preserved ---
    if n_control > 0:
        rng = np.random.default_rng(seed)
        assigned = gene_bin.notna().values
        sizes = info["n_genes"].values
        amps = []
        for _ in range(n_control):
            perm = rng.permutation(np.where(assigned)[0])
            shuffled = pd.Series(pd.NA, index=gene_bin.index, dtype="object")
            s = 0
            for bid, sz in zip(bins, sizes):
                shuffled.iloc[perm[s:s + sz]] = bid
                s += sz
            Wc = _bin_counts(X, shuffled, bins)
            pc = Wc[ref].sum(0) / max(Wc[ref].sum(), 1.0)
            lrc = np.log2((Wc + 1.0) / (np.outer(Wc.sum(1), pc) + 1.0))
            amps.append(float(np.mean(np.abs(lrc))))
        res["control_amplitude"] = float(np.mean(amps))
        res["control_amplitude_sd"] = float(np.std(amps))
        res["amplitude_ratio"] = res["amplitude"] / max(res["control_amplitude"], 1e-12)
        log(f"amplitude {res['amplitude']:.4f} / expression-program control"
            f" {res['control_amplitude']:.4f} +/- {res['control_amplitude_sd']:.4f}"
            f" -> ratio {res['amplitude_ratio']:.2f}x")
        # Rule of thumb: with a genuine aneuploidy (Case1, whole-chromosome
        # bins, normal reference) the ratio measured 1.43x. The control
        # (randomly reassigned genes) itself measures 0.131, driven by
        # cell-to-cell differences in expression programs and not removed
        # by more sequencing depth — so roughly 70% of observed amplitude is
        # not CNV in general. Spike-in tests still detected a 1.5x change in
        # 69% of cells, so a low ratio doesn't mean the method is unusable;
        # 1.2x is used as a rough lower bound.
        if res["amplitude_ratio"] < 1.2:
            warn(f"Ratio to control is only {res['amplitude_ratio']:.2f}x."
                 " The observed 'CNV' signal is at a level explainable by"
                 " co-expression or cell-type-specific gene clusters. Try"
                 " larger bins or whole-chromosome resolution (whole_chromosome=True)")

    if inplace:
        adata.obsm[f"X_{key_added}"] = sp.csr_matrix(lr)
        chr_pos, off = {}, 0
        for c, sub in info.groupby("chromosome", sort=False):
            chr_pos[c] = off
            off += len(sub)
        adata.uns[key_added] = {"chr_pos": chr_pos, "method": "bin_composition",
                                "amplitude": res["amplitude"],
                                "control_amplitude": res.get("control_amplitude"),
                                "amplitude_ratio": res.get("amplitude_ratio")}
    return res


def null_calibration(adata, *, bin_genes: int = DEFAULT_BIN_GENES,
                     layer: str | None = "counts",
                     bin_by_options: Sequence[str] = ("genes", "counts", "detected"),
                     dist_options: Sequence[str] = ("binom", "betabinom", "nb"),
                     n_shuffle: int = 5,
                     chrom_key: str = "chromosome", start_key: str = "start",
                     exclude_chromosomes: Sequence[str] = DEFAULT_EXCLUDE,
                     seed: int = 0) -> pd.DataFrame:
    """Choose a bin-splitting strategy and distribution model by calibration.

    Computes z on a null genome (genes randomly reassigned to bins) — a
    correctly calibrated model should have SD(z) ~= 1 and P(|z| > 1.96) ~=
    0.05; a model that exceeds 1 overstates significance and shouldn't be used.

    Also measures **detection power uniformity**: the CV of per-bin null-z
    SD, and the CV of per-bin expected counts (mean(n) * p_b). Fixed-gene-
    count bins have widely varying expected counts due to gene density and
    expression differences, giving non-uniform power; count-equalized bins
    are meant to reduce this.
    """
    X = adata.layers[layer] if (layer and layer in adata.layers) else adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsr()
    rng = np.random.default_rng(seed)
    rows = []
    for bb in bin_by_options:
        gw = None
        if bb == "counts":
            gw = np.asarray(X.sum(0)).ravel().astype(float)
        elif bb == "detected":
            gw = np.asarray((X > 0).sum(0)).ravel().astype(float)
        gene_bin, info = genomic_bins(adata.var, bin_genes, chrom_key, start_key,
                                      exclude_chromosomes, False,
                                      gene_weight=gw, bin_by=bb)
        bins = info.index
        # Spread of expected counts across the actual bins (= spread of power)
        W_real = _bin_counts(X, gene_bin, bins)
        n_i = W_real.sum(1)
        p_real = W_real.sum(0) / max(W_real.sum(), 1.0)
        exp_counts = float(np.mean(n_i)) * p_real
        cv_power = float(np.std(exp_counts) / max(np.mean(exp_counts), 1e-12))

        for dd in dist_options:
            sds, exceed, cv_z = [], [], []
            for _ in range(n_shuffle):
                gb = _shuffle_gene_bin(gene_bin, info, rng)
                W = _bin_counts(X, gb, bins)
                n = W.sum(1)
                p = W.sum(0) / max(W.sum(), 1.0)
                z, _ = _dispersion_z(W, np.outer(n, p), n, p, dd)
                z = z[np.isfinite(z)]
                sds.append(float(np.std(z)))
                exceed.append(float(np.mean(np.abs(z) > 1.96)))
                zz, _ = _dispersion_z(W, np.outer(n, p), n, p, dd)
                per_bin = np.nanstd(zz, axis=0)
                cv_z.append(float(np.nanstd(per_bin) / max(np.nanmean(per_bin), 1e-12)))
            rows.append({
                "bin_by": bb, "dist": dd, "n_bins": len(bins),
                "null_z_sd": float(np.mean(sds)),
                "null_frac_gt196": float(np.mean(exceed)),
                "cv_expected_counts": cv_power,
                "cv_per_bin_null_z_sd": float(np.mean(cv_z)),
            })
    out = pd.DataFrame(rows)
    log("calibration (null = genes randomly reassigned to bins, bin sizes preserved):")
    for _, r in out.iterrows():
        log(f"  {r.bin_by:<9} {r.dist:<10} bins {int(r.n_bins):4d}"
            f" / null SD(z) {r.null_z_sd:5.2f}"
            f" / P(|z|>1.96) {r.null_frac_gt196:6.3f} (target 0.050)"
            f" / expected-count CV {r.cv_expected_counts:5.2f}"
            f" / per-bin SD(z) CV {r.cv_per_bin_null_z_sd:5.3f}")
    return out


def power_by_bin(adata, *, effect: float = 1.5, frac_cells: float = 0.20,
                 target_chromosome: str | None = None,
                 bin_genes: int = DEFAULT_BIN_GENES,
                 bin_by_options: Sequence[str] = ("genes", "counts", "detected"),
                 dist: str = DEFAULT_DIST, layer: str | None = "counts",
                 chrom_key: str = "chromosome", start_key: str = "start",
                 exclude_chromosomes: Sequence[str] = DEFAULT_EXCLUDE,
                 fdr: float = BIN_FDR, seed: int = 0) -> pd.DataFrame:
    """Spike in an artificial CNV and measure whether detection power is uniform across bins.

    Multiplies one chromosome's gene counts by `effect` in a subset of
    cells, then measures the fraction of those cells detected at the given
    FDR, per bin. Used to empirically check whether a fixed gene-count
    window gives non-uniform power across regions of differing gene density
    and expression level.
    """
    import anndata as ad_mod

    X = adata.layers[layer] if (layer and layer in adata.layers) else adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsr().astype(float)
    order = ordered_genes(adata.var, chrom_key, start_key, exclude_chromosomes)
    chrom = target_chromosome or max(order, key=lambda c: len(order[c]))
    cols = order[chrom]
    rng = np.random.default_rng(seed)
    hit = rng.random(adata.n_obs) < frac_cells

    # Multiply only the selected-cells x selected-genes submatrix by effect.
    # diags(row mask) @ X @ diags(col mask) is exactly that submatrix.
    colind = np.zeros(X.shape[1]); colind[cols] = 1.0
    Xp = X + (effect - 1.0) * (sp.diags(hit.astype(float)) @ X @ sp.diags(colind))
    Xp.data = np.rint(Xp.data)
    Xp.eliminate_zeros()

    rows = []
    for bb in bin_by_options:
      for ref_tag, rmask in (("no_reference", None), ("normal_reference", ~hit)):
        tmp = ad_mod.AnnData(Xp.copy(), var=adata.var.copy())
        tmp.obs_names = adata.obs_names
        tmp.layers["counts"] = Xp.copy()
        r = bin_composition_cnv(tmp, bin_genes=bin_genes, bin_by=bb, dist=dist,
                                chrom_key=chrom_key, start_key=start_key,
                                exclude_chromosomes=exclude_chromosomes,
                                reference_mask=rmask,
                                n_control=0, n_null_shuffle=3, fdr=fdr,
                                inplace=False, seed=seed)
        info = r["bins"]
        target = info.index[info["chromosome"] == chrom]
        jj = [info.index.get_loc(b) for b in target]
        det = (r["q"][np.ix_(hit, jj)] < fdr).mean(0)
        # Control: false positive rate on chromosomes that were not spiked
        other = [info.index.get_loc(b) for b in info.index if info.loc[b, "chromosome"] != chrom]
        fpr = float((r["q"][np.ix_(hit, other)] < fdr).mean())
        rows.append({"bin_by": bb, "reference": ref_tag,
                     "chromosome": chrom, "n_target_bins": len(jj),
                     "power_median": float(np.median(det)),
                     "power_min": float(det.min()), "power_max": float(det.max()),
                     "power_cv": float(np.std(det) / max(np.mean(det), 1e-12)),
                     "false_positive_rate_other_chrom": fpr,
                     "n_genes_cv": float(info.loc[target, "n_genes"].std()
                                         / max(info.loc[target, "n_genes"].mean(), 1e-12))})
    out = pd.DataFrame(rows)
    log(f"Spike-in power ({chrom} at {effect}x, {int(hit.sum()):,} cells):")
    for _, x in out.iterrows():
        log(f"  {x.bin_by:<9} {x.reference:<16} target bins {int(x.n_target_bins):3d}"
            f" / power median {x.power_median:.3f}"
            f" [{x.power_min:.3f}, {x.power_max:.3f}] / CV {x.power_cv:.3f}"
            f" / other-chrom FPR {x.false_positive_rate_other_chrom:.4f}")
    return out


def contiguity_test(adata, *, bin_genes: int = DEFAULT_BIN_GENES,
                    n_perm: int = DEFAULT_N_CONTROL,
                    layer: str | None = "counts",
                    chrom_key: str = "chromosome", start_key: str = "start",
                    exclude_chromosomes: Sequence[str] = DEFAULT_EXCLUDE,
                    seed: int = 0) -> dict:
    """Test whether sub-chromosomal structure is supported.

    Permutes gene order **within each chromosome**, preserving the
    chromosome mean. If there is segmental CNV, between-bin variance should
    shrink under this permutation; if it doesn't shrink, only whole-
    chromosome aneuploidy is supported (this was the case for Case1).
    """
    X = adata.layers[layer] if (layer and layer in adata.layers) else adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    gene_bin, info = genomic_bins(adata.var, bin_genes, chrom_key, start_key,
                                  exclude_chromosomes, False)
    order = ordered_genes(adata.var, chrom_key, start_key, exclude_chromosomes)
    bins = info.index

    def within_chrom_spread(gb: pd.Series) -> float:
        W = _bin_counts(X, gb, bins)
        n = W.sum(1)
        p = W.sum(0) / max(W.sum(), 1.0)
        lr = np.log2((W + 1.0) / (np.outer(n, p) + 1.0))
        d = pd.DataFrame(lr.T, index=bins)
        d["chromosome"] = info["chromosome"].values
        # SD of residuals after subtracting the within-chromosome bin mean =
        # magnitude of sub-chromosomal structure
        resid = d.groupby("chromosome", sort=False).transform(lambda s: s - s.mean())
        return float(np.nanstd(resid.values))

    obs = within_chrom_spread(gene_bin)
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_perm):
        gb = gene_bin.copy()
        for c, idx in order.items():
            perm = rng.permutation(idx)
            gb.iloc[idx] = gene_bin.iloc[perm].values
        null.append(within_chrom_spread(gb))
    null = np.asarray(null)
    p = float((np.sum(null >= obs) + 1) / (len(null) + 1))
    out = {"observed": obs, "null_mean": float(null.mean()),
           "null_sd": float(null.std()), "p_value": p,
           "ratio": obs / max(float(null.mean()), 1e-12),
           "supports_segmental": bool(p < 0.05 and obs > null.mean())}
    log(f"Genomic continuity: observed {obs:.4f} / null {null.mean():.4f} +/- {null.std():.4f}"
        f" -> ratio {out['ratio']:.2f}, p={p:.4f}"
        f" -> sub-chromosomal structure {'supported' if out['supports_segmental'] else 'not supported'}")
    return out


# ---------------------------------------------------------------------------
# Embedding and clustering (replaces infercnvpy's thin scanpy wrappers)
# ---------------------------------------------------------------------------

def cnv_embedding(adata, key: str = "cnv", n_comps: int = 50,
                  n_neighbors: int = 15, resolution: float = 1.0,
                  cluster_key: str | None = None, seed: int = 0) -> None:
    """Run PCA -> neighbors -> Leiden on X_cnv.

    infercnvpy's `cnv.tl.pca` / `cnv.pp.neighbors` / `cnv.tl.leiden` are just
    wrappers that run scanpy on a temporary AnnData exposing obsm['X_cnv'] as
    X; this does the same thing directly.
    """
    import anndata as ad
    import scanpy as sc

    if f"X_{key}" not in adata.obsm:
        raise KeyError(f"obsm['X_{key}'] not found. Run CNV estimation first")
    Xc = adata.obsm[f"X_{key}"]
    tmp = ad.AnnData(sp.csr_matrix(Xc) if not sp.issparse(Xc) else Xc.copy())
    tmp.obs_names = adata.obs_names
    n_comps = int(max(2, min(n_comps, min(tmp.shape) - 1)))
    sc.pp.pca(tmp, n_comps=n_comps, zero_center=False, random_state=seed)
    adata.obsm[f"X_{key}_pca"] = tmp.obsm["X_pca"]
    sc.pp.neighbors(tmp, use_rep="X_pca", n_neighbors=int(min(n_neighbors, tmp.n_obs - 1)))
    ck = cluster_key or f"{key}_leiden"
    try:
        sc.tl.leiden(tmp, resolution=resolution, key_added=ck,
                     flavor="igraph", n_iterations=2, directed=False,
                     random_state=seed)
    except TypeError:
        sc.tl.leiden(tmp, resolution=resolution, key_added=ck, random_state=seed)
    adata.obs[ck] = tmp.obs[ck].astype(str).values
    log(f"CNV embedding: PCA {n_comps} dims -> Leiden {adata.obs[ck].nunique()} clusters")
