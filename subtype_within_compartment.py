#!/usr/bin/env python3
"""
subtype_within_compartment.py
=============================

Tests whether CNV-driven expression changes are breaking subtype
classification, then splits metacells into normal / CNV(malignant)
compartments and classifies subtypes within each compartment.

Usage:
    conda activate scrna-cnv
    python subtype_within_compartment.py --results-dir <sample>/seacell
    # also pass single-cell data when available -- needed to check whether
    # metacells are mixing lineages

This script runs four tests:

  Test A  Does CNV shift expression genome-wide?
          Compares each chromosome's mean log2FC (malignant vs normal)
          against that chromosome's CNV difference. Correlating cnv_score
          with PCs is avoided because it would be circular -- cnv_score is
          derived from the same expression matrix.

  Test B  Is the marker score confounded by depth (cell count / UMI)?
          Compares a raw z-mean score against a control-gene-matched
          score_genes score by their correlation with depth.

  Test C  Once split into two compartments, do subtypes separate when
          markers are compared within each one? Re-embeds within the
          compartment, clusters with leiden, and scores within-compartment;
          judged by whether exclusive-marker correlations turn negative and
          the top1/top2 score gap widens.

  Test D  Do metacells mix lineages in the first place? Assigns a lineage
          per single cell, then reports per-metacell composition and
          normalized entropy. High entropy means metacell-level subtype
          classification is not possible in principle, and metacells need
          to be rebuilt smaller within each compartment.

See README.md for the rationale behind each test.

Required inputs:
    metacell_obs.csv            putative_malignant / n_cells / cnv_score
    cnv_metacells.h5ad          layers['counts'] (or X), obsm['X_cnv'], uns['cnv']
    de_malignant_vs_normal.csv  used by Test A (Test A is skipped if absent)
    singlecells_qc.h5ad         used by Test D (skipped if absent)
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C_BLUE, C_ORANGE, C_AQUA, C_RED = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
C_GRAY, C_INK, C_INK2 = "#b9b8b2", "#0b0b0b", "#52514e"
SEQ_CMAP, DIV_CMAP = "Blues", "RdBu_r"

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 200, "savefig.bbox": "tight", "font.size": 9,
    "axes.edgecolor": C_INK2, "axes.labelcolor": C_INK, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.grid": True, "grid.color": "#e6e5e0",
    "grid.linewidth": 0.6, "axes.axisbelow": True, "xtick.color": C_INK2,
    "ytick.color": C_INK2, "legend.frameon": False,
})


def _setup_font() -> bool:
    import matplotlib.font_manager as fm
    for name in ("Noto Sans CJK JP", "IPAexGothic", "IPAGothic", "Hiragino Sans",
                 "Yu Gothic", "TakaoGothic", "VL PGothic"):
        try:
            fm.findfont(name, fallback_to_default=False)
        except Exception:
            continue
        plt.rcParams["font.family"] = name
        return True
    return False


HAS_CJK = _setup_font()


def T(ja: str, en: str) -> str:
    return ja if HAS_CJK else en


REPORT: list[str] = []


def say(line: str = "") -> None:
    print(line)
    REPORT.append(line)


def section(title: str) -> None:
    say("")
    say("=" * 72)
    say(title)
    say("=" * 72)


def _despine(ax) -> None:
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ---------------------------------------------------------------------------
# Marker panels
# ---------------------------------------------------------------------------
# It matters that these panels are mutually exclusive (non-overlapping).
# Mixing sets that share genes (e.g. Myeloid and Macrophage) would make
# "correlation between exclusive markers" impossible to measure.
DEFAULT_PANEL: dict[str, list[str]] = {
    "T/NK": ["CD3D", "CD3E", "CD3G", "CD2", "LCK", "ITK", "THEMIS", "SKAP1", "CD247",
             "IL7R", "CD28", "GZMA", "GZMB", "NKG7", "GNLY", "KLRD1", "PRF1", "EOMES"],
    "B/Plasma": ["CD79A", "CD79B", "MS4A1", "PAX5", "BANK1", "EBF1", "BLK", "CD22",
                 "FCRL1", "JCHAIN", "MZB1", "XBP1", "DERL3", "TNFRSF17", "SDC1"],
    "MonoMac": ["CD68", "CSF1R", "MRC1", "CD163", "C1QA", "C1QB", "C1QC", "MSR1",
                "AIF1", "TYROBP", "LYZ", "CD14", "ITGAM", "MARCO", "VSIG4", "FCER1G",
                "MNDA", "IRF8", "ZBTB46", "BATF3", "FLT3", "S100A8", "S100A9", "MMP9"],
    "Endothelial": ["PECAM1", "CDH5", "VWF", "KDR", "CLDN5", "EGFL7", "TEK", "ERG",
                    "FLT1", "ESAM"],
    "Fibro/Mural": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FBN1", "POSTN",
                    "THY1", "FAP", "PDGFRB", "ACTA2", "RGS5", "TAGLN", "NOTCH3", "MYH11"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "SFN", "KRT5", "KRT14"],
    "Mast": ["KIT", "CPA3", "MS4A2", "CMA1", "GATA2", "HDC"],
}

#: Sets representing *state* rather than lineage. Not used for lineage
#: assignment; reported for reference only.
STATE_PANEL: dict[str, list[str]] = {
    "Proliferating": ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "BIRC5", "TYMS", "RRM2"],
    "Hypoxia": ["VEGFA", "SLC2A1", "LDHA", "PGK1", "CA9", "ADM", "NDRG1"],
    "Interferon": ["ISG15", "IFI6", "MX1", "OAS1", "STAT1", "IRF7", "IFIT3"],
    "EMT": ["VIM", "FN1", "SNAI2", "ZEB1", "TWIST1", "CDH2", "SPARC"],
}


_NC_ACC = re.compile(r"^(?:chr)?(N[CTWZ]_\d+)(?:\.\d+)?$")
_NON_CHROM = {"chrM", "chrMT", "chrMito", "M", "MT"}


def load_chromosome_map(path: Path | None) -> dict[str, str]:
    """Map accession -> chromosome name, from an NCBI assembly report or a 2-column TSV."""
    if path is None or not Path(path).exists():
        return {}
    m: dict[str, str] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) >= 7:
                mol, acc = f[2].strip(), f[6].strip()
                if acc and acc not in ("na", "-"):
                    m[acc] = mol if mol not in ("na", "-", "") else f[0].strip()
            elif len(f) >= 2:
                m[f[0].strip()] = f[1].strip()
    return m


def chromosome_labels(names, chromosome_map: dict | None = None) -> tuple[dict, bool]:
    """Turn accessions into short chromosome names. Returns (label_map, inferred)."""
    chromosome_map = chromosome_map or {}
    lab, unresolved = {}, []
    for n in names:
        n = str(n)
        mt = _NC_ACC.match(n)
        key = n[3:] if n.startswith("chr") else n
        hit = chromosome_map.get(key) or (chromosome_map.get(mt.group(1)) if mt else None)
        if hit:
            lab[n] = str(hit)
        elif mt:
            unresolved.append(n)
        else:
            lab[n] = n.replace("chr", "", 1) or n
    chrom_like = [n for n in unresolved if _NC_ACC.match(n).group(1).startswith("NC_")]
    for n in unresolved:
        if n not in chrom_like:
            lab[n] = _NC_ACC.match(n).group(1)
    inferred = False
    if chrom_like:
        inferred = True
        order = sorted(chrom_like, key=lambda x: _NC_ACC.match(x).group(1))
        for i, n in enumerate(order):
            lab[n] = "X*" if (i == len(order) - 1 and len(order) >= 10) else f"{i + 1}*"
    return lab, inferred


def is_placed_chromosome(name: str) -> bool:
    """Excludes unplaced scaffolds (NW_/NT_) and mtDNA."""
    n = str(name)
    if n in _NON_CHROM:
        return False
    m = _NC_ACC.match(n)
    if m:
        return m.group(1).startswith("NC_")
    return True


def load_panel(path: Path | None) -> tuple[dict, dict]:
    """Use the JSON passed via --marker-panel, if given.

    Accepted shapes: {"lineage": {...}, "state": {...}}, or {"type": [gene, ...]}.
    """
    if path is None:
        return dict(DEFAULT_PANEL), dict(STATE_PANEL)
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    if "lineage" in d:
        return dict(d["lineage"]), dict(d.get("state", STATE_PANEL))
    return dict(d), dict(STATE_PANEL)


def restrict_panel(panel: dict, var_names) -> tuple[dict, pd.DataFrame]:
    """Restrict to genes present in the data; also return an availability table."""
    vs = set(map(str, var_names))
    out, rows = {}, []
    for k, gs in panel.items():
        have = [g for g in gs if g in vs]
        rows.append({"set": k, "n_total": len(gs), "n_found": len(have),
                     "missing": ",".join(g for g in gs if g not in vs) or "-"})
        if len(have) >= 2:
            out[k] = have
    return out, pd.DataFrame(rows).set_index("set")


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def _dense(x):
    import scipy.sparse as sp
    return np.asarray(x.todense()) if sp.issparse(x) else np.asarray(x, dtype=float)


def score_zmean(adata, panel: dict, center_rows: bool = False,
                min_detect: float = 0.10) -> pd.DataFrame:
    """Raw z-mean score, matching the pipeline's previous approach (for comparison).

    With center_rows=True, subtracts each observation's mean z across all
    genes -- a minimal correction for the "everything high / everything
    low" depth-and-complexity axis.
    """
    X = _dense(adata.X)
    det = (X > 0).mean(0) > min_detect
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    if center_rows:
        Z = Z - Z[:, det].mean(1, keepdims=True)
    cols = {}
    for k, gs in panel.items():
        idx = [adata.var_names.get_loc(g) for g in gs
               if g in adata.var_names and det[adata.var_names.get_loc(g)]]
        if idx:
            cols[k] = Z[:, idx].mean(1)
    return pd.DataFrame(cols, index=adata.obs_names)


def score_control_matched(adata, panel: dict, ctrl_size: int = 50) -> pd.DataFrame:
    """scanpy's score_genes (difference against expression-level-matched control genes).

    Control genes absorb whether an observation is globally high or low
    expressing, so depth and metacell-size effects are much smaller than
    with the z-mean.
    """
    import scanpy as sc

    q = adata.copy()
    cols = {}
    for k, gs in panel.items():
        g = [x for x in gs if x in q.var_names]
        if len(g) < 2:
            continue
        try:
            sc.tl.score_genes(q, g, score_name="_s", ctrl_size=ctrl_size, random_state=0)
            cols[k] = q.obs["_s"].to_numpy()
        except Exception as exc:
            say(f"  warning: score computation failed for {k}: {exc}")
    return pd.DataFrame(cols, index=adata.obs_names)


def assign_lineage(S: pd.DataFrame, min_gap: float = 0.05) -> pd.DataFrame:
    """Assign the top-scoring set; mark unassigned if the gap to runner-up is below min_gap."""
    if S.empty or S.shape[1] < 2:
        return pd.DataFrame({"lineage": "unassigned", "top": np.nan, "gap": np.nan},
                            index=S.index)
    top1 = S.max(1)
    top2 = S.apply(lambda r: r.nlargest(2).iloc[1], axis=1)
    gap = top1 - top2
    lin = S.idxmax(1).where(gap >= min_gap, "unassigned")
    return pd.DataFrame({"lineage": lin, "top": top1, "gap": gap})


def exclusivity(S: pd.DataFrame) -> tuple[float, float, pd.DataFrame]:
    """Correlation between exclusive marker sets; should be negative for a pure population."""
    C = S.corr(method="spearman")
    off = C.to_numpy()[np.triu_indices(C.shape[0], 1)]
    off = off[np.isfinite(off)]
    if off.size == 0:
        return np.nan, np.nan, C
    return float(off.mean()), float(np.mean(off < 0)), C


# ---------------------------------------------------------------------------
# Test A: does CNV shift expression at the chromosome level?
# ---------------------------------------------------------------------------

def chromosome_dosage_test(results: Path, cnv_adata, mc_obs: pd.DataFrame,
                           fig_dir: Path, chromosome_map: dict | None = None) -> dict:
    section(T("Test A: does CNV shift expression chromosome-wide?",
              "Test A: does CNV shift expression chromosome-wide?"))
    say("This is the only non-circular test here. cnv_score and X_cnv are both derived")
    say("from the expression matrix, so correlating PCs with cnv_score would be")
    say("circular. Instead this compares two independently estimated quantities:")
    say("the malignant-vs-normal DE log2FC, and that chromosome's CNV difference.")

    de_path = results / "de_malignant_vs_normal.csv"
    if not de_path.exists():
        say("Skipped: de_malignant_vs_normal.csv not found.")
        return {}
    de = pd.read_csv(de_path, index_col=0)
    if "X_cnv" not in getattr(cnv_adata, "obsm", {}) or "chromosome" not in cnv_adata.var:
        say("Skipped: obsm['X_cnv'] / var['chromosome'] not found.")
        return {}

    # Per-chromosome mean CNV (malignant - normal)
    import scipy.sparse as sp
    Xc = _dense(cnv_adata.obsm["X_cnv"])
    chr_pos = sorted(dict(cnv_adata.uns["cnv"]["chr_pos"]).items(), key=lambda kv: int(kv[1]))
    mal = (mc_obs["putative_malignant"].reindex(cnv_adata.obs_names) == "malignant").to_numpy()
    if mal.sum() < 3 or (~mal).sum() < 3:
        say("Skipped: fewer than 3 cells in malignant or normal group.")
        return {}
    dose = {}
    for i, (nm, s0) in enumerate(chr_pos):
        e0 = int(chr_pos[i + 1][1]) if i + 1 < len(chr_pos) else Xc.shape[1]
        if e0 - s0 < 1:
            continue
        if not is_placed_chromosome(nm):
            continue
        blk = Xc[:, s0:e0]
        dose[str(nm)] = float(blk[mal].mean() - blk[~mal].mean())
    d_chr = pd.Series(dose)
    labels, inferred = chromosome_labels(list(d_chr.index), chromosome_map)
    say("")
    say(f"Scope: {len(d_chr)} placed chromosomes (unplaced scaffolds and mtDNA excluded)")
    if inferred:
        say("  Note: chromosome names marked * in the figure are inferred from accession"
            " order; pass --chromosome-map for the real names.")

    g_chr = cnv_adata.var["chromosome"].astype(str)
    df = pd.DataFrame({"lfc": de["log2FoldChange"]}).join(
        pd.DataFrame({"chrom": g_chr})).dropna()
    df["dose"] = df["chrom"].map(d_chr)
    df = df.dropna(subset=["dose"])
    if df["chrom"].nunique() < 5:
        say("Skipped: fewer than 5 chromosomes.")
        return {}

    from scipy.stats import pearsonr, spearmanr
    r_gene = spearmanr(df["lfc"], df["dose"])
    r2_gene = float(pearsonr(df["dose"], df["lfc"]).statistic ** 2)
    g = df.groupby("chrom").agg(mean_lfc=("lfc", "mean"), n=("lfc", "size"))
    g["dose"] = d_chr.reindex(g.index)
    g = g.dropna()
    r_chr = spearmanr(g["mean_lfc"], g["dose"])

    up = d_chr[d_chr > d_chr.quantile(0.75)].index
    dn = d_chr[d_chr < d_chr.quantile(0.25)].index
    s_up, s_dn = df[df["chrom"].isin(up)], df[df["chrom"].isin(dn)]

    say("")
    say(f"Per-gene: Spearman rho = {r_gene.statistic:+.3f} (p={r_gene.pvalue:.2g}), "
        f"R^2 = {r2_gene:.4f}")
    say(f"Per-chromosome: Spearman rho = {r_chr.statistic:+.3f} (p={r_chr.pvalue:.2g}), "
        f"n = {len(g)} chromosomes")
    say(f"CNV-gained chromosomes ({len(up)}, {len(s_up):,} genes): "
        f"mean log2FC {s_up['lfc'].mean():+.3f} / fraction positive {(s_up['lfc'] > 0).mean():.1%}")
    say(f"CNV-lost chromosomes ({len(dn)}, {len(s_dn):,} genes): "
        f"mean log2FC {s_dn['lfc'].mean():+.3f} / fraction positive {(s_dn['lfc'] > 0).mean():.1%}")
    say(f"All genes:                              "
        f"mean log2FC {df['lfc'].mean():+.3f} / fraction positive {(df['lfc'] > 0).mean():.1%}")
    say("")
    if r_chr.statistic > 0.6 and r_chr.pvalue < 0.01:
        say("→ Verdict: CNV consistently shifts expression chromosome-wide (dosage effect present).")
    else:
        say("→ Verdict: no consistent chromosome-level dosage effect detected.")
    if r2_gene < 0.05:
        say(f"→ However, only {r2_gene:.1%} of per-gene variance is explained.")
        say("   The dosage shift is small compared to typical lineage-marker expression")
        say("   differences (2-4 log2), so \"CNV broke marker calling\" is not a likely explanation.")
    gshift = float(df["lfc"].mean())
    if abs(gshift) > 0.2:
        say("")
        say(f"Note: the mean log2FC across all genes is {gshift:+.3f}, with "
            f"{(df['lfc'] > 0).mean():.0%} positive. This is a global offset unrelated to")
        say("  chromosome position, caused by differences in metacell cell count and")
        say("  per-cell UMI (post-normalization complexity); it's larger than the dosage"
            " effect and more harmful -- see Test B.")

    # --- Figure: chromosome-level relationship ---
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0))
    ax = axes[0]
    ax.scatter(g["dose"], g["mean_lfc"], s=np.clip(g["n"] / 8, 8, 90),
               color=C_BLUE, alpha=0.8, linewidth=0)
    ext = pd.concat([g.nlargest(4, "dose"), g.nsmallest(4, "dose")])
    for k, (nm, row) in enumerate(ext.iterrows()):
        ax.annotate(labels.get(str(nm), str(nm)), (row["dose"], row["mean_lfc"]),
                    fontsize=7.5, fontweight="bold", color=C_INK,
                    xytext=(6, 6 if k % 2 == 0 else -12), textcoords="offset points",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=C_INK2))
    if len(g) > 2:
        sl = np.polyfit(g["dose"], g["mean_lfc"], 1)
        xs = np.linspace(g["dose"].min(), g["dose"].max(), 10)
        ax.plot(xs, np.polyval(sl, xs), color=C_RED, lw=1.2)
    ax.axvline(0, color=C_INK2, lw=0.8, ls=":")
    ax.set_xlabel(T("chromosome CNV difference (malignant - normal)", "chromosome CNV difference (malignant - normal)"))
    ax.set_ylabel(T("mean log2FC of genes on that chromosome", "mean log2FC of genes on that chromosome"))
    ax.set_title(T(f"Test A  chromosome-level rho = {r_chr.statistic:+.3f}",
                   f"Test A  chromosome-level rho = {r_chr.statistic:+.3f}"))
    _despine(ax)

    ax = axes[1]
    bins = np.linspace(np.nanpercentile(df["lfc"], 0.5), np.nanpercentile(df["lfc"], 99.5), 60)
    ax.hist(s_dn["lfc"], bins=bins, color=C_BLUE, alpha=0.65, density=True,
            label=T(f"lost chromosomes ({len(dn)})", f"lost chromosomes ({len(dn)})"))
    ax.hist(s_up["lfc"], bins=bins, color=C_RED, alpha=0.65, density=True,
            label=T(f"gained chromosomes ({len(up)})", f"gained chromosomes ({len(up)})"))
    ax.axvline(0, color=C_INK2, lw=0.9)
    ax.axvline(df["lfc"].mean(), color=C_INK, lw=1.2, ls="--",
               label=T(f"mean over all genes {df['lfc'].mean():+.2f}",
                       f"mean over all genes {df['lfc'].mean():+.2f}"))
    ax.set_xlabel("log2 fold change (malignant / normal)")
    ax.set_ylabel(T("density", "density"))
    ax.set_title(T("the dosage effect appears as a shift of the distribution",
                   "the dosage effect appears as a shift of the distribution"))
    ax.legend(fontsize=7.5)
    _despine(ax)
    fig.suptitle(T("Fig S1  CNV dosage effect and the global offset",
                   "Fig S1  CNV dosage effect and the global offset"),
                 y=1.03, fontsize=11, fontweight="bold")
    fig.savefig(fig_dir / "S1_dosage_test.png")
    plt.close(fig)
    say(f"→ figures/S1_dosage_test.png")
    g.insert(0, "chromosome", [labels.get(str(i), str(i)) for i in g.index])
    g.round(4).to_csv(results / "chromosome_dosage_vs_lfc.csv")
    say("→ chromosome_dosage_vs_lfc.csv")
    return {"rho_gene": float(r_gene.statistic), "r2_gene": r2_gene,
            "rho_chrom": float(r_chr.statistic), "global_shift": gshift}


# ---------------------------------------------------------------------------
# Test B: is the marker score confounded by depth?
# ---------------------------------------------------------------------------

def depth_confound_test(mc, panel: dict, mc_obs: pd.DataFrame, fig_dir: Path) -> dict:
    section(T("Test B: is the marker score driven by depth (cell count / UMI)?",
              "Test B: is the marker score driven by depth?"))
    from scipy.stats import spearmanr

    depth = None
    for key in ("n_cells", "total_counts", "n_counts"):
        if key in mc_obs:
            depth = mc_obs[key].reindex(mc.obs_names).to_numpy(dtype=float)
            depth_name = key
            break
    if depth is None:
        depth = np.asarray(_dense(mc.X).sum(1), dtype=float)
        depth_name = "sum(X)"

    methods = {
        "zmean": (T("raw z-mean (previous)", "raw z-mean (previous)"), score_zmean(mc, panel)),
        "zmean_centered": (T("z-mean + row centering", "z-mean + row centering"),
                           score_zmean(mc, panel, True)),
        "score_genes": (T("score_genes (control-matched)", "score_genes (control-matched)"),
                        score_control_matched(mc, panel)),
    }
    say("")
    say(f"Depth metric: {depth_name}")
    say(f"{'method':<34} {'depth corr':>12} {'exclusivity':>10} {'frac neg':>9} {'median gap':>13} {'n types':>5}")
    out = {}
    for key, (nm, S) in methods.items():
        if S.empty:
            continue
        A = assign_lineage(S)
        rd = abs(spearmanr(S.max(1), depth).statistic)
        mean_off, frac_neg, _ = exclusivity(S)
        say(f"{nm:<34} {rd:12.3f} {mean_off:+10.3f} {frac_neg:9.0%} "
            f"{A['gap'].median():13.3f} {A['lineage'].nunique():5d}")
        out[key] = {"label": nm, "depth_corr": rd, "mean_exclusivity": mean_off,
                    "frac_negative": frac_neg, "median_gap": float(A["gap"].median())}
    say("")
    say("A method with a high correlation to depth is mostly capturing \"deeper metacells")
    say("score higher on every marker,\" so the types it assigns are nearly meaningless.")
    say("If score_genes or row-centering lowers this correlation, the previous score was"
        " artifactual to that same degree.")

    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    names = list(out)
    ax.barh(range(len(names)), [out[n]["depth_corr"] for n in names],
            color=[C_RED if out[n]["depth_corr"] > 0.25 else C_BLUE for n in names], height=0.55)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([out[n]["label"] for n in names], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0.25, color=C_INK2, ls=":", lw=0.9)
    ax.set_xlabel(T("|Spearman rho|: max marker score vs depth",
                    "|Spearman rho|: max marker score vs depth"))
    ax.set_title(T("Fig S2  depth confounding", "Fig S2  depth confounding"))
    _despine(ax)
    fig.savefig(fig_dir / "S2_depth_confound.png")
    plt.close(fig)
    say("→ figures/S2_depth_confound.png")
    return out


# ---------------------------------------------------------------------------
# Test C: subtype classification within each compartment
# ---------------------------------------------------------------------------

def cluster_within(adata, n_hvg: int, n_pcs: int, resolution: float,
                   restrict_to: list[str] | None = None):
    """Within a compartment: HVG -> PCA -> kNN -> leiden. If restrict_to is
    given, use only that gene space."""
    import scanpy as sc

    t = adata.copy()
    # If a highly_variable column is still present in var, scanpy's pca()
    # picks it up automatically via mask_var='highly_variable', which
    # defeats an explicit gene-space restriction. Drop it first.
    t.var = t.var.drop(columns=[c for c in ("highly_variable",) if c in t.var.columns])
    if restrict_to:
        keep = [g for g in restrict_to if g in t.var_names]
        if len(keep) < 10:
            return None, None
        t = t[:, keep].copy()
    else:
        sc.pp.highly_variable_genes(t, n_top_genes=min(n_hvg, t.n_vars - 1))
        t = t[:, t.var.highly_variable].copy()
        t.var = t.var.drop(columns=["highly_variable"])
    sc.pp.scale(t, max_value=10)
    npc = int(min(n_pcs, t.n_vars - 1, t.n_obs - 1))
    if npc < 2:
        return None, None
    try:
        sc.tl.pca(t, n_comps=npc, mask_var=None)
    except Exception:
        npc = int(min(npc, min(t.n_obs, t.n_vars) - 1))
        if npc < 2:
            return None, None
        sc.tl.pca(t, n_comps=npc, mask_var=None, svd_solver="randomized")
    npc = int(min(npc, t.obsm["X_pca"].shape[1]))
    sc.pp.neighbors(t, n_neighbors=int(min(15, max(2, t.n_obs - 1))), n_pcs=npc)
    sc.tl.leiden(t, resolution=resolution, key_added="cl",
                 flavor="igraph", n_iterations=2, directed=False)
    return t.obs["cl"].astype(str).to_numpy(), t.obsm["X_pca"]



def _relative_threshold_check(cluster_scores: pd.DataFrame, best_label: pd.Series,
                              sizes: pd.Series, min_zscore: float = 1.0) -> None:
    """Apply the same relative threshold the pipeline's annotate_coarse_celltype
    uses ("cluster-to-cluster z >= min_zscore") and report how many clusters
    end up with a type.

    Because the threshold is relative, only roughly the top fifth of k
    clusters can structurally satisfy z >= 1, regardless of the data. Even
    when every cluster is the same lineage with a uniformly high score,
    all but the top 1-2 fall to "Other" -- this can be a major cause of
    annotation failures on its own.
    """
    cs = cluster_scores
    sd = cs.std(ddof=0).replace(0, np.nan)
    z = ((cs - cs.mean()) / sd).fillna(0.0)
    best = cs.max(1)
    bz = pd.Series([z.loc[i, l] for i, l in best_label.items()], index=cs.index)
    keep = (best > 0) & (bz >= min_zscore)
    n_keep = int(keep.sum())
    n_obs_keep = int(sizes.reindex(cs.index).fillna(0)[keep.to_numpy()].sum())
    say(f"    Applying the pipeline's relative threshold (score>0 and cluster-to-cluster"
        f" z >= {min_zscore}): only {n_keep}/{len(cs)} clusters get a type "
        f"({n_obs_keep}/{int(sizes.sum())} metacells). The rest fall to Other.")
    if n_keep < len(cs) and best_label.nunique() == 1:
        say(f"    Note: all {len(cs)} clusters are the same type ({best_label.iloc[0]}), "
            "with uniformly positive scores, yet most still fall to Other under the"
            " relative threshold.")
        say("       This isn't a missing-marker or CNV-broken-expression problem --"
            " it's the calling rule itself.")
        say(f"       Fix: lower min_zscore, use --min-score for an absolute threshold, or "
            "inspect the score table and set normal_clusters by hand.")


def subtype_within_compartments(mc, panel: dict, state_panel: dict, mc_obs: pd.DataFrame,
                                results: Path, fig_dir: Path, n_hvg: int, n_pcs: int,
                                resolution: float, min_gap: float,
                                compartment_key: str) -> pd.DataFrame:
    section(T("Test C: subtype classification within each compartment",
              "Test C: subtype classification within each compartment"))
    comp = mc_obs[compartment_key].astype(str).reindex(mc.obs_names)
    say(f"Compartments ({compartment_key}): {comp.value_counts().to_dict()}")

    lin_genes = sorted({g for v in panel.values() for g in v})
    rows, per_obs = [], []
    for name in sorted(comp.dropna().unique()):
        sub = mc[(comp == name).to_numpy()].copy()
        if sub.n_obs < 6:
            say(f"\n--- {name} (n={sub.n_obs}) --- skipped, too few")
            continue
        say("")
        say(f"--- {name} compartment (n={sub.n_obs}) ---")
        S = score_control_matched(sub, panel)              # within-compartment score (depth-corrected)
        A = assign_lineage(S, min_gap=min_gap)
        mean_off, frac_neg, C = exclusivity(S)
        say(f"Exclusive-marker correlation: mean {mean_off:+.3f} / fraction negative {frac_neg:.0%}")
        say(f"Top1-Top2 gap: median {A['gap'].median():.3f} / "
            f"fraction with gap >= {min_gap}: {(A['gap'] >= min_gap).mean():.0%}")
        say(f"Assignment: {A['lineage'].value_counts().to_dict()}")

        for space, restrict in ((T("HVG space", "HVG space"), None),
                                (T("lineage-marker space", "lineage-marker space"), lin_genes)):
            lab, _ = cluster_within(sub, n_hvg, n_pcs, resolution, restrict)
            if lab is None:
                continue
            m = S.groupby(lab).mean()
            asg = m.idxmax(1)
            gap = m.max(1) - m.apply(lambda r: r.nlargest(2).iloc[1], axis=1)
            say(f"  {space}: {len(set(lab))} clusters → {asg.nunique()} types "
                f"{asg.value_counts().to_dict()}")
            say(f"    Cluster-mean-score top1-top2 gap: median {gap.median():.3f} / "
                f"max {gap.max():.3f}")
            rows.append({"compartment": name, "space": space,
                         "n_clusters": len(set(lab)), "n_types": int(asg.nunique()),
                         "median_gap": float(gap.median()), "max_gap": float(gap.max())})
            if restrict is None:
                sub.obs["cluster_within"] = lab
                _relative_threshold_check(
                    m, asg, pd.Series(lab).value_counts().sort_index(), min_zscore=1.0)

        SS = score_control_matched(sub, state_panel) if state_panel else pd.DataFrame()
        if not SS.empty:
            top = SS.mean().sort_values(ascending=False)
            say(f"  Compartment-mean state scores (reference): "
                + ", ".join(f"{k}={v:+.3f}" for k, v in top.items()))

        d = A.copy()
        d["compartment"] = name
        d["cluster_within"] = sub.obs.get("cluster_within", pd.Series(index=sub.obs_names))
        for k in S.columns:
            d[f"score_{k}"] = S[k]
        per_obs.append(d)

        # --- Figure: cluster x marker-set score ---
        if "cluster_within" in sub.obs:
            m = S.groupby(sub.obs["cluster_within"].to_numpy()).mean()
            fig, ax = plt.subplots(figsize=(1.0 + 0.62 * m.shape[1], 0.9 + 0.34 * m.shape[0]))
            v = float(np.nanmax(np.abs(m.to_numpy()))) or 1e-3
            im = ax.imshow(m.to_numpy(), cmap=DIV_CMAP, vmin=-v, vmax=v, aspect="auto")
            ax.set_xticks(range(m.shape[1]))
            ax.set_xticklabels(m.columns, rotation=45, ha="right", fontsize=8)
            ax.set_yticks(range(m.shape[0]))
            ax.set_yticklabels([f"cl{i}" for i in m.index], fontsize=8)
            for i in range(m.shape[0]):
                j = int(np.argmax(m.to_numpy()[i]))
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                           edgecolor=C_INK, lw=1.4))
            ax.grid(False)
            fig.colorbar(im, ax=ax, label=T("within-compartment score", "within-compartment score"), shrink=0.8)
            ax.set_title(T(f"Fig S3  marker scores of clusters within {name}",
                           f"Fig S3  marker scores of clusters within {name}"), fontsize=9.5)
            fig.savefig(fig_dir / f"S3_marker_scores_{name}.png")
            plt.close(fig)
            say(f"  → figures/S3_marker_scores_{name}.png")

    if rows:
        say("")
        say("Summary:")
        for line in pd.DataFrame(rows).round(3).to_string(index=False).split("\n"):
            say("  " + line)
    if per_obs:
        out = pd.concat(per_obs)
        out.round(4).to_csv(results / "metacell_subtype_within_compartment.csv")
        say("→ metacell_subtype_within_compartment.csv")
        return out
    return pd.DataFrame()


# ---------------------------------------------------------------------------
# Test D: do metacells mix lineages? (requires single-cell data)
# ---------------------------------------------------------------------------

def metacell_purity_test(sc_adata, panel: dict, mc_obs: pd.DataFrame, results: Path,
                         fig_dir: Path, min_gap: float, metacell_key: str = "SEACell",
                         min_cells: int = 20) -> dict:
    section(T("Test D: do metacells mix lineages?", "Test D: do metacells mix lineages?"))
    import scanpy as sc

    if metacell_key not in sc_adata.obs:
        say(f"Skipped: obs['{metacell_key}'] not found.")
        return {}
    s = sc_adata.copy()
    if "counts" in getattr(s, "layers", {}):
        s.X = s.layers["counts"].copy()
        s.uns.pop("log1p", None)
        sc.pp.normalize_total(s, target_sum=1e4)
        sc.pp.log1p(s)
    umi = np.asarray(_dense(s.layers["counts"]).sum(1)) if "counts" in s.layers else None
    if umi is not None:
        say(f"Single cells n={s.n_obs:,} / median UMI {np.median(umi):,.0f}")
        if np.median(umi) < 800:
            say("Note: median UMI per cell is below 800. Per-cell marker calling becomes")
            say("  unstable due to dropout; the composition estimates below may be underestimates.")
    else:
        say(f"Single cells n={s.n_obs:,}")

    S = score_control_matched(s, panel)
    A = assign_lineage(S, min_gap=min_gap)
    say(f"Cell lineage assignment: {A['lineage'].value_counts().to_dict()}")
    mean_off, frac_neg, _ = exclusivity(S)
    say(f"Cell-level exclusive-marker correlation: mean {mean_off:+.3f} / fraction negative {frac_neg:.0%}")

    d = pd.DataFrame({"mc": s.obs[metacell_key].astype(str).to_numpy(),
                      "lin": A["lineage"].to_numpy()}, index=s.obs_names)
    d = d[d["lin"] != "unassigned"]
    if d.empty:
        say("Skipped: no assigned cells.")
        return {}
    ct = pd.crosstab(d["mc"], d["lin"])
    frac = ct.div(ct.sum(1), axis=0)
    big = ct.sum(1) >= min_cells
    f = frac.replace(0, np.nan)
    H = -(f * np.log(f)).sum(1) / np.log(max(ct.shape[1], 2))
    dom = frac.max(1)

    say("")
    say(f"Metacells with >= {min_cells} assigned cells: {int(big.sum())}")
    say(f"Dominant-lineage occupancy: median {dom[big].median():.1%} / "
        f">70% in {np.mean(dom[big] > 0.7):.0%} / >90% in {np.mean(dom[big] > 0.9):.0%}")
    say(f"Normalized entropy (0=single lineage, 1=fully mixed): median {H[big].median():.2f}")
    say("")
    if H[big].median() > 0.4:
        say("→ Verdict: metacells mix lineages. Metacell-level subtype classification is")
        say("   not possible in principle; rebuild smaller metacells within each compartment")
        say("   (e.g. drop cells-per-metacell to 20-30 for the normal group).")
    else:
        say("→ Verdict: metacells are largely single-lineage. Metacell-level subtype")
        say("   analysis is reasonable.")

    comp = mc_obs["putative_malignant"].astype(str) if "putative_malignant" in mc_obs else None
    out = pd.DataFrame({"n_assigned": ct.sum(1), "dominant": frac.idxmax(1),
                        "dominant_frac": dom.round(3), "entropy": H.round(3)}).join(
        frac.round(3).add_prefix("frac_"))
    if comp is not None:
        out.insert(0, "compartment", comp.reindex(out.index))
    out.to_csv(results / "metacell_lineage_composition.csv")
    say("→ metacell_lineage_composition.csv")

    # --- Figure ---
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.9),
                             gridspec_kw={"width_ratios": [1.55, 1]})
    ax = axes[0]
    sel = frac[big].copy()
    order = sel.idxmax(1).astype(str) + "_" + (1 - sel.max(1)).round(3).astype(str)
    sel = sel.loc[order.sort_values().index]
    bottom = np.zeros(len(sel))
    cols = [C_BLUE, C_ORANGE, C_AQUA, C_RED, "#8e6bbf", "#c9a227", "#4c8f8b", "#b0567a"]
    for i, c in enumerate(sel.columns):
        ax.bar(range(len(sel)), sel[c].to_numpy(), bottom=bottom, width=1.0,
               color=cols[i % len(cols)], label=str(c), linewidth=0)
        bottom += sel[c].to_numpy()
    ax.set_xlim(-0.5, len(sel) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_ylabel(T("lineage composition", "lineage composition"))
    ax.set_xlabel(T(f"metacells (n={len(sel)}, grouped by dominant lineage)",
                    f"metacells (n={len(sel)}, grouped by dominant lineage)"))
    ax.legend(fontsize=7, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.42))
    ax.grid(False)
    ax.set_title(T("lineage composition of each metacell", "lineage composition of each metacell"))

    ax = axes[1]
    ax.hist(H[big].dropna(), bins=20, color=C_BLUE, alpha=0.85)
    ax.axvline(float(H[big].median()), color=C_RED, lw=1.4,
               label=T(f"median {H[big].median():.2f}", f"median {H[big].median():.2f}"))
    ax.axvline(0.4, color=C_INK2, lw=0.9, ls=":",
               label=T("mixing threshold 0.4", "mixing threshold 0.4"))
    ax.set_xlabel(T("normalised entropy", "normalised entropy"))
    ax.set_ylabel(T("metacells", "metacells"))
    ax.set_title(T("close to 0 = single lineage", "close to 0 = single lineage"))
    ax.legend(fontsize=7.5)
    _despine(ax)
    fig.suptitle(T("Fig S4  lineage purity of metacells", "Fig S4  lineage purity of metacells"),
                 y=1.04, fontsize=11, fontweight="bold")
    fig.savefig(fig_dir / "S4_metacell_purity.png")
    plt.close(fig)
    say("→ figures/S4_metacell_purity.png")

    A.assign(metacell=s.obs[metacell_key].astype(str)).round(4).to_csv(
        results / "cell_lineage_assignment.csv")
    say("→ cell_lineage_assignment.csv")
    return {"median_entropy": float(H[big].median()),
            "median_dominant": float(dom[big].median()),
            "n_metacells": int(big.sum())}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Subtype classification within CNV-based compartments, and tests of its assumptions")
    p.add_argument("--results-dir", default="results")
    p.add_argument("--out-dir", default=None, help="Figure output directory (default: <results>/figures)")
    p.add_argument("--marker-panel", default=None,
                   help="Marker panel JSON. {\"lineage\":{...},\"state\":{...}} or {type:[genes]}")
    p.add_argument("--compartment-key", default="putative_malignant",
                   help="Column in metacell_obs.csv used to split into two compartments (default putative_malignant)")
    p.add_argument("--min-gap", type=float, default=0.05,
                   help="Mark unassigned if the top1-top2 score gap is below this")
    p.add_argument("--resolution", type=float, default=1.0, help="Leiden resolution used within each compartment")
    p.add_argument("--n-hvg", type=int, default=2000)
    p.add_argument("--n-pcs", type=int, default=20)
    p.add_argument("--min-cells-per-metacell", type=int, default=20,
                   help="Minimum assigned cells for a metacell to be included in Test D's composition")
    p.add_argument("--chromosome-map", default=None,
                   help="NCBI assembly report; resolves real chromosome names in Fig S1")
    p.add_argument("--skip-singlecell", action="store_true",
                   help="Do not read singlecells_qc.h5ad (skip Test D)")
    args = p.parse_args(argv)

    results = Path(args.results_dir)
    if not results.exists():
        print(f"Directory does not exist: {results}", file=sys.stderr)
        return 1
    fig_dir = Path(args.out_dir) if args.out_dir else results / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    import scanpy as sc

    say("CNV compartments x within-compartment subtype classification  (subtype_within_compartment.py 1.0)")
    say(f"Input: {results.resolve()}")
    say(f"Figures: {fig_dir.resolve()}")

    obs_path = results / "metacell_obs.csv"
    if not obs_path.exists():
        print("metacell_obs.csv not found", file=sys.stderr)
        return 1
    mc_obs = pd.read_csv(obs_path, index_col=0)
    if args.compartment_key not in mc_obs:
        print(f"metacell_obs.csv has no '{args.compartment_key}' column", file=sys.stderr)
        return 1

    cnv_adata = None
    for nm in ("cnv_metacells.h5ad", "metacells.h5ad"):
        if (results / nm).exists():
            cnv_adata = ad.read_h5ad(results / nm)
            say(f"metacell: {nm} {cnv_adata.shape}")
            break
    if cnv_adata is None:
        print("cnv_metacells.h5ad / metacells.h5ad not found", file=sys.stderr)
        return 1

    # Build a log1p(CPM) metacell matrix
    mc = cnv_adata.copy()
    if "counts" in getattr(cnv_adata, "layers", {}):
        mc.X = cnv_adata.layers["counts"].copy()
        mc.uns.pop("log1p", None)
        sc.pp.normalize_total(mc, target_sum=1e4)
        sc.pp.log1p(mc)
        say("Expression: normalized layers['counts'] to CPM and log1p")
    else:
        say("Note: layers['counts'] not found; using X as-is (assumed already normalized)")
    for k in ("n_cells", "total_counts"):
        if k in mc_obs:
            mc.obs[k] = mc_obs[k].reindex(mc.obs_names).values

    panel_raw, state_raw = load_panel(Path(args.marker_panel) if args.marker_panel else None)
    panel, avail = restrict_panel(panel_raw, mc.var_names)
    state_panel, _ = restrict_panel(state_raw, mc.var_names)
    section(T("marker panel availability", "marker panel availability"))
    for line in avail.to_string().split("\n"):
        say("  " + line)
    if len(panel) < 2:
        say("Note: fewer than 2 usable marker sets. Match gene names via --marker-panel.")
        return 1
    say("")
    say("If markers are present but no type gets assigned, the cause is usually not")
    say("gene names, but scoring (Test B) or metacell lineage mixing (Test D).")

    chrom_map = load_chromosome_map(Path(args.chromosome_map) if args.chromosome_map else None)
    resA = chromosome_dosage_test(results, cnv_adata, mc_obs, fig_dir, chrom_map)
    resB = depth_confound_test(mc, panel, mc_obs, fig_dir)
    subtype_within_compartments(mc, panel, state_panel, mc_obs, results, fig_dir,
                                args.n_hvg, args.n_pcs, args.resolution,
                                args.min_gap, args.compartment_key)

    resD = {}
    sc_path = results / "singlecells_qc.h5ad"
    if not args.skip_singlecell and sc_path.exists():
        try:
            sc_adata = ad.read_h5ad(sc_path)
            sc_panel, _ = restrict_panel(panel_raw, sc_adata.var_names)
            resD = metacell_purity_test(sc_adata, sc_panel, mc_obs, results, fig_dir,
                                        args.min_gap,
                                        min_cells=args.min_cells_per_metacell)
        except Exception as exc:
            say(f"Note: Test D failed: {exc}")
            plt.close("all")
    elif not args.skip_singlecell:
        say("")
        say("Note: singlecells_qc.h5ad not found; cannot run Test D (metacell lineage purity).")
        say("  Whether metacells mix lineages can only be determined by this test, so be")
        say("  sure to run it wherever a single-cell file is available.")

    # --- Overall verdict ---
    section(T("overall verdict", "overall verdict"))
    if resA:
        if resA["rho_chrom"] > 0.6 and resA["r2_gene"] < 0.05:
            say("1. CNV's dosage effect is chromosome-level consistent "
                f"(rho={resA['rho_chrom']:+.2f}), but explains only "
                f"{resA['r2_gene']:.1%} of per-gene variance.")
            say("   → \"CNV broke marker calling\" is not the main cause.")
        elif resA["rho_chrom"] > 0.6:
            say(f"1. CNV's dosage effect is strong (rho={resA['rho_chrom']:+.2f}, "
                f"R^2={resA['r2_gene']:.1%}). Worth considering a dosage correction.")
        else:
            say("1. No chromosome-level dosage effect was detected.")
        if abs(resA["global_shift"]) > 0.2:
            say(f"   Note: the mean log2FC across all genes is {resA['global_shift']:+.2f}, "
                "a global offset unrelated to chromosome position, caused by depth differences.")
    if resB:
        worst = max(resB, key=lambda k: resB[k]["depth_corr"])
        best = min(resB, key=lambda k: resB[k]["depth_corr"])
        say(f"2. Marker-score-vs-depth correlation: worst {resB[worst]['label']} = "
            f"{resB[worst]['depth_corr']:.2f} / best {resB[best]['label']} = "
            f"{resB[best]['depth_corr']:.2f}.")
        cur = resB.get("zmean")
        if best != "zmean" and cur and cur["depth_corr"] - resB[best]["depth_corr"] > 0.15:
            say(f"   → Switching scoring to \"{resB[best]['label']}\" is the top-priority fix. "
                "As long as the raw z-mean is used, some assigned types simply reflect depth.")
        elif best == "zmean":
            say("   → The plain z-mean is the least depth-confounded here; replacing the "
                "scoring method is not the top priority for this sample.")
        if resB[best]["depth_corr"] > 0.25:
            say(f"   Note: even the best method still has a depth correlation of "
                f"{resB[best]['depth_corr']:.2f}. Metacell size may be strongly skewed "
                "(rebuild with matched cells-per-metacell across compartments).")
        sg = resB.get("score_genes")
        if sg and sg["frac_negative"] > (cur or sg)["frac_negative"]:
            say(f"   Also, the fraction of exclusive-marker pairs with negative correlation "
                f"improves from {(cur or sg)['frac_negative']:.0%} to {sg['frac_negative']:.0%} "
                "with score_genes -- weigh this alongside the depth correlation when choosing.")
    if resD:
        if resD["median_entropy"] > 0.4:
            say(f"3. Median metacell lineage entropy {resD['median_entropy']:.2f}: lineages are")
            say("   mixed, so subtypes won't emerge without rebuilding smaller metacells per compartment.")
        else:
            say(f"3. Median metacell lineage entropy {resD['median_entropy']:.2f}, dominant-lineage "
                f"occupancy {resD['median_dominant']:.0%}. Metacells are usable as-is.")
    else:
        say("3. Metacell lineage purity not verified (requires single-cell data).")

    out = results / "subtype_report.txt"
    out.write_text("\n".join(REPORT) + "\n", encoding="utf-8")
    print()
    print(f"[subtype] report: {out}")
    print(f"[subtype] figures: {sorted(q.name for q in fig_dir.glob('S*.png'))}")
    return 0


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    sys.exit(main())
