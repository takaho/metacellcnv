#!/usr/bin/env python3
"""
Interpret metacellcnv.py output and generate a report plus figures
==============================================================

Version: 1.4 (expects the results/ directory produced by pipeline v2.3)

Usage:
    conda activate scrna-cnv
    python metacellcnv_visualize.py --results-dir results
    # Figures go to results/figures/, the summary to results/interpretation_report.txt

Outputs:
    figures/01_qc.png                QC flag breakdown and distributions
    figures/02_metacell_quality.png  compactness x separation (2D quality view)
    figures/03_cnv_heatmap_*.png     Chromosome-ordered CNV heatmap, by group
    figures/04_cnv_score.png         cnv_score distribution per clone and the malignant/normal boundary
    figures/05_cnv_umap.png          UMAP of CNV space (small multiples per clone)
    figures/06_de_volcano.png        DE volcano / MA plots
    figures/07_mito_balance.png      Mitochondrial/nuclear expression balance and anomalous metacells
    figures/08_cnv_chromosome_heatmap.png
                                     Chromosome x metacell gain/loss heatmap. Metacells are
                                     ordered by hierarchical clustering of the CNV profile, not by ID
    figures/09_deg_swarm.png         Top-N DEG metacell expression as a swarm plot, by cell type group
    cnv_chromosome_means.csv         Values behind figure 8 (mean CNV per chromosome and cluster number)
    de_swarm_top_genes.csv           Genes plotted in figure 9 and their per-group medians
    de_top_genes.csv                 Top list of significant genes
    interpretation_report.txt        Numeric summary and interpretation notes

Design rationale (visualization):
- CNV is polarity data (gain / neutral / loss), so a diverging colormap (2 colors +
  neutral gray) is used; rainbow colormaps are avoided.
- Magnitude data such as cnv_score uses a single-hue sequential colormap.
- Identity data such as clone labels is shown with small multiples and direct labels
  rather than relying on color alone in scatter plots (for color-vision accessibility
  and black-and-white printing).
- Dual-axis plots (different scales on left/right) are not used.

See README.md for design rationale and validation results.
"""

from __future__ import annotations

import argparse
import re
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- Color palette (validated) ----------------------------------------------
C_BLUE = "#2a78d6"      # categorical slot 1 / primary hue for sequential
C_ORANGE = "#eb6834"
C_AQUA = "#1baf7a"
C_RED = "#e34948"
C_GRAY = "#b9b8b2"
C_INK = "#0b0b0b"
C_INK2 = "#52514e"
SEQ_CMAP = "Blues"       # single hue
DIV_CMAP = "RdBu_r"      # two hues + neutral (CNV gain/loss)

plt.rcParams.update(
    {
        "figure.dpi": 130,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.edgecolor": C_INK2,
        "axes.labelcolor": C_INK,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.grid": True,
        "grid.color": "#e6e5e0",
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "xtick.color": C_INK2,
        "ytick.color": C_INK2,
        "legend.frameon": False,
    }
)

def _setup_font() -> bool:
    """Use a CJK font if one is available; fall back to English figure labels otherwise.

    Remote Linux servers often lack CJK fonts, in which case CJK labels would render
    as tofu boxes and make the figure unreadable. Detect font availability so the
    label language can adapt to the environment.
    """
    import matplotlib.font_manager as fm

    available = {f.name for f in fm.fontManager.ttflist}
    for name in (
        "Hiragino Sans", "Hiragino Kaku Gothic ProN", "Yu Gothic", "YuGothic",
        "Noto Sans CJK JP", "Noto Sans JP", "Source Han Sans JP",
        "IPAexGothic", "IPAGothic", "TakaoGothic", "VL Gothic", "MS Gothic",
    ):
        if name in available:
            plt.rcParams["font.family"] = [name, "DejaVu Sans"]
            plt.rcParams["axes.unicode_minus"] = False
            return True
    return False


HAS_CJK = _setup_font()


def T(ja: str, en: str) -> str:
    """Figure label text. Falls back to English when no CJK font is available."""
    return ja if HAS_CJK else en


REPORT: list[str] = []


def warn(line: str = "") -> None:
    """Print a warning, tagged so the reason a figure was skipped doesn't get lost."""
    print(f"[visualize][WARN] {line}", flush=True)


def say(line: str = "") -> None:
    print(line, flush=True)
    REPORT.append(line)


def section(title: str) -> None:
    say()
    say("=" * 72)
    say(title)
    say("=" * 72)


def _despine(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


# ---------------------------------------------------------------------------
# 1. QC
# ---------------------------------------------------------------------------

def report_qc(results: Path, fig_dir: Path, sc_adata=None) -> None:
    path = results / "qc_metrics.csv"
    if not path.exists():
        say("qc_metrics.csv not found (skipping)")
        return
    qc = pd.read_csv(path, index_col=0)

    section("1. QC (qc_metrics.csv)")
    say(f"Input cells: {len(qc):,}")
    say("")
    say("Cells per flag (qc_fail = low_count | low_feature | high_pctmt):")
    for col in qc.columns:
        n = int(qc[col].sum())
        say(f"  {col:24s} {n:7,d}  ({n / len(qc):6.2%})")
    say("")
    say(f"Cells passing QC: {int((~qc['qc_fail']).sum()):,}")
    if "outlier_high_count" in qc:
        n_high = int(qc["outlier_high_count"].sum())
        say(
            f"High-count outliers: {n_high:,} cells. These are not excluded by design "
            "(malignant/aneuploid cells inherently carry more RNA; see spec sec. 5-1)."
        )
        if n_high == 0:
            say(
                "  -> 0 such cells means the high-RNA population formed its own coarse "
                "cluster and so was not flagged as an outlier within it -- evidence that "
                "stratification worked as intended."
            )
    if "outlier_high_pctmt" in qc and int(qc["outlier_high_pctmt"].sum()) == 0:
        say(
            "  [!] outlier_high_pctmt is 0. Mitochondrial gene identification may have "
            "failed. Pipeline v1.7+ identifies mito genes by GTF sequence ID, so re-run "
            "with --mito-chromosome set to the mtDNA accession "
            "(e.g. NC_002008.4 for CanFam6 / Dog10K_Boxer_Tasha)."
        )

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.1))
    counts = qc.drop(columns=["qc_fail"], errors="ignore").sum().sort_values()
    ax = axes[0]
    bars = ax.barh(counts.index, counts.values, color=C_BLUE, height=0.55)
    for rect, v in zip(bars, counts.values):
        ax.text(
            rect.get_width() + max(counts.values) * 0.02,
            rect.get_y() + rect.get_height() / 2,
            f"{int(v):,}",
            va="center",
            fontsize=8,
            color=C_INK2,
        )
    ax.set_title(T("Cells per QC flag", "Cells per QC flag"))
    ax.set_xlabel("cells")
    ax.set_xlim(0, max(max(counts.values), 1) * 1.25)
    ax.grid(axis="y", visible=False)
    _despine(ax)

    if sc_adata is not None and "total_counts" in sc_adata.obs:
        ax = axes[1]
        ax.hist(np.log10(sc_adata.obs["total_counts"] + 1), bins=60, color=C_BLUE)
        ax.set_title(T("Total UMI of passing cells (log10)", "Total UMI of passing cells (log10)"))
        ax.set_xlabel("log10(total_counts+1)")
        ax.set_ylabel("cells")
        _despine(ax)
        ax = axes[2]
        if "pct_counts_mt" in sc_adata.obs:
            ax.hist(sc_adata.obs["pct_counts_mt"], bins=60, color=C_ORANGE)
            ax.set_title(T("pctMT of passing cells", "pctMT of passing cells"))
            ax.set_xlabel("pct_counts_mt (%)")
            ax.set_ylabel("cells")
        _despine(ax)
    else:
        for ax in axes[1:]:
            ax.axis("off")
    fig.suptitle(T("Fig 1  QC breakdown", "Fig 1  QC breakdown"), y=1.04, fontsize=11, fontweight="bold")
    fig.savefig(fig_dir / "01_qc.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 2. Metacell quality
# ---------------------------------------------------------------------------

def report_metacell_quality(results: Path, fig_dir: Path, mc_obs: pd.DataFrame | None) -> None:
    path = results / "metacell_metrics.csv"
    if not path.exists():
        say("metacell_metrics.csv not found (skipping)")
        return
    m = pd.read_csv(path)

    section("2. Metacell quality (metacell_metrics.csv)")
    say(f"Metacells: {len(m):,}")
    say("")
    say("compactness = variance within the metacell in PCA space (lower = more homogeneous)")
    say("separation  = distance to the nearest metacell (higher = more distinct)")
    say("Neither has meaning in absolute terms; use them for relative comparison within this distribution.")
    say("")
    for col in ("compactness", "separation"):
        if col not in m:
            continue
        s = m[col]
        say(
            f"  {col:12s} min={s.min():.3g} q25={s.quantile(.25):.3g} "
            f"median={s.median():.3g} q75={s.quantile(.75):.3g} q95={s.quantile(.95):.3g} "
            f"max={s.max():.3g}"
        )
    if {"compactness", "separation"} <= set(m.columns) and len(m) > 2:
        rho = m[["compactness", "separation"]].corr(method="spearman").iloc[0, 1]
        say("")
        say(f"Rank correlation (Spearman) between compactness and separation: {rho:.3f}")
        say(
            "  When this correlation is high, do not cut on compactness alone. A metacell "
            "in a sparse region of the manifold is both 'loose inside' and 'far from its "
            "neighbors', so both metrics rise together -- often the result of a rare "
            "population being folded into a single metacell."
        )
        comp_cut = m["compactness"].quantile(0.95)
        sep_cut = m["separation"].median()
        bad = m[(m["compactness"] > comp_cut) & (m["separation"] < sep_cut)]
        say(
            f"  Drop candidates (compactness > q95={comp_cut:.3g} and separation < median={sep_cut:.3g}): "
            f"{len(bad)} metacells"
        )
        if len(bad):
            say("    " + ", ".join(bad["SEACell"].astype(str).tolist()[:20]))
        else:
            say("    -> No metacells need to be dropped; --drop-bad-metacells is unnecessary.")

    if "cell_type_purity" in m and m["cell_type_purity"].nunique() == 1:
        say("")
        say(
            f"[!] cell_type_purity is uniformly {m['cell_type_purity'].iloc[0]:.2f}. "
            "This does not reflect good quality -- it's trivially that value because "
            "cell_type has only one label, so the metric carries no information here."
        )

    n_panels = 2 + (1 if mc_obs is not None and "n_cells" in mc_obs else 0)
    fig, axes = plt.subplots(1, n_panels, figsize=(4.0 * n_panels, 3.4))
    axes = np.atleast_1d(axes)

    ax = axes[0]
    if {"compactness", "separation"} <= set(m.columns):
        comp_cut = m["compactness"].quantile(0.95)
        sep_cut = m["separation"].median()
        bad_mask = (m["compactness"] > comp_cut) & (m["separation"] < sep_cut)
        rare_mask = (m["compactness"] > comp_cut) & ~bad_mask
        ax.scatter(
            m.loc[~(bad_mask | rare_mask), "compactness"],
            m.loc[~(bad_mask | rare_mask), "separation"],
            s=14, color=C_GRAY, edgecolor="white", linewidth=0.5, label=T("typical", "typical"),
        )
        ax.scatter(
            m.loc[rare_mask, "compactness"], m.loc[rare_mask, "separation"],
            s=26, color=C_AQUA, edgecolor="white", linewidth=0.6,
            label=T("sparse region (keep)", "sparse region (keep)"),
        )
        ax.scatter(
            m.loc[bad_mask, "compactness"], m.loc[bad_mask, "separation"],
            s=30, color=C_RED, edgecolor="white", linewidth=0.6,
            label=T("drop candidate", "drop candidate"),
        )
        ax.axvline(comp_cut, color=C_INK2, lw=1.0, ls="--")
        ax.axhline(sep_cut, color=C_INK2, lw=1.0, ls="--")
        ax.text(comp_cut, ax.get_ylim()[1], T(" compactness q95", " compactness q95"), fontsize=7,
                color=C_INK2, va="top", ha="left")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(T("compactness (lower = more homogeneous)", "compactness (lower = more homogeneous)"))
        ax.set_ylabel(T("separation (higher = more distinct)", "separation (higher = more distinct)"))
        ax.set_title(T("Judge metacell quality in 2D", "Judge metacell quality in 2D"))
        ax.legend(fontsize=7.5, loc="lower right")
        _despine(ax)

    ax = axes[1]
    if "compactness" in m:
        ax.hist(np.log10(m["compactness"]), bins=30, color=C_BLUE)
        ax.set_xlabel("log10(compactness)")
        ax.set_ylabel(T("metacells", "metacells"))
        ax.set_title(T("compactness is heavy-tailed", "compactness is heavy-tailed"))
        _despine(ax)

    if n_panels == 3:
        ax = axes[2]
        ax.hist(mc_obs["n_cells"], bins=30, color=C_AQUA)
        ax.axvline(mc_obs["n_cells"].median(), color=C_INK, lw=1.2)
        ax.set_xlabel(T("cells per metacell", "cells per metacell"))
        ax.set_ylabel(T("metacells", "metacells"))
        ax.set_title(T("x", "x") if False else (f"Cells per metacell (median {mc_obs['n_cells'].median():.0f})" if HAS_CJK else f"Cells per metacell (median {mc_obs['n_cells'].median():.0f})"))
        _despine(ax)

    fig.suptitle(T("Fig 2  Metacell quality", "Fig 2  Metacell quality"), y=1.03, fontsize=11, fontweight="bold")
    fig.savefig(fig_dir / "02_metacell_quality.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 3-5. CNV
# ---------------------------------------------------------------------------

def report_cnv(cnv_adata, fig_dir: Path, malignant_labels: pd.Series | None = None) -> None:
    import infercnvpy as cnv

    section("3. CNV estimation (cnv_metacells.h5ad)")
    say(f"Metacells: {cnv_adata.n_obs:,} / genes: {cnv_adata.n_vars:,}")

    if "X_cnv" not in cnv_adata.obsm:
        say("[!] obsm['X_cnv'] is missing. CNV cannot be interpreted.")
        return

    x = cnv_adata.obsm["X_cnv"]
    say(f"X_cnv shape: {x.shape} (metacell x genomic window)")
    if "cnv" in cnv_adata.uns and "chr_pos" in cnv_adata.uns["cnv"]:
        chr_pos = cnv_adata.uns["cnv"]["chr_pos"]
        say(f"Chromosomes covered: {len(chr_pos)} -> {sorted(chr_pos)[:8]}{' ...' if len(chr_pos) > 8 else ''}")
        say(
            "  Any chromosome not listed here was dropped from the CNV computation. "
            "If the count is lower than expected, check the GTF chromosome naming and "
            "the exclude settings."
        )

    if "cell_type" in cnv_adata.obs:
        vc = cnv_adata.obs["cell_type"].value_counts()
        say("")
        say(f"cell_type (CNV normal reference): {vc.to_dict()}")
        if len(vc) == 1:
            say(
                "  [!] Only one cell_type. This means the run had no normal reference "
                "(baseline is the mean across all metacells), which attenuates the "
                "malignant signal. Usable for relative comparison between clones, but "
                "do not read it as 'CNV relative to normal'."
            )

    if "cnv_leiden" in cnv_adata.obs:
        say("")
        say("cnv_leiden = cluster based on CNV profile (candidate malignant subclones):")
        grp = cnv_adata.obs.groupby("cnv_leiden", observed=True)["cnv_score"]
        tbl = pd.DataFrame(
            {"n_metacells": grp.size(), "cnv_score_median": grp.median().round(4)}
        ).sort_values("cnv_score_median", ascending=False)
        for line in tbl.to_string().split("\n"):
            say("  " + line)
        say("")
        say(
            "cnv_score is the L2 norm of X_cnv (a proxy for aneuploidy burden). "
            "Its absolute value has no meaning; the rank across clusters does."
        )
        say(
            "The top cluster is the malignant-clone candidate, the bottom cluster the "
            "normal-cell candidate. The first sanity check is whether this ranking "
            "agrees with cell_type (e.g. immune cells)."
        )

    # --- Fig 3: chromosome-ordered heatmap ---
    for groupby in ("cnv_leiden", "cell_type"):
        if groupby not in cnv_adata.obs:
            continue
        if cnv_adata.obs[groupby].nunique() < 2:
            continue
        try:
            cnv.pl.chromosome_heatmap(
                cnv_adata, groupby=groupby, cmap=DIV_CMAP, show=False, save=False
            )
            fig = plt.gcf()
            fig.suptitle(
                (f"Fig 3  Chromosome CNV heatmap by {groupby}" if HAS_CJK else f"Fig 3  Chromosome CNV heatmap by {groupby}"),
                y=1.01, fontsize=11, fontweight="bold",
            )
            fig.savefig(fig_dir / f"03_cnv_heatmap_{groupby}.png", dpi=180, bbox_inches="tight")
            plt.close("all")
            say(f"-> figures/03_cnv_heatmap_{groupby}.png")
        except Exception as exc:
            say(f"[!] Failed to render chromosome_heatmap({groupby}): {exc}")
            plt.close("all")

    # --- Fig 4: cnv_score distribution and the malignant/normal cutoff ---
    if "cnv_score" in cnv_adata.obs:
        obs = cnv_adata.obs
        has_leiden = "cnv_leiden" in obs and obs["cnv_leiden"].nunique() >= 2
        fig, axes = plt.subplots(1, 2 if has_leiden else 1, figsize=(9.5 if has_leiden else 5, 3.5))
        axes = np.atleast_1d(axes)

        ax = axes[0]
        order = (
            obs.groupby("cnv_leiden", observed=True)["cnv_score"].median().sort_values().index
            if has_leiden
            else None
        )
        if has_leiden:
            data = [obs.loc[obs["cnv_leiden"] == g, "cnv_score"].values for g in order]
            ax.boxplot(
                data, labels=[str(g) for g in order], patch_artist=True,
                widths=0.55, showfliers=False,
                medianprops=dict(color=C_INK, lw=1.4),
                boxprops=dict(facecolor="#e9eef6", edgecolor=C_INK2, lw=0.8),
                whiskerprops=dict(color=C_INK2, lw=0.8),
                capprops=dict(color=C_INK2, lw=0.8),
            )
            for i, g in enumerate(order, start=1):
                y = obs.loc[obs["cnv_leiden"] == g, "cnv_score"].values
                ax.scatter(
                    np.random.default_rng(0).normal(i, 0.06, len(y)), y,
                    s=6, color=C_BLUE, alpha=0.45, linewidth=0,
                )
            ax.set_xlabel(T("cnv_leiden (by ascending median cnv_score)", "cnv_leiden (by ascending median cnv_score)"))
            ax.set_ylabel("cnv_score")
            ax.set_title(T("Aneuploidy burden per candidate clone", "Aneuploidy burden per candidate clone"))
            _despine(ax)

        ax = axes[-1]
        vals = np.sort(obs["cnv_score"].values)
        ax.plot(vals, np.linspace(0, 1, len(vals)), color=C_BLUE, lw=2)
        # Draw the boundary actually used for the malignant/normal call (not the median)
        boundary, boundary_label = None, ""
        if malignant_labels is not None:
            mal = obs.index[malignant_labels.reindex(obs.index) == "malignant"]
            if len(mal):
                boundary = float(obs.loc[mal, "cnv_score"].min())
                boundary_label = T(" min cnv_score called malignant", " min cnv_score called malignant")
        if boundary is None:
            boundary = float(np.nanmedian(obs["cnv_score"]))
            boundary_label = T(" median (reference only)", " median (reference only)")
        ax.axvline(boundary, color=C_RED, lw=1.4, ls="--")
        ax.text(boundary, 0.06, boundary_label, color=C_RED, fontsize=7.5)
        ax.set_xlabel("cnv_score")
        ax.set_ylabel(T("cumulative fraction", "cumulative fraction"))
        ax.set_title(
            T("Malignant / normal boundary", "Malignant / normal boundary")
        )
        _despine(ax)

        fig.suptitle("Fig 4  Distribution of cnv_score", y=1.03, fontsize=11, fontweight="bold")
        fig.savefig(fig_dir / "04_cnv_score.png")
        plt.close(fig)

    # --- Fig 5: UMAP of CNV space (small multiples) ---
    try:
        if "X_cnv_umap" not in cnv_adata.obsm:
            if "X_cnv_pca" not in cnv_adata.obsm:
                cnv.tl.pca(cnv_adata)
            if "cnv_neighbors" not in cnv_adata.uns:
                cnv.pp.neighbors(cnv_adata)
            cnv.tl.umap(cnv_adata)
        emb = cnv_adata.obsm["X_cnv_umap"]
        groups = (
            cnv_adata.obs["cnv_leiden"].astype(str)
            if "cnv_leiden" in cnv_adata.obs
            else pd.Series(["all"] * cnv_adata.n_obs, index=cnv_adata.obs_names)
        )
        levels = sorted(groups.unique(), key=lambda v: (len(v), v))
        n_panels = len(levels) + 2
        ncol = min(4, n_panels)
        nrow = int(np.ceil(n_panels / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.9 * nrow))
        axes = np.atleast_1d(axes).ravel()

        # (a) all points + direct labels at cluster centroids (not relying on color alone)
        ax = axes[0]
        ax.scatter(emb[:, 0], emb[:, 1], s=10, color=C_GRAY, linewidth=0)
        for g in levels:
            mask = (groups == g).values
            cx, cy = emb[mask, 0].mean(), emb[mask, 1].mean()
            ax.text(
                cx, cy, str(g), fontsize=9, fontweight="bold", color=C_INK,
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=C_INK2, lw=0.6, alpha=0.9),
            )
        ax.set_title(T("CNV-UMAP: clone layout", "CNV-UMAP: clone layout"))

        # (b) cnv_score with a single-hue sequential colormap
        ax = axes[1]
        sc_ = ax.scatter(
            emb[:, 0], emb[:, 1], s=12, c=cnv_adata.obs["cnv_score"].values,
            cmap=SEQ_CMAP, linewidth=0.3, edgecolor="white",
        )
        fig.colorbar(sc_, ax=ax, shrink=0.8, label="cnv_score")
        ax.set_title(T("Aneuploidy burden", "Aneuploidy burden"))

        # (c) small multiples per cluster
        for ax, g in zip(axes[2:], levels):
            mask = (groups == g).values
            ax.scatter(emb[~mask, 0], emb[~mask, 1], s=8, color="#e6e5e0", linewidth=0)
            ax.scatter(
                emb[mask, 0], emb[mask, 1], s=14, color=C_BLUE,
                linewidth=0.4, edgecolor="white",
            )
            ax.set_title(f"cnv_leiden = {g}  (n={int(mask.sum())})", fontsize=9)
        for ax in axes[len(levels) + 2:]:
            ax.axis("off")
        for ax in axes:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
        fig.suptitle(T("Fig 5  UMAP of CNV space", "Fig 5  UMAP of CNV space"), y=1.0, fontsize=11, fontweight="bold")
        fig.savefig(fig_dir / "05_cnv_umap.png")
        plt.close(fig)
        say("-> figures/05_cnv_umap.png")
    except Exception as exc:
        say(f"[!] Failed to render the CNV-UMAP: {exc}")
        plt.close("all")


# ---------------------------------------------------------------------------
# 5b. Mitochondrial / nuclear expression balance
# ---------------------------------------------------------------------------

def report_mito_qc(results: Path, fig_dir: Path, mc_obs: pd.DataFrame | None) -> None:
    """Report mitochondrial anomalies from metacell_mito_qc.csv (or metacell_obs.csv)."""
    df = None
    path = results / "metacell_mito_qc.csv"
    if path.exists():
        df = pd.read_csv(path, index_col=0)
    elif mc_obs is not None and "mito_anomaly" in mc_obs:
        df = mc_obs
    if df is None or "mt_frac" not in df:
        say("")
        say(
            "metacell_mito_qc.csv not found (generated by pipeline v1.8+). "
            "Skipping the mitochondrial/nuclear balance check."
        )
        return

    section("6. Mitochondrial / nuclear expression balance (metacell_mito_qc.csv)")
    say("Four axes compare each metacell to the others (all MAD-based):")
    say("  mito_high          high mtDNA fraction -> dead/apoptotic cells aggregated together")
    say("  mito_low           low mtDNA fraction -> nucleus-only capture / ambient RNA dominant")
    say("  ratio_outlier      log2 ratio of mtDNA-encoded to nuclear-encoded OXPHOS is off")
    say("  mito_heterogeneous per-cell pctMT is inconsistent within the metacell -> aggregation failure")
    say("")
    say(
        "The third one matters most. Because OXPHOS combines mtDNA-encoded and "
        "nuclear-encoded subunits in fixed stoichiometry, their ratio should stay "
        "roughly constant. A metacell where only this ratio is off can have a pctMT "
        "that looks perfectly normal, so a conventional pctMT filter can never catch "
        "it. It's a candidate for mitochondrial dysfunction or mtDNA copy-number change."
    )
    say(
        "  The metric is log2(mtDNA-encoded counts / nuclear-encoded machinery counts). "
        "Comparing two 'fractions' that share the same total-count denominator would "
        "absorb any real signal into their structural negative correlation (one rising "
        "forces the other down), so instead we look directly at the log-ratio between "
        "the two parts (a log-ratio of compositional data)."
    )
    say("")
    pct = df["mt_frac"] * 100
    say(
        f"Metacell pctMT: median {pct.median():.2f}% / "
        f"q05 {pct.quantile(.05):.2f}% / q95 {pct.quantile(.95):.2f}% / max {pct.max():.2f}%"
    )
    for col in ("mito_high", "mito_low", "ratio_outlier", "mito_heterogeneous", "mito_anomaly"):
        if col in df:
            n = int(df[col].astype(bool).sum())
            say(f"  {col:20s} {n:4d}  ({n / len(df):.1%})")

    if "mito_anomaly" in df and df["mito_anomaly"].astype(bool).any():
        cols = [
            c
            for c in ("mt_frac", "mito_nuc_log2ratio", "mt_frac_sd_within", "mito_anomaly_reason")
            if c in df
        ]
        show = df.loc[df["mito_anomaly"].astype(bool), cols].sort_values(
            "mt_frac", ascending=False
        )
        say("")
        say("Metacells flagged as anomalous:")
        for line in show.head(20).round(4).to_string().split("\n"):
            say("  " + line)

    if mc_obs is not None and "putative_malignant" in mc_obs and "mito_anomaly" in df:
        common = df.index.intersection(mc_obs.index)
        if len(common) > 4:
            tab = pd.crosstab(
                df.loc[common, "mito_anomaly"].astype(bool),
                mc_obs.loc[common, "putative_malignant"],
            )
            say("")
            say("Confounding check against the malignant label:")
            for line in tab.to_string().split("\n"):
                say("  " + line)
            if tab.shape == (2, 2):
                try:
                    from scipy.stats import fisher_exact

                    odds, pval = fisher_exact(tab.values)
                    say(f"  Fisher's exact test: odds ratio {odds:.3g}, p = {pval:.3g}")
                    if pval < 0.05:
                        say(
                            "  [!] Significant skew. The difference DE detects may reflect "
                            "cell state (dead/stressed cells) rather than tumor biology. "
                            "Re-run with --exclude-mito-anomaly and check whether the top "
                            "genes hold up."
                        )
                    else:
                        say("  The skew is not significant; confounding of DE is likely limited.")
                except Exception:
                    pass

    # --- Fig 7 ---
    has_ratio = "nuc_mito_frac" in df and df["nuc_mito_frac"].notna().any()
    has_sd = "mt_frac_sd_within" in df and df["mt_frac_sd_within"].notna().any()
    n_panels = 1 + int(has_ratio) + int(has_sd)
    fig, axes = plt.subplots(1, n_panels, figsize=(4.2 * n_panels, 3.6))
    axes = np.atleast_1d(axes)
    bad = df["mito_anomaly"].astype(bool) if "mito_anomaly" in df else pd.Series(False, index=df.index)

    ax = axes[0]
    if has_ratio:
        ax.scatter(df.loc[~bad, "nuc_mito_frac"] * 100, df.loc[~bad, "mt_frac"] * 100,
                   s=16, color=C_GRAY, edgecolor="white", linewidth=0.5,
                   label=T("typical", "typical"))
        ax.scatter(df.loc[bad, "nuc_mito_frac"] * 100, df.loc[bad, "mt_frac"] * 100,
                   s=34, color=C_RED, edgecolor="white", linewidth=0.6,
                   label=T("anomalous", "anomalous"))
        # Directly label only the anomalous points (not relying on color alone)
        for name, row in df.loc[bad].iterrows():
            if np.isfinite(row.get("nuc_mito_frac", np.nan)):
                ax.annotate(str(name).replace("SEACell-", "MC"),
                            (row["nuc_mito_frac"] * 100, row["mt_frac"] * 100),
                            fontsize=6.5, color=C_INK, xytext=(4, 3),
                            textcoords="offset points")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(T("nuclear-encoded OXPHOS fraction (%)", "nuclear-encoded OXPHOS fraction (%)"))
        ax.set_ylabel(T("mtDNA-encoded fraction (%)", "mtDNA-encoded fraction (%)"))
        ax.set_title(T("mtDNA vs nuclear stoichiometry", "mtDNA vs nuclear stoichiometry"))
        # Anomalies tend to fall lower-right / upper-left, so place the legend lower-left
        # (avoids overlapping the direct labels)
        ax.legend(fontsize=7.5, loc="lower left")
    else:
        ax.hist(df["mt_frac"] * 100, bins=40, color=C_BLUE)
        ax.set_xlabel("pctMT (%)"); ax.set_ylabel(T("metacells", "metacells"))
        ax.set_title(T("pctMT per metacell", "pctMT per metacell"))
    _despine(ax)

    i = 1
    if has_ratio:
        ax = axes[i]; i += 1
        r = df["mito_nuc_log2ratio"].dropna()
        ax.hist(r, bins=40, color=C_BLUE)
        med = np.median(r); mad = np.median(np.abs(r - med)) * 1.4826
        for k in (-3.5, 3.5):
            ax.axvline(med + k * mad, color=C_RED, lw=1.2, ls="--")
        ax.set_xlabel(T("log2(mtDNA-encoded / nuclear-encoded)", "log2(mtDNA-encoded / nuclear-encoded)"))
        ax.set_ylabel(T("metacells", "metacells"))
        ax.set_title(T("Stoichiometry log-ratio (±3.5 MAD)", "Stoichiometry log-ratio (±3.5 MAD)"))
        _despine(ax)

    if has_sd:
        ax = axes[i]
        ax.scatter(df.loc[~bad, "mt_frac"] * 100, df.loc[~bad, "mt_frac_sd_within"],
                   s=16, color=C_GRAY, edgecolor="white", linewidth=0.5)
        ax.scatter(df.loc[bad, "mt_frac"] * 100, df.loc[bad, "mt_frac_sd_within"],
                   s=34, color=C_RED, edgecolor="white", linewidth=0.6)
        ax.set_xlabel(T("pctMT per metacell (%)", "pctMT per metacell (%)"))
        ax.set_ylabel(T("SD of pctMT among member cells", "SD of pctMT among member cells"))
        ax.set_title(T("Uniform mean, mixed contents", "Uniform mean, mixed contents"))
        _despine(ax)

    fig.suptitle(T("Fig 7  Mito / nuclear balance", "Fig 7  Mito / nuclear balance"),
                 y=1.03, fontsize=11, fontweight="bold")
    fig.savefig(fig_dir / "07_mito_balance.png")
    plt.close(fig)
    say("")
    say("-> figures/07_mito_balance.png")


# ---------------------------------------------------------------------------
# 6. DE
# ---------------------------------------------------------------------------

def report_de(
    results: Path, fig_dir: Path, mc_obs: pd.DataFrame | None, top_n: int = 25
) -> None:
    path = results / "de_malignant_vs_normal.csv"
    if not path.exists():
        say("de_malignant_vs_normal.csv not found (skipping)")
        return
    de = pd.read_csv(path, index_col=0)

    section("4. DE analysis (de_malignant_vs_normal.csv)")
    say("Column meanings:")
    say("  baseMean        mean normalized count across all metacells (expression level)")
    say("  log2FoldChange  log2 ratio of malignant / normal; positive = higher in malignant")
    say("  lfcSE           standard error of log2FoldChange; larger = less stable estimate")
    say("  stat            Wald statistic (log2FoldChange / lfcSE)")
    say("  pvalue / padj   Wald test p-value and BH-adjusted q-value; use padj")
    say("")

    n_total = len(de)
    n_tested = int(de["padj"].notna().sum())
    say(f"Genes tested: {n_total:,} (with a padj value: {n_tested:,})")
    say(
        "  Genes with padj = NaN were excluded by pyDESeq2's independent filtering or "
        "Cook's distance and were not tested. Do not include them in the denominator."
    )
    for thr in (0.05, 0.01):
        n_sig = int((de["padj"] < thr).sum())
        say(f"  padj < {thr}: {n_sig:,}  ({n_sig / max(n_tested, 1):.1%} of tested)")
    strong = de[(de["padj"] < 0.05) & (de["log2FoldChange"].abs() > 1)]
    say(f"  padj < 0.05 and |log2FC| > 1: {len(strong):,}")
    say(f"    higher in malignant: {int((strong['log2FoldChange'] > 0).sum()):,}")
    say(f"    lower in malignant: {int((strong['log2FoldChange'] < 0).sum()):,}")

    if n_tested and (de["padj"] < 0.05).sum() / n_tested > 0.5:
        say("")
        say(
            "[!] More than half of the tested genes have padj < 0.05. This is not "
            "necessarily a biological finding -- it can indicate that the group "
            "definition itself derives from the expression profile (circular "
            "reasoning). Be sure to read the 'How the malignant label is built' note below."
        )

    sig = de[de["padj"] < 0.05].copy()
    if len(sig):
        sig["rank_score"] = -np.log10(sig["pvalue"].clip(lower=1e-300)) * sig["log2FoldChange"].abs()
        top = pd.concat(
            [
                sig.sort_values("log2FoldChange", ascending=False).head(top_n),
                sig.sort_values("log2FoldChange").head(top_n),
            ]
        )
        top.to_csv(results / "de_top_genes.csv")
        say("")
        say("Top 10 genes higher in malignant (padj < 0.05, descending log2FC):")
        for line in (
            sig.sort_values("log2FoldChange", ascending=False)
            .head(10)[["baseMean", "log2FoldChange", "padj"]]
            .round({"baseMean": 1, "log2FoldChange": 2})
            .to_string()
            .split("\n")
        ):
            say("  " + line)
        say("")
        say("-> saved the top/bottom %d genes each to de_top_genes.csv" % top_n)
        n_loc = int(sig.index.astype(str).str.startswith("LOC").sum())
        if n_loc:
            say(
                f"  Note: of the {len(sig):,} significant genes, {n_loc:,} "
                f"({n_loc / len(sig):.0%}) are unnamed LOC* IDs, which are hard to "
                "interpret functionally (see spec sec. 2.2.2). They need mapping to "
                "human orthologs."
            )

    # --- Fig 6: volcano + MA ---
    d = de.dropna(subset=["padj", "log2FoldChange"]).copy()
    if not len(d):
        return
    d["nlp"] = -np.log10(d["pvalue"].clip(lower=1e-300))
    up = (d["padj"] < 0.05) & (d["log2FoldChange"] > 1)
    down = (d["padj"] < 0.05) & (d["log2FoldChange"] < -1)
    ns = ~(up | down)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    ax = axes[0]
    ax.scatter(d.loc[ns, "log2FoldChange"], d.loc[ns, "nlp"], s=5, color=C_GRAY, linewidth=0,
               label=f"n.s. ({int(ns.sum()):,})")
    ax.scatter(d.loc[up, "log2FoldChange"], d.loc[up, "nlp"], s=9, color=C_RED, linewidth=0,
               label=(f"up in malignant ({int(up.sum()):,})" if HAS_CJK else f"up in malignant ({int(up.sum()):,})"))
    ax.scatter(d.loc[down, "log2FoldChange"], d.loc[down, "nlp"], s=9, color=C_BLUE, linewidth=0,
               label=(f"down in malignant ({int(down.sum()):,})" if HAS_CJK else f"down in malignant ({int(down.sum()):,})"))
    ax.axhline(-np.log10(0.05), color=C_INK2, lw=0.8, ls=":")
    for x0 in (-1, 1):
        ax.axvline(x0, color=C_INK2, lw=0.8, ls=":")
    # Directly label only the top genes (not every point)
    labelled = pd.concat(
        [
            d[up].sort_values("nlp", ascending=False).head(6),
            d[down].sort_values("nlp", ascending=False).head(6),
        ]
    )
    for gene, row in labelled.iterrows():
        ax.annotate(
            str(gene), (row["log2FoldChange"], row["nlp"]),
            fontsize=7, color=C_INK, xytext=(3, 3), textcoords="offset points",
        )
    ax.set_xlabel(T("log2 fold change (malignant / normal)", "log2 fold change (malignant / normal)"))
    ax.set_ylabel("-log10(p value)")
    ax.set_title("Volcano")
    ax.legend(fontsize=7.5, loc="upper left")
    _despine(ax)

    ax = axes[1]
    ax.scatter(d.loc[ns, "baseMean"], d.loc[ns, "log2FoldChange"], s=5, color=C_GRAY, linewidth=0)
    ax.scatter(d.loc[up, "baseMean"], d.loc[up, "log2FoldChange"], s=9, color=C_RED, linewidth=0)
    ax.scatter(d.loc[down, "baseMean"], d.loc[down, "log2FoldChange"], s=9, color=C_BLUE, linewidth=0)
    ax.axhline(0, color=C_INK2, lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel(T("baseMean (expression level)", "baseMean (expression level)"))
    ax.set_ylabel("log2 fold change")
    ax.set_title(T("MA: wider log2FC at low expression is expected", "MA: wider log2FC at low expression is expected"))
    _despine(ax)

    fig.suptitle(T("Fig 6  DE: malignant vs normal metacells", "Fig 6  DE: malignant vs normal metacells"), y=1.02, fontsize=11, fontweight="bold")
    fig.savefig(fig_dir / "06_de_volcano.png")
    plt.close(fig)

    # --- Interpretation notes ---
    section("5. Interpretation notes (read these)")
    if mc_obs is not None and "putative_malignant" in mc_obs:
        vc = mc_obs["putative_malignant"].value_counts().to_dict()
        say(f"putative_malignant breakdown: {vc}")
    say(
        "1. Always check how the malignant label was built. The rationale is recorded "
        "in malignant_call.txt."
    )
    say(
        "   Pipeline v1.4 and earlier split on the median cnv_score. By construction "
        "that splits roughly 50:50 and is not a tumor-cell detection -- this DE is then "
        "just 'high aneuploidy-burden metacells vs low aneuploidy-burden metacells'. If "
        "malignant is around 50% and malignant_call.txt is absent (or says the median "
        "method was used), do not read the result as biological malignant vs normal. "
        "Re-run with v1.5 or later."
    )
    say(
        "   v1.5+ calls malignancy per clone (cnv_leiden). With a normal reference the "
        "threshold is 'normal-reference cnv_score mean + 3SD'; without one it splits at "
        "the largest gap between clone medians, and abstains from calling when "
        "bimodality is not clear. Check figure 4 for whether the boundary falls in a "
        "valley between clusters, and figure 3's heatmap for whether the two groups' "
        "CNV patterns genuinely differ."
    )
    say(
        "2. With a single specimen, sample_id drops out of the design. The resulting "
        "p-values reflect 'variation among metacells within this tumor' and do not "
        "generalize across patients (the pseudoreplication issue remains in principle; "
        "see spec sec. 5-4)."
    )
    say(
        "3. infercnvpy is explicitly marked experimental by its developers. Cross-check "
        "any important conclusion with an independent method such as CopyKAT or SCEVAN "
        "(see spec sec. 5-6)."
    )
    say(
        "4. Without a secured normal reference (cell_type), CNV is relative to the mean "
        "across all metacells. In a specimen where tumor cells dominate, that mean "
        "itself carries aneuploidy, so the CNV amplitude is underestimated."
    )


# ---------------------------------------------------------------------------
# 7. Chromosome-level CNV heatmap (metacells ordered by CNV clustering)
# ---------------------------------------------------------------------------

_NC_ACC = re.compile(r"^(?:chr)?(N[CTWZ]_\d+)(?:\.\d+)?$")

#: Known non-chromosome labels to exclude when laying out the figure's x-axis
_NON_CHROM = {"chrM", "chrMT", "chrMito"}


def load_chromosome_map(path: Path | None) -> dict[str, str]:
    """Build an accession -> chromosome name mapping from an NCBI assembly report
    (or a 2-column TSV).

    Assembly report columns are
        Sequence-Name  Sequence-Role  Assigned-Molecule ... RefSeq-Accession ...
    A 2-column TSV is treated as "accession<TAB>chromosome name".
    """
    if path is None or not Path(path).exists():
        return {}
    mapping: dict[str, str] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) >= 7:
                # assembly report columns: 0=Sequence-Name, 2=Assigned-Molecule, 6=RefSeq-Accn
                mol, acc = f[2].strip(), f[6].strip()
                if acc and acc not in ("na", "-"):
                    if mol and mol not in ("na", "-"):
                        mapping[acc] = mol
                    else:
                        mapping[acc] = f[0].strip()
            elif len(f) >= 2:
                mapping[f[0].strip()] = f[1].strip()
    return mapping


def chromosome_display_names(
    chrom_names: list[str], chromosome_map: dict[str, str] | None = None
) -> tuple[dict[str, str], bool]:
    """Replace accessions with short chromosome names for the heatmap's axis labels.

    Returns (label_map, inferred). inferred=True means the numbers were
    "inferred from accession order" rather than measured.
    """
    chromosome_map = chromosome_map or {}
    labels: dict[str, str] = {}
    unresolved: list[str] = []
    for name in chrom_names:
        m = _NC_ACC.match(str(name))
        key_full = str(name)[3:] if str(name).startswith("chr") else str(name)
        hit = chromosome_map.get(key_full) or (chromosome_map.get(m.group(1)) if m else None)
        if hit:
            labels[name] = hit if str(hit).lower().startswith(("chr",)) else str(hit)
        elif m:
            unresolved.append(name)
        else:
            labels[name] = str(name).replace("chr", "", 1) or str(name)

    # Unplaced scaffolds (NW_/NT_) have no chromosome number, so exclude them from inference.
    chrom_like = [n for n in unresolved if _NC_ACC.match(str(n)).group(1).startswith("NC_")]
    for name in unresolved:
        if name not in chrom_like:
            labels[name] = _NC_ACC.match(str(name)).group(1)

    inferred = False
    if chrom_like:
        # Numbers can only be assigned under the assumption that accession order =
        # chromosome number order. Mark them as "inferred" in the figure (a * suffix).
        inferred = True
        ordered = sorted(chrom_like, key=lambda s: _NC_ACC.match(str(s)).group(1))
        n = len(ordered)
        for i, name in enumerate(ordered):
            # Treat the last one as X by convention for autosomes + X assemblies.
            # If the count is small, a subset was likely passed in, so keep it numeric.
            labels[name] = "X*" if (i == n - 1 and n >= 10) else f"{i + 1}*"
    return labels, inferred


def chromosome_cnv_matrix(cnv_adata, exclude_unplaced: bool = True):
    """Aggregate obsm['X_cnv'] into per-chromosome values using uns['cnv']['chr_pos'] boundaries.

    Returns (mean_df, gain_df, loss_df), all metacell x chromosome.
      mean_df : mean CNV of the windows within the chromosome (amplitude)
      gain_df : fraction of windows above +threshold (breadth of gain)
      loss_df : fraction of windows below -threshold (breadth of loss)
    """
    import scipy.sparse as sp

    if "cnv" not in cnv_adata.uns or "chr_pos" not in cnv_adata.uns["cnv"]:
        raise ValueError("uns['cnv']['chr_pos'] is missing")
    x = cnv_adata.obsm["X_cnv"]
    x = np.asarray(x.todense()) if sp.issparse(x) else np.asarray(x)

    chr_pos = dict(cnv_adata.uns["cnv"]["chr_pos"])
    items = sorted(chr_pos.items(), key=lambda kv: int(kv[1]))
    bounds = []
    for i, (name, start) in enumerate(items):
        end = int(items[i + 1][1]) if i + 1 < len(items) else x.shape[1]
        bounds.append((str(name), int(start), int(end)))

    keep = []
    for name, s, e in bounds:
        if e - s < 1:
            continue
        if name in _NON_CHROM:
            continue
        if exclude_unplaced and _NC_ACC.match(name) and not _NC_ACC.match(name).group(1).startswith("NC_"):
            continue
        keep.append((name, s, e))

    # The threshold is for judging the breadth of "nonzero" windows. X_cnv is a
    # smoothed relative value, not an absolute copy number, so the threshold is
    # derived from the overall distribution via MAD.
    nz = x[x != 0]
    if nz.size:
        thr = float(1.4826 * np.median(np.abs(nz - np.median(nz))))
        thr = max(thr, 1e-3)
    else:
        thr = 1e-3

    names = [n for n, _, _ in keep]
    mean = np.zeros((x.shape[0], len(keep)))
    gain = np.zeros_like(mean)
    loss = np.zeros_like(mean)
    for j, (_, s, e) in enumerate(keep):
        blk = x[:, s:e]
        mean[:, j] = blk.mean(axis=1)
        gain[:, j] = (blk > thr).mean(axis=1)
        loss[:, j] = (blk < -thr).mean(axis=1)
    idx = pd.Index(cnv_adata.obs_names, name="SEACell")
    return (
        pd.DataFrame(mean, index=idx, columns=names),
        pd.DataFrame(gain, index=idx, columns=names),
        pd.DataFrame(loss, index=idx, columns=names),
        thr,
    )


def cluster_order(mat: np.ndarray, n_clusters: int | None = None):
    """Hierarchically cluster CNV profiles and return the row order, linkage, and cluster labels.

    - Distance is Euclidean, linkage is Ward. Correlation distance is avoided because
      it becomes NaN for a zero-variance metacell.
    - optimal_leaf_ordering reorders leaves so neighbors are as similar as possible.
    """
    from scipy.cluster.hierarchy import dendrogram, fcluster, linkage, optimal_leaf_ordering
    from scipy.spatial.distance import pdist

    if mat.shape[0] < 3:
        return np.arange(mat.shape[0]), None, np.zeros(mat.shape[0], dtype=int)
    d = pdist(mat, metric="euclidean")
    z = linkage(d, method="ward")
    try:
        z = optimal_leaf_ordering(z, d)
    except Exception:
        pass
    order = dendrogram(z, no_plot=True)["leaves"]
    labels = np.zeros(mat.shape[0], dtype=int)
    if n_clusters and n_clusters >= 2:
        labels = fcluster(z, t=int(n_clusters), criterion="maxclust")
    return np.asarray(order), z, labels


def _annotation_strip(ax, values: pd.Series, palette: dict, title: str) -> None:
    """Draw a single-column categorical annotation strip (pass rows in the same order as the heatmap)."""
    rgba = np.zeros((len(values), 1, 4))
    for i, v in enumerate(values.astype(str).values):
        rgba[i, 0] = matplotlib.colors.to_rgba(palette.get(v, C_GRAY))
    ax.imshow(rgba, aspect="auto", interpolation="nearest", origin="upper")
    ax.set_xticks([0])
    ax.set_xticklabels([title], rotation=90, fontsize=7)
    ax.set_yticks([])
    ax.grid(False)
    for s in ax.spines.values():
        s.set_linewidth(0.4)


_STRIP_COLORS = [C_BLUE, C_ORANGE, C_AQUA, C_RED, "#8e6bbf", "#c9a227", "#4c8f8b", "#b0567a",
                 "#5b7fa6", "#d98c5f", "#7aa64c", "#a3527a", "#6f6f6f", "#2f6f4f"]


def report_chromosome_cnv_heatmap(
    cnv_adata,
    fig_dir: Path,
    results: Path,
    malignant_labels: pd.Series | None = None,
    chromosome_map_path: Path | None = None,
    n_clusters: int | None = None,
    exclude_unplaced: bool = True,
) -> None:
    """Fig 8: chromosome x metacell gain/loss heatmap. Metacells are ordered by CNV clustering."""
    section("6. Chromosome-level CNV heatmap (fig 8)")

    if "X_cnv" not in cnv_adata.obsm:
        say("[!] Cannot build this figure: obsm['X_cnv'] is missing.")
        return
    try:
        mean_df, gain_df, loss_df, thr = chromosome_cnv_matrix(
            cnv_adata, exclude_unplaced=exclude_unplaced
        )
    except Exception as exc:
        say(f"[!] Failed to aggregate to chromosome level: {exc}")
        return
    if mean_df.shape[1] < 2:
        say("[!] One or fewer chromosomes were aggregated. Check the GTF chromosome naming.")
        return

    chrom_map = load_chromosome_map(chromosome_map_path)
    labels, inferred = chromosome_display_names(list(mean_df.columns), chrom_map)

    # --- Column (chromosome) order: natural order of display names ---
    def _chrkey(name: str):
        lab = labels[name].rstrip("*")
        try:
            return (0, float(lab.replace("chr", "")))
        except ValueError:
            return (1, lab)

    col_order = sorted(mean_df.columns, key=_chrkey)
    mean_df = mean_df[col_order]
    gain_df = gain_df[col_order]
    loss_df = loss_df[col_order]

    obs = cnv_adata.obs
    if n_clusters is None:
        n_clusters = int(obs["cnv_leiden"].nunique()) if "cnv_leiden" in obs else 4
    order, z, cut = cluster_order(mean_df.to_numpy(), n_clusters=n_clusters)

    mean_ord = mean_df.iloc[order]
    row_names = mean_ord.index
    cut_ord = pd.Series(cut, index=mean_df.index).iloc[order]

    say(f"Chromosomes aggregated: {mean_df.shape[1]} / metacells: {mean_df.shape[0]}")
    say(f"Gain/loss call threshold (per window): |X_cnv| > {thr:.4f} (auto-set from the MAD of nonzero windows)")
    if inferred:
        say(
            "[!] Chromosome numbers marked with * in the figure were inferred from "
            "accession order. Pass --chromosome-map with an NCBI assembly report to "
            "get the actual numbers."
        )
    if exclude_unplaced:
        say("Unplaced scaffolds such as NW_* and chrM are excluded from the figure.")
    say(
        "Row order: not by ID, but by hierarchical clustering of per-chromosome mean "
        "CNV (Euclidean + Ward, optimal leaf ordering), placing similar metacells next "
        "to each other."
    )

    # --- Prepare annotation columns ---
    strips: list[tuple[str, pd.Series, dict]] = []
    cut_lab = cut_ord.astype(str)
    strips.append(
        (
            T("CNV cluster", "CNV cluster"),
            cut_lab,
            {v: _STRIP_COLORS[i % len(_STRIP_COLORS)] for i, v in enumerate(sorted(cut_lab.unique(), key=lambda s: (len(s), s)))},
        )
    )
    if "cnv_leiden" in obs:
        v = obs["cnv_leiden"].astype(str).reindex(row_names)
        strips.append(("cnv_leiden", v, {g: _STRIP_COLORS[i % len(_STRIP_COLORS)] for i, g in enumerate(sorted(v.unique(), key=lambda s: (len(s), s)))}))
    if "cell_type" in obs:
        v = obs["cell_type"].astype(str).reindex(row_names)
        strips.append(("cell_type", v, {g: _STRIP_COLORS[i % len(_STRIP_COLORS)] for i, g in enumerate(sorted(v.unique()))}))
    if malignant_labels is not None:
        v = malignant_labels.reindex(row_names).astype(str)
        strips.append((T("malignant call", "malignant call"), v, {"malignant": C_RED, "normal": C_BLUE, "nan": C_GRAY, "undetermined": C_GRAY}))
    if "mito_anomaly" in obs:
        v = obs["mito_anomaly"].reindex(row_names).astype(str)
        strips.append((T("mt anomaly", "mt anomaly"), v, {"True": C_ORANGE, "False": "#eeeeea"}))

    # --- Layout ---
    n_row, n_col = mean_ord.shape
    h = max(4.0, min(14.0, 0.028 * n_row + 2.4))
    fig = plt.figure(figsize=(2.0 + 0.20 * n_col + 0.22 * len(strips), h + 1.9))
    gs = fig.add_gridspec(
        2, 2 + len(strips),
        width_ratios=[0.9] + [0.20 * n_col] + [0.22] * len(strips),
        height_ratios=[h, 1.5],
        wspace=0.06, hspace=0.14,
        left=0.055, right=0.895, top=0.955, bottom=0.105,
    )

    # (a) dendrogram
    ax_d = fig.add_subplot(gs[0, 0])
    if z is not None:
        from scipy.cluster.hierarchy import dendrogram
        with plt.rc_context({"lines.linewidth": 0.6}):
            dendrogram(
                z, orientation="left", no_labels=True, color_threshold=0,
                above_threshold_color=C_INK2, ax=ax_d,
            )
        ax_d.invert_yaxis()
    ax_d.set_xticks([])
    ax_d.set_yticks([])
    ax_d.grid(False)
    for s in ax_d.spines.values():
        s.set_visible(False)
    ax_d.set_title(T("CNV\nclustering", "CNV\nclustering"), fontsize=8)
    ax_d.set_ylabel(
        T(f"metacells (n={n_row})", f"metacells (n={n_row})"), fontsize=8, labelpad=2
    )

    # (b) main heatmap
    ax = fig.add_subplot(gs[0, 1])
    vmax = float(np.nanpercentile(np.abs(mean_ord.to_numpy()), 99))
    vmax = max(vmax, 1e-4)
    im = ax.imshow(
        mean_ord.to_numpy(), aspect="auto", cmap=DIV_CMAP,
        vmin=-vmax, vmax=vmax, interpolation="nearest", origin="upper",
    )
    ax.set_xticks(range(n_col))
    ax.set_xticklabels([labels[c] for c in mean_ord.columns], fontsize=7, rotation=90)
    ax.set_yticks([])
    ax.grid(False)
    for x0 in np.arange(0.5, n_col - 0.5):
        ax.axvline(x0, color="white", lw=0.35)
    # Horizontal line at each cluster boundary
    changes = np.flatnonzero(cut_ord.values[1:] != cut_ord.values[:-1])
    for y0 in changes:
        ax.axhline(y0 + 0.5, color=C_INK, lw=0.7)
    ax.set_title(
        T("Mean CNV per chromosome (red=gain / blue=loss)", "Mean CNV per chromosome (red=gain / blue=loss)"),
        fontsize=9,
    )

    # (c) annotation strips
    for k, (title, values, palette) in enumerate(strips):
        axs = fig.add_subplot(gs[0, 2 + k])
        _annotation_strip(axs, values, palette, title)

    # (d) bottom row: breadth of gain/loss per chromosome
    ax2 = fig.add_subplot(gs[1, 1], sharex=ax)
    if malignant_labels is not None and (malignant_labels.reindex(mean_df.index) == "malignant").any():
        sel = (malignant_labels.reindex(mean_df.index) == "malignant").values
        sel_label = T("metacells called malignant", "metacells called malignant")
    else:
        sel = np.ones(mean_df.shape[0], dtype=bool)
        sel_label = T("all metacells", "all metacells")
    g = gain_df[mean_ord.columns].to_numpy()[sel].mean(axis=0)
    l = loss_df[mean_ord.columns].to_numpy()[sel].mean(axis=0)
    xpos = np.arange(n_col)
    ax2.bar(xpos, g, color=C_RED, width=0.7, label=T("fraction of windows gained", "fraction of windows gained"))
    ax2.bar(xpos, -l, color=C_BLUE, width=0.7, label=T("fraction of windows lost", "fraction of windows lost"))
    ax2.axhline(0, color=C_INK2, lw=0.8)
    ax2.set_xlim(-0.5, n_col - 0.5)
    mx = float(max(np.nanmax(g) if g.size else 0.0, np.nanmax(l) if l.size else 0.0, 1e-3))
    ax2.set_ylim(-mx * 1.30, mx * 1.30)
    ax2.set_ylabel(T("window fraction", "window fraction"), fontsize=8)
    ax2.set_xlabel(T("chromosome", "chromosome"), fontsize=8)
    ax2.legend(fontsize=7, ncol=2, loc="upper right")
    ax2.set_title(sel_label + T(": breadth of gains / losses", ": breadth of gains / losses"), fontsize=8.5)
    ax2.tick_params(labelbottom=True)
    _despine(ax2)

    # (e) colorbar and legend
    cax = fig.add_axes([0.925, 0.60, 0.011, 0.24])
    fig.colorbar(im, cax=cax, label=T("mean CNV", "mean CNV"))
    handles = []
    for title, values, palette in strips:
        seen = [v for v in palette if v in set(values.astype(str))]
        for v in seen:
            handles.append(
                matplotlib.patches.Patch(facecolor=palette[v], edgecolor="none", label=f"{title}: {v}")
            )
    if handles:
        fig.legend(
            handles=handles, loc="lower center", ncol=min(6, max(2, (len(handles) + 3) // 4)),
            fontsize=6.8, bbox_to_anchor=(0.5, 0.002), frameon=False,
        )

    fig.suptitle(
        T("Fig 8  Chromosome-level CNV (metacells ordered by CNV clustering)",
          "Fig 8  Chromosome-level CNV (metacells ordered by CNV clustering)"),
        y=0.992, fontsize=11, fontweight="bold",
    )
    out = fig_dir / "08_cnv_chromosome_heatmap.png"
    fig.savefig(out, dpi=190)
    plt.close(fig)
    say(f"-> figures/{out.name}")

    # --- Also keep the numbers ---
    tidy = mean_df.copy()
    tidy.columns = [labels[c] for c in tidy.columns]
    tidy.insert(0, "cnv_cluster_from_heatmap", pd.Series(cut, index=mean_df.index))
    tidy.insert(1, "heatmap_row_order", pd.Series(np.argsort(order), index=mean_df.index[order]).reindex(mean_df.index))
    for c in ("cnv_leiden", "cell_type", "cnv_score"):
        if c in obs:
            tidy.insert(2, c, obs[c].reindex(mean_df.index).astype(str) if c != "cnv_score" else obs[c].reindex(mean_df.index))
    if malignant_labels is not None:
        tidy.insert(2, "putative_malignant", malignant_labels.reindex(mean_df.index))
    tidy.to_csv(results / "cnv_chromosome_means.csv")
    say("-> cnv_chromosome_means.csv (mean CNV per chromosome and cluster number)")

    # --- Quantify cross-cluster agreement (whether the tumor looks single-origin) ---
    if malignant_labels is not None:
        mal = (malignant_labels.reindex(mean_df.index) == "malignant").values
        if mal.sum() >= 2:
            prof = mean_df.to_numpy()[mal].mean(axis=0)
            aff = np.flatnonzero(np.abs(prof) > np.percentile(np.abs(prof), 60))
            if len(aff):
                sub = mean_df.to_numpy()[mal][:, aff]
                agree = float(np.mean(np.sign(sub) == np.sign(prof[aff])[None, :]))
                say("")
                say(
                    f"Sign agreement among malignant metacells (top {len(aff)} chromosomes by "
                    f"magnitude of change): {agree:.1%}"
                )
                say(
                    "  If the tumor is single-origin, the chromosomes that gain/lose should "
                    "be shared across metacells. A low value here (roughly < 0.8) suggests "
                    "multiple lineages or contamination by false positives."
                )


# ---------------------------------------------------------------------------
# 8. Top-N DEG swarm plot, by cell type group
# ---------------------------------------------------------------------------

def _lineage_sort_key(label: str):
    """'Myeloid_2' -> ('Myeloid', 2). Sort key that groups by lineage name and orders by branch number."""
    pref = ["Other", "Myeloid", "T/NK", "B", "Leukocyte"]
    s = str(label)
    base, num = s, -1
    if "_" in s:
        head, _, tail = s.rpartition("_")
        if tail.isdigit():
            base, num = head, int(tail)
    rank = pref.index(base) if base in pref else len(pref)
    return (rank, base, num)


def _beeswarm_offsets(y: np.ndarray, width: float = 0.34, n_bins: int | None = None) -> np.ndarray:
    """Bin y and spread points within the same bin left/right (a simple beeswarm, no seaborn dependency)."""
    y = np.asarray(y, dtype=float)
    if y.size == 0:
        return np.zeros(0)
    if n_bins is None:
        # Coarse bins look like horizontal stripes, so scale bin count with point count
        n_bins = int(np.clip(y.size / 3.5, 18, 90))
    lo, hi = np.nanmin(y), np.nanmax(y)
    if not np.isfinite(lo) or hi <= lo:
        bins = np.zeros(y.size, dtype=int)
    else:
        bins = np.clip(((y - lo) / (hi - lo) * n_bins).astype(int), 0, n_bins - 1)
    off = np.zeros(y.size)
    for b in np.unique(bins):
        idx = np.flatnonzero(bins == b)
        k = idx.size
        if k == 1:
            continue
        # Place points alternately outward from the center
        pos = np.arange(k) - (k - 1) / 2.0
        step = min(width / max((k - 1) / 2.0, 1.0), 0.075)
        off[idx[np.argsort(y[idx])]] = pos * step
    return np.clip(off, -width, width)


def _metacell_expression(results: Path, cnv_adata):
    """Return the metacell x gene expression matrix as (DataFrame, description).

    Priority order:
      1. Raw counts from metacells.h5ad -> CPM -> log1p
      2. layers['counts'] from cnv_metacells.h5ad -> CPM -> log1p
      3. X of cnv_metacells.h5ad (already normalized)
    """
    import scipy.sparse as sp

    def _cpm_log1p(mat):
        mat = np.asarray(mat.todense()) if sp.issparse(mat) else np.asarray(mat, dtype=float)
        tot = mat.sum(axis=1, keepdims=True)
        tot[tot == 0] = 1.0
        return np.log1p(mat / tot * 1e6)

    p = results / "metacells.h5ad"
    if p.exists():
        try:
            import anndata as ad
            a = ad.read_h5ad(p)
            src = a.layers["counts"] if "counts" in a.layers else a.X
            return (
                pd.DataFrame(_cpm_log1p(src), index=a.obs_names, columns=a.var_names),
                T("log1p(CPM) from raw counts in metacells.h5ad",
                  "log1p(CPM) from raw counts in metacells.h5ad"),
            )
        except Exception:
            pass
    if cnv_adata is None:
        return None, ""
    if "counts" in getattr(cnv_adata, "layers", {}):
        return (
            pd.DataFrame(
                _cpm_log1p(cnv_adata.layers["counts"]),
                index=cnv_adata.obs_names, columns=cnv_adata.var_names,
            ),
            T("log1p(CPM) from layers['counts'] in cnv_metacells.h5ad",
              "log1p(CPM) from layers['counts'] in cnv_metacells.h5ad"),
        )
    x = cnv_adata.X
    x = np.asarray(x.todense()) if sp.issparse(x) else np.asarray(x)
    return (
        pd.DataFrame(x, index=cnv_adata.obs_names, columns=cnv_adata.var_names),
        T("X of cnv_metacells.h5ad (already normalised by the pipeline)",
          "X of cnv_metacells.h5ad (already normalised by the pipeline)"),
    )


def select_top_degs(
    de: pd.DataFrame, top_n: int = 20, rank: str = "balanced",
    padj_max: float = 0.05, min_abs_lfc: float = 0.5, min_base_mean: float = 5.0,
) -> pd.DataFrame:
    """Select the genes shown in the swarm plot.

    rank:
      balanced  half up-regulated, half down-regulated (default). Because many genes
                tie at padj=0, ranking uses -log10(p) * |log2FC|.
      lfc       descending |log2FC|
      padj      ascending padj (ties broken by |log2FC|)
    """
    d = de.dropna(subset=["padj", "log2FoldChange"]).copy()
    d = d[(d["padj"] < padj_max) & (d["log2FoldChange"].abs() >= min_abs_lfc)]
    if "baseMean" in d:
        d = d[d["baseMean"] >= min_base_mean]
    if not len(d):
        return d
    d["nlp"] = -np.log10(d["pvalue"].clip(lower=1e-300))
    d["rank_score"] = d["nlp"] * d["log2FoldChange"].abs()
    if rank == "lfc":
        return d.reindex(d["log2FoldChange"].abs().sort_values(ascending=False).index).head(top_n)
    if rank == "padj":
        return d.sort_values(["padj", "rank_score"], ascending=[True, False]).head(top_n)
    up = d[d["log2FoldChange"] > 0].sort_values("rank_score", ascending=False)
    dn = d[d["log2FoldChange"] < 0].sort_values("rank_score", ascending=False)
    k = top_n // 2
    picked = pd.concat([up.head(k), dn.head(top_n - k)])
    if len(picked) < top_n:  # if one side runs short, fill in from the other
        rest = d.drop(index=picked.index).sort_values("rank_score", ascending=False)
        picked = pd.concat([picked, rest.head(top_n - len(picked))])
    return picked.sort_values("log2FoldChange", ascending=False)


def report_deg_swarm(
    results: Path,
    fig_dir: Path,
    cnv_adata,
    mc_obs: pd.DataFrame | None,
    top_n: int = 20,
    rank: str = "balanced",
    group_key: str = "cell_type",
) -> None:
    """Fig 9: draw a swarm plot of top-N DEG metacell expression, by cell type group."""
    section("7. Top DEG swarm plot by group (fig 9)")

    path = results / "de_malignant_vs_normal.csv"
    if not path.exists():
        say("de_malignant_vs_normal.csv not found; skipping.")
        return
    de = pd.read_csv(path, index_col=0)
    top = select_top_degs(de, top_n=top_n, rank=rank)
    if not len(top):
        say("[!] No significant genes meet the criteria (padj < 0.05, |log2FC| >= 0.5, baseMean >= 5).")
        return

    expr, expr_note = _metacell_expression(results, cnv_adata)
    if expr is None:
        say("[!] No metacell expression matrix found.")
        return

    obs = mc_obs if mc_obs is not None else (cnv_adata.obs if cnv_adata is not None else None)
    if obs is None or group_key not in obs:
        say(f"[!] Column '{group_key}' not found; skipping.")
        return
    groups = obs[group_key].astype(str).reindex(expr.index)
    keep = groups.notna()
    expr, groups = expr.loc[keep], groups[keep]
    levels = sorted(groups.unique(), key=_lineage_sort_key)

    mal = None
    if "putative_malignant" in obs:
        mal = obs["putative_malignant"].astype(str).reindex(expr.index)

    missing = [g for g in top.index if g not in expr.columns]
    genes = [g for g in top.index if g in expr.columns]
    if missing:
        say(f"[!] {len(missing)} genes were not in the expression matrix: {missing[:6]}")
    if not genes:
        say("[!] No genes left to plot.")
        return

    say(f"Gene selection: rank='{rank}' (top {top_n} among padj < 0.05, |log2FC| >= 0.5, baseMean >= 5)")
    say(f"Expression values: {expr_note}")
    say(f"Groups ({group_key}): " + ", ".join(f"{g} (n={int((groups == g).sum())})" for g in levels))
    if len(levels) == 1:
        say(
            "  [!] Only one group. cell_type is not split by lineage here (v2.3+ "
            "per-cluster labels split into Myeloid_1, Myeloid_2, ...), so this is not "
            "a between-group comparison."
        )
    if mal is not None:
        say("Point color: red = metacell called malignant / blue = normal. Gray = no call.")

    ncol = 5 if len(genes) >= 10 else max(1, min(len(genes), 3))
    nrow = int(np.ceil(len(genes) / ncol))
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(2.25 * ncol + 0.6, 2.15 * nrow + 0.9), squeeze=False
    )
    axes = axes.ravel()
    rng = np.random.default_rng(0)

    for ax, gene in zip(axes, genes):
        vals_all = expr[gene].to_numpy(dtype=float)
        for i, lv in enumerate(levels):
            m = (groups == lv).to_numpy()
            y = vals_all[m]
            if not y.size:
                continue
            xoff = _beeswarm_offsets(y)
            if mal is not None:
                cats = mal.to_numpy()[m]
                colors = np.where(cats == "malignant", C_RED, np.where(cats == "normal", C_BLUE, C_GRAY))
            else:
                colors = np.full(y.size, C_BLUE)
            ax.scatter(
                i + xoff + rng.normal(0, 0.006, y.size), y,
                s=7, c=colors, linewidth=0.15, edgecolor="white", alpha=0.85, zorder=3,
            )
            if mal is not None:
                for sub, col, (x0, x1) in (
                    ("malignant", C_RED, (i - 0.40, i - 0.06)),
                    ("normal", C_BLUE, (i + 0.06, i + 0.40)),
                ):
                    ys = y[cats == sub]
                    if ys.size >= 3:  # a median from n<3 is misleading, so skip it
                        m = float(np.median(ys))
                        ax.plot([x0, x1], [m, m], color=col, lw=2.4, alpha=0.9, zorder=5)
            med = float(np.median(y))
            ax.plot([i - 0.42, i + 0.42], [med, med], color=C_INK, lw=1.2, zorder=6)
        lfc = float(top.loc[gene, "log2FoldChange"])
        padj = float(top.loc[gene, "padj"])
        ptxt = "padj<1e-300" if padj < 1e-300 else f"padj={padj:.1e}"
        ax.set_title(f"{gene}\nlog2FC={lfc:+.2f}  {ptxt}", fontsize=7.8, fontweight="bold")
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels(levels, rotation=45, ha="right", fontsize=7)
        ax.set_xlim(-0.6, len(levels) - 0.4)
        ax.tick_params(axis="y", labelsize=7)
        _despine(ax)
    for ax in axes[len(genes):]:
        ax.axis("off")

    for i, ax in enumerate(axes[: len(genes)]):
        if i % ncol == 0:
            ax.set_ylabel(T("log1p(CPM)", "log1p(CPM)"), fontsize=8)

    handles = [
        matplotlib.lines.Line2D([], [], marker="o", ls="none", ms=5, color=C_RED,
                                label=T("called malignant", "called malignant")),
        matplotlib.lines.Line2D([], [], marker="o", ls="none", ms=5, color=C_BLUE,
                                label=T("called normal", "called normal")),
        matplotlib.lines.Line2D([], [], color=C_RED, lw=2.4, label=T("median of malignant", "median of malignant")),
        matplotlib.lines.Line2D([], [], color=C_BLUE, lw=2.4, label=T("median of normal", "median of normal")),
        matplotlib.lines.Line2D([], [], color=C_INK, lw=1.2, label=T("median of the whole group", "median of the whole group")),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.012), frameon=False)
    fig.suptitle(
        T(f"Fig 9  Top {len(genes)} DEGs across {group_key} groups (metacell level)",
          f"Fig 9  Top {len(genes)} DEGs across {group_key} groups (metacell level)"),
        y=1.0, fontsize=11, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0.015, 1, 0.975))
    out = fig_dir / "09_deg_swarm.png"
    fig.savefig(out, dpi=190, bbox_inches="tight")
    plt.close(fig)
    say(f"-> figures/{out.name}")

    # --- Per-group median table and CSV ---
    med = pd.DataFrame(
        {lv: expr.loc[(groups == lv).to_numpy(), genes].median(axis=0) for lv in levels}
    )
    med.insert(0, "log2FoldChange", top.loc[genes, "log2FoldChange"].round(3))
    med.insert(1, "padj", top.loc[genes, "padj"])
    med.round(3).to_csv(results / "de_swarm_top_genes.csv")
    say("-> de_swarm_top_genes.csv (selected genes and their per-group medians)")
    say("")
    say(
        "How to read this: a gene that rises similarly across all groups may reflect "
        "metacell size or sequencing depth rather than 'malignant vs normal'. "
        "Conversely, a gene that rises only among the points called malignant within a "
        "group suggests tumor cells are mixed into that group."
    )
    if mal is not None and len(levels) > 1:
        cross = pd.crosstab(groups, mal)
        say("")
        say("Group x malignant-call cross table:")
        for line in cross.to_string().split("\n"):
            say("  " + line)
        say(
            "  If a normal lineage (e.g. Myeloid) contains many 'malignant' calls, "
            "suspect false positives in the CNV call. Check this alongside the sign "
            "agreement in figure 8."
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_plotly(args, fig_dir: Path) -> None:
    """Plotly-based figures. Combined into one HTML, with PNGs written too when possible."""
    import metacellcnv_plotly as PL

    res = Path(args.results_dir)
    figs: dict[str, object] = {}
    obs = PL.load_metacell_obs(res)

    prep = Path(args.prep_dir) if args.prep_dir else (res / "scanpy")
    if (prep / "umap3d.tsv.gz").exists() and (res / "cell_to_metacell.csv").exists():
        try:
            figs["Cell embedding"] = PL.fig_scatter_panels(prep, res)
        except Exception as exc:
            warn(f"Skipping scatter panels: {exc}")
    else:
        say(f"{prep}/umap3d.tsv.gz not found; skipping scatter panels")

    qc = prep / "qc_metrics.tsv.gz"
    if qc.exists():
        try:
            figs["Per-cell QC"] = PL.fig_qc(pd.read_csv(qc, sep="\t", index_col=0))
        except Exception as exc:
            warn(f"Skipping the QC figure: {exc}")
    try:
        figs["Metacell quality"] = PL.fig_metacell_quality(obs)
    except Exception as exc:
        warn(f"Skipping the metacell quality figure: {exc}")

    chrm = res / "cnv_chromosome_means.csv"
    if chrm.exists():
        try:
            figs["Chromosome-level CNV"] = PL.fig_cnv_heatmap(
                pd.read_csv(chrm, index_col=0), obs)
        except Exception as exc:
            warn(f"Skipping the CNV heatmap: {exc}")

    nwk = res / "cnv_lineage.nwk"
    if nwk.exists():
        try:
            import json
            import cnv_lineage as LIN
            rep = json.loads((res / "cnv_lineage_report.json").read_text(encoding="utf-8")) \
                if (res / "cnv_lineage_report.json").exists() else {}
            ev = pd.read_csv(res / "cnv_lineage_events.csv", index_col=0).astype(bool)
            root, _ = LIN.dollo_tree(ev, verbose=False)
            figs["Tumour lineage"] = PL.fig_lineage_tree(
                root, has_structure=bool(rep.get("has_structure", False)))
        except Exception as exc:
            warn(f"Skipping the lineage tree figure: {exc}")

    out = fig_dir / "report.html"
    PL.write_report(figs, out, title=f"metacellcnv report — {res.name}",
                    png_dir=fig_dir)
    say(f"Plotly report: {out}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Interpret metacellcnv.py output and generate a report and figures"
    )
    parser.add_argument("--results-dir", default="results", help="Pipeline output directory")
    parser.add_argument("--out-dir", default=None, help="Figure output directory (default: <results>/figures)")
    parser.add_argument("--top-n", type=int, default=25, help="Number of top/bottom genes kept in de_top_genes.csv")
    parser.add_argument(
        "--skip-singlecell",
        action="store_true",
        help="Do not load singlecells_qc.h5ad (large file)",
    )
    parser.add_argument(
        "--chromosome-map",
        default=None,
        help="NCBI assembly report (or a 2-column TSV of accession<TAB>chromosome name). "
        "Turns figure 8's chromosome labels into actual numbers. If omitted, numbers "
        "are inferred from accession order",
    )
    parser.add_argument(
        "--chr-heatmap-clusters",
        type=int,
        default=None,
        help="Number of clusters used to split metacells in figure 8 (default: number of cnv_leiden clusters)",
    )
    parser.add_argument(
        "--include-unplaced",
        action="store_true",
        help="Include unplaced scaffolds such as NW_* in figure 8",
    )
    parser.add_argument(
        "--swarm-top-n", type=int, default=20, help="Number of DEGs shown in the figure 9 swarm plot (default 20)"
    )
    parser.add_argument(
        "--swarm-rank",
        choices=("balanced", "lfc", "padj"),
        default="balanced",
        help="How genes are chosen for figure 9. balanced=half up/half down (default) / lfc=by |log2FC| / padj=by padj",
    )
    parser.add_argument(
        "--engine",
        choices=("plotly", "matplotlib", "both"),
        default="plotly",
        help="Plotting engine (default plotly). plotly combines everything into one HTML "
        "with all-English figure text; matplotlib produces the legacy PNGs",
    )
    parser.add_argument(
        "--prep-dir",
        default=None,
        help="Output of metacellcnv_scanpy (umap3d.tsv.gz etc). "
        "When given, also draws the UMAP scatter panels. Defaults to <results>/scanpy if it exists",
    )
    parser.add_argument(
        "--swarm-group-key",
        default="auto",
        help="metacell column used to group figure 9 (default cell_type; used as-is if split into Myeloid_1, _2 ...)",
    )
    args = parser.parse_args(argv)

    results = Path(args.results_dir)
    if not results.exists():
        print(f"Directory does not exist: {results}", file=sys.stderr)
        return 1
    fig_dir = Path(args.out_dir) if args.out_dir else results / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Swarm group key: when 'auto', prefer anno_label (the evidence-backed label)
    if getattr(args, "swarm_group_key", "auto") == "auto":
        try:
            import metacellcnv_plotly as _PL
            _o = _PL.load_metacell_obs(results)   # also merges annotations from a separate file
            args.swarm_group_key = ("anno_label" if "anno_label" in _o.columns
                                    else "cell_type")
        except Exception:
            args.swarm_group_key = "cell_type"
        say(f"Figure 9 group key: {args.swarm_group_key} (override with --swarm-group-key)")

    if args.engine in ("plotly", "both"):
        try:
            run_plotly(args, fig_dir)
        except Exception as exc:
            warn(f"Plotly rendering failed: {exc}")
        if args.engine == "plotly":
            return 0

    import anndata as ad

    say("metacellcnv output interpretation report")
    say(f"Input: {results.resolve()}")
    say(f"Figure output: {fig_dir.resolve()}")
    say("")
    say("File roles:")
    say("  qc_metrics.csv               Per-cell QC outlier flags (step 1b)")
    say("  metacell_metrics.csv         Metacell quality metrics (step 2)")
    say("  metacell_obs.csv             Metacell attributes (member cell count, cnv_score, malignant label)")
    say("  singlecells_qc.h5ad          Single cells after QC/doublet removal (starting point for re-analysis)")
    say("  metacells.h5ad               Raw counts after metacell aggregation (DE input)")
    say("  cnv_metacells.h5ad           CNV matrix obsm['X_cnv'] and clone classification")
    say("  de_malignant_vs_normal.csv   pyDESeq2 results")
    say("  environment.lock.txt         Package list for reproducibility")

    mc_obs = None
    malignant_labels = None
    p = results / "metacell_obs.csv"
    if p.exists():
        mc_obs = pd.read_csv(p, index_col=0)
        if "putative_malignant" in mc_obs:
            malignant_labels = mc_obs["putative_malignant"].astype(str)
    call_note = results / "malignant_call.txt"
    if call_note.exists():
        say("")
        say("Rationale for the malignant call (malignant_call.txt):")
        say("  " + call_note.read_text(encoding="utf-8").strip())

    sc_adata = None
    if not args.skip_singlecell and (results / "singlecells_qc.h5ad").exists():
        try:
            sc_adata = ad.read_h5ad(results / "singlecells_qc.h5ad")
        except Exception as exc:
            say(f"[!] Failed to load singlecells_qc.h5ad: {exc}")

    report_qc(results, fig_dir, sc_adata)
    report_metacell_quality(results, fig_dir, mc_obs)

    cnv_adata = None
    if (results / "cnv_metacells.h5ad").exists():
        try:
            cnv_adata = ad.read_h5ad(results / "cnv_metacells.h5ad")
            report_cnv(cnv_adata, fig_dir, malignant_labels)
        except Exception as exc:
            say(f"[!] Failed to process cnv_metacells.h5ad: {exc}")

    report_mito_qc(results, fig_dir, mc_obs)
    report_de(results, fig_dir, mc_obs, top_n=args.top_n)

    if cnv_adata is not None:
        try:
            report_chromosome_cnv_heatmap(
                cnv_adata,
                fig_dir,
                results,
                malignant_labels=malignant_labels,
                chromosome_map_path=Path(args.chromosome_map) if args.chromosome_map else None,
                n_clusters=args.chr_heatmap_clusters,
                exclude_unplaced=not args.include_unplaced,
            )
        except Exception as exc:
            say(f"[!] Failed to build figure 8 (chromosome-level CNV heatmap): {exc}")
            plt.close("all")
    try:
        report_deg_swarm(
            results,
            fig_dir,
            cnv_adata,
            mc_obs,
            top_n=args.swarm_top_n,
            rank=args.swarm_rank,
            group_key=args.swarm_group_key,
        )
    except Exception as exc:
        say(f"[!] Failed to build figure 9 (DEG swarm plot): {exc}")
        plt.close("all")

    report_path = results / "interpretation_report.txt"
    report_path.write_text("\n".join(REPORT) + "\n", encoding="utf-8")
    print()
    print(f"[interpret] Report saved: {report_path}")
    print(f"[interpret] Figures: {sorted(p.name for p in fig_dir.glob('*.png'))}")
    return 0


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    sys.exit(main())
