#!/usr/bin/env python3
"""Draw Figs. 3-5 of the guide from the validation CSVs in docs/data (white background).

Usage (from the repository root):
    python docs/make_figures.py --lang ja   # -> docs/img/
    python docs/make_figures.py --lang en   # -> docs/img/en/
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))

T = {
    "ja": dict(
        groups=["外部参照のみ\n(公開データ)", "外部参照 +\n遺伝子性質の補正", "同一実験の\n他株の中央値"],
        schemes=["位置基準の窓", "位置基準 + 遺伝子数≥50% の窓"],
        f3_y="DepMap の相対コピー数との相関 r\n(9 株・20 Mb 窓)",
        f4_labels=None,
        f4_thr="CNV 判定の閾値 0.3",
        f4_x="差プロファイルの標準偏差(20 Mb 窓, log2)",
        f5_x="細胞数(粗クラスタに付いた細胞型ラベルごと)\n[ ] 内は、参照型マーカー遺伝子のうちデータに存在した数/全数",
        f5_legend="付いたラベル",
        f5_note="参照型(T/NK・B・Myeloid など)と判定された細胞\n0 / {tot}",
    ),
    "en": dict(
        groups=["External reference only\n(public data)", "External reference +\ngene-property correction",
                "Median of other lines\nin the same experiment"],
        schemes=["Position-based windows", "Position-based + gene count ≥50%"],
        f3_y="Correlation r with DepMap relative copy number\n(9 lines, 20 Mb windows)",
        f4_labels=["Within sample, same type (random halves)",
                   "Within sample, same type (Case23 Myeloid halves)",
                   "Within sample, same type (Case23 Other halves)",
                   "Within sample, different type (Case23 Other vs Myeloid)",
                   "Common model (Case1 Myeloid vs dog PBMC)",
                   "Common model (Case23 Myeloid vs dog PBMC)",
                   "Common model (Case23 Other vs dog PBMC)",
                   "Across samples, different type (Case23 Other vs Case1 Myeloid)",
                   "Across samples, same type (Case23 Myeloid vs Case1 Myeloid)"],
        f4_thr="CNV call threshold 0.3",
        f4_x="SD of the difference profile (20 Mb windows, log2)",
        f5_x="Cells (by cell-type label given to the coarse clusters)\n[ ] = reference-type marker genes found in the data / total",
        f5_legend="Assigned label",
        f5_note="Cells called a reference type (T/NK, B, Myeloid, ...)\n0 / {tot}",
    ),
}


def set_style():
    avail = {f.name for f in fm.fontManager.ttflist}
    for n in ["Noto Sans CJK JP", "IPAexGothic", "IPAGothic", "Hiragino Sans", "Yu Gothic", "Meiryo"]:
        if n in avail:
            plt.rcParams["font.family"] = [n, "DejaVu Sans"]
            break
    plt.rcParams.update({"axes.unicode_minus": False, "figure.facecolor": "white", "axes.facecolor": "white",
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 10})


def fig3(t, D, out):
    """Correlation with DepMap copy number, by reference type (20 Mb windows)."""
    pv = pd.read_csv(D("pos_vs_count.csv"))
    pv = pv[pv.W == 20]
    profiles = ["A0", "A-共変量補正", "B"]
    schemes = [("位置基準(Mb固定, ≥3遺伝子)", "#9ecae1"),
               ("位置+遺伝子数(Mb固定 かつ ≥50%遺伝子)", "#08519c")]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    w = 0.36
    for j, (sch, col) in enumerate(schemes):
        vals = [float(pv[(pv.profile == g) & (pv.scheme == sch)].r.iloc[0]) for g in profiles]
        xs = [i + (j - 0.5) * w for i in range(len(profiles))]
        ax.bar(xs, vals, w, color=col, label=t["schemes"][j])
        for x, v in zip(xs, vals):
            ax.text(x, v + 0.01, f"{v:.2f}", ha="center", fontsize=9)
    ax.set_xticks(range(len(profiles)))
    ax.set_xticklabels(t["groups"])
    ax.set_ylabel(t["f3_y"])
    ax.set_ylim(0, 0.85)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig3_reference_accuracy.png"), dpi=150)
    plt.close(fig)


def fig4(t, D, out):
    """SD of the difference profile between normal cell sets."""
    cc = pd.read_csv(D("cross_case2_null.csv"))
    sd_col = [c for c in cc.columns if "20Mb" in c and "標準偏差" in c][0]
    # (CSV row, bar colour), top to bottom
    rows = [(0, "#7f7f7f"), (1, "#7f7f7f"), (2, "#7f7f7f"), (5, "#bdbdbd"), (7, "#fd8d3c"),
            (8, "#fd8d3c"), (9, "#fd8d3c"), (6, "#d62728"), (3, "#d62728")]
    if t["f4_labels"] is None:
        labels = [cc.iloc[i, 0] for i, _ in rows]
    else:
        labels = t["f4_labels"]   # already in plotting order
    vals = [float(cc.iloc[i][sd_col]) for i, _ in rows]
    cols = [c for _, c in rows]
    fig, ax = plt.subplots(figsize=(9.6, 4.4))
    ax.barh(labels[::-1], vals[::-1], color=cols[::-1])
    for y, v in enumerate(vals[::-1]):
        ax.text(v + 0.01, y, f"{v:.2f}", va="center", fontsize=9)
    ax.axvline(0.3, color="k", ls="--", lw=1)
    ax.text(0.31, len(rows) - 0.55, t["f4_thr"], fontsize=8.5, va="bottom")
    ax.set_xlabel(t["f4_x"])
    ax.set_xlim(0, 0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig4_cross_sample_null.png"), dpi=150)
    plt.close(fig)


def fig5(t, D, out):
    """Cells called a reference type in tumor-only samples."""
    cl = pd.read_csv(D("ref_count_clusters.csv"))
    tab = cl.pivot_table(index="sample", columns="label", values="n_cells", aggfunc="sum", fill_value=0)
    tab = tab.reindex(columns=[c for c in ["Epithelial", "Fibroblast", "Other"] if c in tab.columns])
    colors = {"Epithelial": "#6baed6", "Fibroblast": "#fd8d3c", "Other": "#bdbdbd"}
    cov = pd.read_csv(D("geo_marker_coverage.csv")).set_index("sample")
    tab.index = [f"{n}  [{int(cov.loc[n].ref_panel_genes_present)}/{int(cov.loc[n].ref_panel_genes_total)}]"
                 for n in tab.index]
    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    left = np.zeros(len(tab))
    for c in tab.columns:
        ax.barh(tab.index, tab[c].values, left=left, color=colors[c], label=c)
        left += tab[c].values
    tot = f"{int(tab.values.sum()):,}"
    ax.set_xlabel(t["f5_x"])
    ax.legend(frameon=False, title=t["f5_legend"], loc="lower right")
    ax.text(0.98, 0.5, t["f5_note"].format(tot=tot), transform=ax.transAxes, ha="right",
            fontsize=10, color="#d62728")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig5_tumor_only_reference_count.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lang", choices=["ja", "en"], default="ja")
    a = ap.parse_args()
    D = lambda f: os.path.join(HERE, "data", f)
    out = os.path.join(HERE, "img") if a.lang == "ja" else os.path.join(HERE, "img", "en")
    os.makedirs(out, exist_ok=True)
    set_style()
    t = T[a.lang]
    fig3(t, D, out)
    fig4(t, D, out)
    fig5(t, D, out)
    print("wrote figures to", out)
