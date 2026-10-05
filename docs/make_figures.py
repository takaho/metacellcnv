#!/usr/bin/env python3
"""docs/data の検証結果 CSV から、ガイドの図 3〜5 を作る(白背景)。
実行: python docs/make_figures.py   (リポジトリのルートで)
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))
D = lambda f: os.path.join(HERE, "data", f)
OUT = lambda f: os.path.join(HERE, "img", f)
avail = {f.name for f in fm.fontManager.ttflist}
for n in ["Noto Sans CJK JP", "Noto Serif CJK JP", "IPAexGothic", "IPAGothic", "Hiragino Sans", "Yu Gothic", "Meiryo"]:
    if n in avail:
        plt.rcParams["font.family"] = [n, "DejaVu Sans"]
        break
plt.rcParams.update({"axes.unicode_minus": False, "figure.facecolor": "white", "axes.facecolor": "white",
                     "axes.spines.top": False, "axes.spines.right": False, "font.size": 10})

# ---- 図3: DepMap を正解にした CNV 検出の相関(20 Mb 窓) ----
pv = pd.read_csv(D("pos_vs_count.csv"))
pv = pv[pv.W == 20]
groups = [("A0", "外部参照のみ\n(公開データ)"), ("A-共変量補正", "外部参照 +\n遺伝子性質の補正"), ("B", "同一実験の\n他株の中央値")]
schemes = [("位置基準(Mb固定, ≥3遺伝子)", "位置基準の窓", "#9ecae1"), ("位置+遺伝子数(Mb固定 かつ ≥50%遺伝子)", "位置基準 + 遺伝子数≥50% の窓", "#08519c")]
fig, ax = plt.subplots(figsize=(7.2, 3.8))
w = 0.36
for j, (sch, lab, col) in enumerate(schemes):
    vals = [float(pv[(pv.profile == g) & (pv.scheme == sch)].r.iloc[0]) for g, _ in groups]
    xs = [i + (j - 0.5) * w for i in range(len(groups))]
    bars = ax.bar(xs, vals, w, color=col, label=lab)
    for x, v in zip(xs, vals):
        ax.text(x, v + 0.01, f"{v:.2f}", ha="center", fontsize=9)
ax.set_xticks(range(len(groups)))
ax.set_xticklabels([g[1] for g in groups])
ax.set_ylabel("DepMap の相対コピー数との相関 r\n(9 株・20 Mb 窓)")
ax.set_ylim(0, 0.85)
ax.legend(frameon=False, loc="upper left", fontsize=9)
fig.tight_layout()
fig.savefig(OUT("fig3_reference_accuracy.png"), dpi=150)
plt.close(fig)

# ---- 図4: 正常細胞どうしの差プロファイルの大きさ(20 Mb 窓) ----
cc = pd.read_csv(D("cross_case2_null.csv"))
sd_col = [c for c in cc.columns if "20Mb" in c and "標準偏差" in c][0]
name_col = cc.columns[0]
rows = [(cc.iloc[i][name_col], i, col) for i, col in
        [(0, "#7f7f7f"), (1, "#7f7f7f"), (2, "#7f7f7f"), (5, "#bdbdbd"), (7, "#fd8d3c"), (8, "#fd8d3c"),
         (9, "#fd8d3c"), (6, "#d62728"), (3, "#d62728")]]
fig, ax = plt.subplots(figsize=(9.6, 4.4))
labels = [r[0] for r in rows][::-1]
vals = [float(cc.iloc[r[1]][sd_col]) for r in rows][::-1]
cols = [r[2] for r in rows][::-1]
ax.barh(labels, vals, color=cols)
for y, v in enumerate(vals):
    ax.text(v + 0.01, y, f"{v:.2f}", va="center", fontsize=9)
ax.axvline(0.3, color="k", ls="--", lw=1)
ax.text(0.31, len(rows) - 0.55, "CNV 判定の閾値 0.3", fontsize=8.5, va="bottom")
ax.set_xlabel("差プロファイルの標準偏差(20 Mb 窓, log2)")
ax.set_xlim(0, 0.8)
fig.tight_layout()
fig.savefig(OUT("fig4_cross_sample_null.png"), dpi=150)
plt.close(fig)

# ---- 図5: 腫瘍のみのサンプルで「参照」と判定された細胞 ----
cl = pd.read_csv(D("ref_count_clusters.csv"))
t = cl.pivot_table(index="sample", columns="label", values="n_cells", aggfunc="sum", fill_value=0)
t = t.reindex(columns=[c for c in ["Epithelial", "Fibroblast", "Other"] if c in t.columns])
colors = {"Epithelial": "#6baed6", "Fibroblast": "#fd8d3c", "Other": "#bdbdbd"}
cov = pd.read_csv(D("geo_marker_coverage.csv")).set_index("sample")
t.index = [f"{n}  [{int(cov.loc[n].ref_panel_genes_present)}/{int(cov.loc[n].ref_panel_genes_total)}]" for n in t.index]
fig, ax = plt.subplots(figsize=(8.4, 4.4))
left = None
import numpy as np
left = np.zeros(len(t))
for c in t.columns:
    ax.barh(t.index, t[c].values, left=left, color=colors[c], label=c)
    left += t[c].values
tot = int(t.values.sum())
ax.set_xlabel("細胞数(粗クラスタに付いた細胞型ラベルごと)\n[ ] 内は、参照型マーカー遺伝子のうちデータに存在した数/全数")
ax.legend(frameon=False, title="付いたラベル", loc="lower right")
ax.text(0.98, 0.5, f"参照型(T/NK・B・Myeloid など)と判定された細胞\n0 / {tot:,}", transform=ax.transAxes, ha="right", fontsize=10, color="#d62728")
fig.tight_layout()
fig.savefig(OUT("fig5_tumor_only_reference_count.png"), dpi=150)
plt.close(fig)
print("done")
