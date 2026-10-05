#!/usr/bin/env python3
"""Draw the "how are normal metacells decided" figure (guide Fig. 6 and 7).

Inputs (from one pipeline run):
  --cnv-h5ad  cnv_metacells.h5ad   (obsm['X_cnv'], uns['cnv']['chr_pos'])
  --obs       metacell_obs.csv     (cell_type, cnv_score, cnv_leiden, putative_malignant)

If anndata cannot read the h5ad (null uns/log1p), no fix is needed here:
this script reads the file with h5py directly.

Example:
  python docs/make_call_figure.py --cnv-h5ad results/Case1/cnv_metacells.h5ad \
      --obs results/Case1/metacell_obs.csv --name Case1 --out-dir docs/img/en --lang en
"""
import argparse
import os
import re

import h5py
import matplotlib
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd
import scipy.sparse as sp

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Cell types used as the normal reference (same as the marker tables).
REFERENCE_TYPES = ["T/NK", "B", "Plasma", "Myeloid", "Macrophage", "Mast", "Endothelial", "SmoothMuscle", "Leukocyte"]
LAG = 10  # windows are 100 genes wide with step 10, so a shift of 10 windows does not overlap
RED, BLUE, GREEN = "#d62728", "#1f77b4", "#2ca02c"

TXT = {
    "ja": {
        "ref": "正常参照(既知の正常細胞型, n={})", "other": "それ以外(n={})",
        "thr": " 閾値=平均+3SD\n={:.2f}", "xa": "cnv_score(CNVプロファイルのL2ノルム)", "ya": "metacell数",
        "ta": "(a) 判定の基準:参照の cnv_score 分布から閾値を決める",
        "xb": "CNVクローン(cnv_leiden)", "yb": "クローン内 cnv_score 中央値", "tb": "(b) クローン単位で判定(赤=悪性 / 青=正常)",
        "yc": "ブロック性(10窓ずらしの自己相関)", "tc": "(c) スコアとは独立な確認:染色体内で値が連続しているか",
        "call": "{}コール(参照以外)", "refl": "正常参照", "side": "判定",
        "td": "(d) CNVプロファイル(行=metacell、クローン順。左の帯=判定、緑点=正常参照、点線=推定CNV変化点)",
        "te": "(e) 平均CNVプロファイル(悪性コール群 − 正常コール群 から変化点を推定)",
        "mal_mean": "悪性コール平均", "nor_mean": "正常コール平均",
        "tf": "(f) 隣接(非重複)窓の発現由来CNV値のmetacell間相関:全体 / 正常コール / 悪性コール(平滑化)",
        "all": "全metacell", "nor": "正常コール群", "mal": "悪性コール群", "r": "相関 r",
        "tg": "(g) クローン別の隣接窓相関(行=クローン、悪性=赤字/正常=青字、8 metacell以上)",
        "th": "(h) 各metacellと悪性コール群平均プロファイルの相関(クローン別)", "yh": "相関",
        "ti": "(i) ブロック性(クローン別)", "yi": "10窓ずらし自己相関",
        "xhi": "CNVクローン(赤=悪性コール/青=正常コール)", "tj": "(j) 同・細胞型(サブクラスタ)別",
        "title": "{}:「正常」metacell の決め方と、判定が CNV の構造に支えられているか",
    },
    "en": {
        "ref": "Normal reference (known normal cell types, n={})", "other": "Others (n={})",
        "thr": " threshold = mean + 3SD\n= {:.2f}", "xa": "cnv_score (L2 norm of the CNV profile)", "ya": "metacells",
        "ta": "(a) Rule: the threshold comes from the reference cnv_score distribution",
        "xb": "CNV clone (cnv_leiden)", "yb": "median cnv_score in clone", "tb": "(b) Call per clone (red = malignant, blue = normal)",
        "yc": "blockiness (autocorrelation at 10-window shift)", "tc": "(c) Check independent of the score: are values continuous along chromosomes?",
        "call": "{} call (non-reference)", "refl": "Normal reference", "side": "call",
        "td": "(d) CNV profile (rows = metacells in clone order; left bar = call, green dot = reference, dotted = estimated CNV breakpoint)",
        "te": "(e) Mean CNV profile (breakpoints estimated from malignant-call minus normal-call mean)",
        "mal_mean": "malignant-call mean", "nor_mean": "normal-call mean",
        "tf": "(f) Correlation across metacells between non-overlapping windows: all / normal call / malignant call (smoothed)",
        "all": "all metacells", "nor": "normal-call group", "mal": "malignant-call group", "r": "correlation r",
        "tg": "(g) Adjacent-window correlation per clone (rows = clones, red = malignant, blue = normal; >= 8 metacells)",
        "th": "(h) Correlation of each metacell with the malignant-call mean profile (per clone)", "yh": "correlation",
        "ti": "(i) Blockiness (per clone)", "yi": "autocorrelation at 10-window shift",
        "xhi": "CNV clone (red = malignant call, blue = normal call)", "tj": "(j) Same, per cell type (sub-cluster)",
        "title": "{}: how \"normal\" metacells are decided, and whether the call is backed by CNV structure",
    },
}


def set_font():
    avail = {f.name for f in fm.fontManager.ttflist}
    for name in ["Noto Sans CJK JP", "Noto Serif CJK JP", "IPAexGothic", "IPAGothic", "Hiragino Sans", "Yu Gothic"]:
        if name in avail:
            plt.rcParams["font.family"] = [name, "DejaVu Sans"]
            break
    plt.rcParams["axes.unicode_minus"] = False


def load(h5_path, obs_path):
    """Dense CNV matrix, chromosome offsets and the metacell table (same row order)."""
    f = h5py.File(h5_path, "r")
    g = f["obsm/X_cnv"]
    idx = f["obs"].attrs["_index"] if "_index" in f["obs"].attrs else "_index"
    names = [x.decode() if isinstance(x, bytes) else x for x in f["obs"][idx][()]]
    X = sp.csr_matrix((g["data"][()], g["indices"][()], g["indptr"][()]), shape=tuple(g.attrs["shape"])).toarray()
    offsets = {k: int(v[()]) for k, v in f["uns/cnv/chr_pos"].items()}
    obs = pd.read_csv(obs_path, index_col=0).loc[names]
    return X, offsets, obs


def chromosome_ranges(offsets, n_windows):
    """(name, start, end) of autosomes and X; mitochondria and unplaced scaffolds are dropped."""
    items = sorted(offsets.items(), key=lambda t: t[1])
    out = []
    for i, (name, start) in enumerate(items):
        end = items[i + 1][1] if i + 1 < len(items) else n_windows
        if re.fullmatch(r"chr(\d+|X)", name):
            out.append((name, start, end))
    return out


def lag_corr(Xg, ranges):
    """Pearson correlation across metacells between window i and window i+LAG (within chromosomes)."""
    r = np.full(Xg.shape[1], np.nan)
    if Xg.shape[0] < 5:
        return r
    sd = Xg.std(0)
    z = (Xg - Xg.mean(0)) / np.where(sd > 1e-9, sd, np.nan)
    for _, a, b in ranges:
        if b - a > LAG:
            r[a:b - LAG] = np.nanmean(z[:, a:b - LAG] * z[:, a + LAG:b], axis=0)
    return r


def find_breakpoints(d, a, b, min_len=40, gain=0.15, shift=0.03, depth=3):
    """Binary segmentation of one chromosome; returns breakpoint window indices."""
    out = []

    def rec(lo, hi, dep):
        if dep == 0 or hi - lo < 2 * min_len:
            return
        x = d[lo:hi]
        total = ((x - x.mean()) ** 2).sum()
        if total <= 0:
            return
        cs, n, best = np.cumsum(x), len(x), None
        for k in range(min_len, n - min_len + 1):
            m1, m2 = cs[k - 1] / k, (cs[-1] - cs[k - 1]) / (n - k)
            sse = ((x[:k] - m1) ** 2).sum() + ((x[k:] - m2) ** 2).sum()
            if best is None or sse < best[0]:
                best = (sse, k, abs(m1 - m2))
        if best and (total - best[0]) / total > gain and best[2] > shift:
            out.append(lo + best[1])
            rec(lo, lo + best[1], dep - 1)
            rec(lo + best[1], hi, dep - 1)

    rec(a, b, depth)
    return sorted(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cnv-h5ad", required=True)
    ap.add_argument("--obs", required=True)
    ap.add_argument("--name", required=True, help="label used in the title and the file names")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--lang", choices=["ja", "en"], default="en")
    args = ap.parse_args()
    T = TXT[args.lang]
    set_font()
    os.makedirs(args.out_dir, exist_ok=True)

    X, offsets, obs = load(args.cnv_h5ad, args.obs)
    ranges = chromosome_ranges(offsets, X.shape[1])
    keep = np.concatenate([np.arange(a, b) for _, a, b in ranges])
    X = X[:, keep]  # keep analysed windows only
    pos, packed = 0, []
    for name, a, b in ranges:
        packed.append((name, pos, pos + (b - a)))
        pos += b - a
    ranges, W = packed, X.shape[1]

    ct = obs["cell_type"].astype(str)
    is_ref = np.array([any(c == k or c.startswith(k + "_") for k in REFERENCE_TYPES) for c in ct])
    call = obs["putative_malignant"].astype(str).values
    clone = obs["cnv_leiden"].astype(str).values
    score = obs["cnv_score"].astype(float).values
    mu, sd = score[is_ref].mean(), score[is_ref].std(ddof=1)
    thr = mu + 3 * sd
    clone_med = pd.Series(score).groupby(clone).median().sort_values()
    is_mal, is_nor = call == "malignant", call == "normal"

    # Breakpoints from the malignant-minus-normal mean profile.
    diff = X[is_mal].mean(0) - X[is_nor].mean(0)
    bp = np.array(sorted(b_ for _, a, b in ranges for b_ in find_breakpoints(diff, a, b)))

    # Blockiness: autocorrelation of each metacell's own profile (10-window shift, within chromosomes).
    blk = np.zeros(len(X))
    for i in range(len(X)):
        num = den = 0.0
        for _, a, b in ranges:
            x = X[i, a:b]
            num += (x[:-LAG] * x[LAG:]).sum()
            den += (x[:-LAG] ** 2).sum() + 1e-12
        blk[i] = num / den

    # Correlation of each metacell with the malignant-call mean profile.
    tp = X[is_mal].mean(0)
    Xc, tpc = X - X.mean(1, keepdims=True), tp - tp.mean()
    tcor = (Xc @ tpc) / (np.linalg.norm(Xc, axis=1) * np.linalg.norm(tpc) + 1e-12)

    r_all, r_nor, r_mal = lag_corr(X, ranges), lag_corr(X[is_nor], ranges), lag_corr(X[is_mal], ranges)
    big = [c for c in clone_med.index if (clone == c).sum() >= 8]
    r_clone = {c: lag_corr(X[clone == c], ranges) for c in big}

    pd.DataFrame(dict(metacell=obs.index, cell_type=ct.values, clone=clone, cnv_score=score, call=call,
                      is_reference=is_ref, blockiness=blk, tumor_pattern_corr=tcor)
                 ).to_csv(os.path.join(args.out_dir, f"{args.name}_metacell_call_detail.csv"), index=False)
    print(f"{args.name}: metacells={len(X)} windows={W} reference={is_ref.sum()} "
          f"mean={mu:.3f} sd={sd:.3f} threshold={thr:.3f} breakpoints={len(bp)}")
    print(clone_med.round(3).to_string())

    # ---------------- figure ----------------
    order = np.lexsort((score, np.array([clone_med[c] for c in clone])))
    fig = plt.figure(figsize=(17, 21))
    gs = fig.add_gridspec(6, 3, height_ratios=[1.5, 3.2, 1.2, 1.2, 1.7, 1.6], hspace=0.45, wspace=0.28)
    col_of = {"malignant": RED, "normal": BLUE}
    chrom_ticks = [(a + b) / 2 for _, a, b in ranges]
    chrom_names = [n.replace("chr", "") for n, _, _ in ranges]

    ax = fig.add_subplot(gs[0, 0])  # (a) threshold rule
    bins = np.linspace(0, score.max() * 1.02, 40)
    ax.hist(score[is_ref], bins, color=GREEN, alpha=.8, label=T["ref"].format(is_ref.sum()))
    ax.hist(score[~is_ref], bins, color="#999", alpha=.7, label=T["other"].format((~is_ref).sum()))
    ax.axvline(thr, color="k", ls="--")
    ax.text(thr, ax.get_ylim()[1] * .92, T["thr"].format(thr), fontsize=8, va="top")
    ax.set_xlabel(T["xa"]); ax.set_ylabel(T["ya"]); ax.legend(fontsize=7.5, loc="upper right")
    ax.set_title(T["ta"], fontsize=9.5, loc="left")

    ax = fig.add_subplot(gs[0, 1])  # (b) clone medians
    xs = np.arange(len(clone_med))
    ax.bar(xs, clone_med.values, color=[col_of[call[clone == c][0]] for c in clone_med.index])
    ax.axhline(thr, color="k", ls="--")
    for x_, c in zip(xs, clone_med.index):
        m = clone == c
        comp = pd.Series(ct.values[m]).str.replace(r"_\d+$", "", regex=True).value_counts()
        ax.text(x_, clone_med[c], f"n={m.sum()}\n" + "/".join(f"{k}{v}" for k, v in comp.head(2).items()),
                ha="center", va="bottom", fontsize=5.5)
    ax.set_xticks(xs); ax.set_xticklabels(clone_med.index)
    ax.set_xlabel(T["xb"]); ax.set_ylabel(T["yb"]); ax.set_title(T["tb"], fontsize=9.5, loc="left")

    ax = fig.add_subplot(gs[0, 2])  # (c) score vs blockiness
    for g_, col in col_of.items():
        m = (call == g_) & ~is_ref
        ax.scatter(score[m], blk[m], s=14, c=col, alpha=.6, label=T["call"].format(g_))
    ax.scatter(score[is_ref], blk[is_ref], s=14, c=GREEN, alpha=.7, label=T["refl"])
    ax.axvline(thr, color="k", ls="--")
    ax.set_xlabel("cnv_score"); ax.set_ylabel(T["yc"]); ax.legend(fontsize=7.5)
    ax.set_title(T["tc"], fontsize=9.5, loc="left")

    sub = gs[1, :].subgridspec(1, 2, width_ratios=[0.025, 1], wspace=0.01)  # (d) heatmap
    ax_side, ax_heat = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
    lim = np.percentile(np.abs(X), 99)
    im = ax_heat.imshow(X[order], aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")
    side = np.array([[0 if call[i] == "normal" else 1] for i in order])
    ax_side.imshow(side, aspect="auto", cmap=ListedColormap([BLUE, RED]), interpolation="nearest")
    ax_side.set_xticks([]); ax_side.set_yticks([]); ax_side.set_title(T["side"], fontsize=7)
    for yi, i in enumerate(order):
        if is_ref[i]:
            ax_side.plot(-0.6, yi, marker=".", color=GREEN, ms=2, clip_on=False)
    clone_sorted = np.array([clone[i] for i in order])
    for y in np.where(clone_sorted[1:] != clone_sorted[:-1])[0] + .5:
        ax_heat.axhline(y, color="k", lw=.6)
    for c in clone_med.index:
        ax_heat.text(W * 1.003, np.where(clone_sorted == c)[0].mean(), f"clone {c}", fontsize=7, va="center")
    for _, a, _ in ranges:
        ax_heat.axvline(a - .5, color="#888", lw=.4)
    for b_ in bp:
        ax_heat.axvline(b_, color="k", lw=.7, ls=":", alpha=.8)
    ax_heat.set_xticks(chrom_ticks); ax_heat.set_xticklabels(chrom_names, fontsize=7); ax_heat.set_yticks([])
    ax_heat.set_title(T["td"], fontsize=9.5, loc="left")
    fig.colorbar(im, cax=fig.add_axes([0.945, 0.62, 0.008, 0.1])).ax.tick_params(labelsize=7)

    def track(row, title):
        a_ = fig.add_subplot(gs[row, :])
        a_.set_xlim(-.5, W - .5)
        for _, a, _ in ranges:
            a_.axvline(a - .5, color="#ccc", lw=.4)
        for b_ in bp:
            a_.axvline(b_, color="k", lw=.7, ls=":", alpha=.8)
        a_.set_xticks(chrom_ticks); a_.set_xticklabels(chrom_names, fontsize=7)
        a_.set_title(title, fontsize=9.5, loc="left")
        return a_

    def smooth(r):
        return pd.Series(r).rolling(15, min_periods=5, center=True).mean().values

    a_ = track(2, T["te"])  # (e) mean profiles
    a_.plot(X[is_mal].mean(0), c=RED, lw=.8, label=T["mal_mean"])
    a_.plot(X[is_nor].mean(0), c=BLUE, lw=.8, label=T["nor_mean"])
    a_.legend(fontsize=7, ncol=2, loc="upper right")

    a_ = track(3, T["tf"])  # (f) adjacent-window correlation
    a_.plot(smooth(r_all), c="k", lw=1, label=T["all"])
    a_.plot(smooth(r_nor), c=BLUE, lw=1, label=T["nor"])
    a_.plot(smooth(r_mal), c=RED, lw=1, label=T["mal"])
    a_.axhline(0, c="#888", lw=.5); a_.set_ylim(-.5, 1); a_.legend(fontsize=7, ncol=3, loc="upper right")
    a_.set_ylabel(T["r"], fontsize=8)

    a_ = track(4, T["tg"])  # (g) per-clone correlation
    M = np.vstack([smooth(r_clone[c]) for c in big])
    im2 = a_.imshow(M, aspect="auto", cmap="viridis", vmin=-.3, vmax=1, extent=(-.5, W - .5, len(big) - .5, -.5),
                    interpolation="nearest")
    a_.set_yticks(range(len(big)))
    a_.set_yticklabels([f"clone {c} (n={int((clone == c).sum())})" for c in big], fontsize=7)
    for t, c in zip(a_.get_yticklabels(), big):
        t.set_color(col_of[call[clone == c][0]])
    fig.colorbar(im2, cax=fig.add_axes([0.945, 0.30, 0.008, 0.07])).ax.tick_params(labelsize=7)

    clones = list(clone_med.index)  # (h)(i) boxplots per clone, (j) per cell type
    for j, (vals, ttl, yl) in enumerate([(tcor, T["th"], T["yh"]), (blk, T["ti"], T["yi"])]):
        a2 = fig.add_subplot(gs[5, j])
        a2.boxplot([vals[clone == c] for c in clones], labels=clones, patch_artist=True,
                   flierprops=dict(ms=2), boxprops=dict(facecolor="#ddd"))
        for t, c in zip(a2.get_xticklabels(), clones):
            t.set_color(col_of[call[clone == c][0]])
        a2.set_title(ttl, fontsize=9.5, loc="left"); a2.set_xlabel(T["xhi"], fontsize=8)
        a2.set_ylabel(yl, fontsize=8); a2.axhline(0, c="#888", lw=.5)
    a2 = fig.add_subplot(gs[5, 2])
    subs = sorted(set(ct.values))
    a2.boxplot([tcor[ct.values == c] for c in subs], labels=subs, patch_artist=True,
               flierprops=dict(ms=2), boxprops=dict(facecolor="#cfe8cf"))
    a2.tick_params(axis="x", labelrotation=60, labelsize=7)
    a2.set_title(T["tj"], fontsize=9.5, loc="left"); a2.axhline(0, c="#888", lw=.5)

    fig.suptitle(T["title"].format(args.name), fontsize=13, y=0.93)
    fig.savefig(os.path.join(args.out_dir, f"{args.name}_normal_call.png"), dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
