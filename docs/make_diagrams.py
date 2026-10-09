#!/usr/bin/env python3
"""Draw Fig. 1 (pipeline flow) and Fig. 2 (choosing a reference) as SVG.

Usage (from the repository root):
    python docs/make_diagrams.py --lang ja   # -> docs/img/
    python docs/make_diagrams.py --lang en   # -> docs/img/en/
"""
import argparse
import os
from xml.sax.saxutils import escape as esc

HERE = os.path.dirname(os.path.abspath(__file__))
FONT = "'Noto Sans CJK JP','Hiragino Sans','Yu Gothic','Meiryo','Helvetica Neue',Arial,sans-serif"

# --- Text in each language ---------------------------------------------------
# Fig. 1: (title, description lines, option lines)
STAGES = {
    "ja": [
        ("① 入力", ["Cell Ranger の filtered_feature_bc_matrix", "+ GTF + 種別のマーカー表"],
         ["--cellranger-dir  --sample-id  --gtf  --gtf-gene-id auto",
          "--chromosome-map  --mito-chromosome  --markers {dog|human|mouse|CSV}",
          "--preflight-only(入力だけ検査して終了)"]),
        ("② 前処理・QC", ["HVG/PCA・粗クラスタ・細胞型注釈", "mtDNA 構成検査・MAD QC・doublet"],
         ["--n-hvg  --n-pcs  --leiden-resolution",
          "--marker-rule  --marker-min-margin  --marker-min-within-z",
          "--no-pctmt-filter  --expected-doublet-rate  --skip-doublet-rescue"]),
        ("③ metacell 構築", ["近い細胞を束ねてノイズを減らす", "(SEACells 非依存の自前実装)"],
         ["--cells-per-metacell  --metacell-seed",
          "--seacell-assignments auto(既存の割り当てを再利用)"]),
        ("④ CNV 推定", ["参照との差を窓移動平均で CNV にする", "参照の選び方が結果の意味を決める(図2)"],
         ["--cnv-reference {celltype|none|external†}",
          "--normal-celltype  --normal-clusters  --cnv-refine",
          "--cnv-method {infercnv|bins}  --cnv-external-reference†"]),
        ("⑤ 悪性の判定", ["CNV クローン(cnv_leiden)単位で", "悪性 / 正常を決める"],
         ["--malignant-call {auto|reference|gap|median}",
          "--malignant-n-sd 3.0  --min-clone-gap 0.25"]),
        ("⑥ 型注釈・系譜・DE", ["metacell 注釈(anno_label)・系譜(Dollo)・DE"],
         ["--lineage  --lineage-level  --lineage-events",
          "--no-metacell-annotation(anno_label は既定で付く。§6.4)",
          "--sample-kind† {auto|tissue|cell-line}(DE を行うか)"]),
        ("⑦ 出力と可視化", ["metacell_obs.csv ほか(表3)、scanpy/(UMAP)", "metacellcnv_visualize.py で図とレポート"],
         ["--results-dir  --engine {plotly|matplotlib|both}",
          "--no-scanpy  --scanpy-dir(UMAP・QC 表は <out-dir>/scanpy に自動出力)"]),
    ],
    "en": [
        ("1. Input", ["Cell Ranger filtered_feature_bc_matrix", "+ GTF + species marker table"],
         ["--cellranger-dir  --sample-id  --gtf  --gtf-gene-id auto",
          "--chromosome-map  --mito-chromosome  --markers {dog|human|mouse|CSV}",
          "--preflight-only (check the inputs, then stop)"]),
        ("2. Preprocessing / QC", ["HVG/PCA, coarse clusters, cell types", "mtDNA check, MAD QC, doublets"],
         ["--n-hvg  --n-pcs  --leiden-resolution",
          "--marker-rule  --marker-min-margin  --marker-min-within-z",
          "--no-pctmt-filter  --expected-doublet-rate  --skip-doublet-rescue"]),
        ("3. Metacells", ["Group similar cells to reduce noise", "(own implementation, no SEACells package)"],
         ["--cells-per-metacell  --metacell-seed",
          "--seacell-assignments auto (reuse existing assignments)"]),
        ("4. CNV inference", ["Windowed running mean of the difference", "from the reference (the reference sets the meaning; Fig. 2)"],
         ["--cnv-reference {celltype|none|external†}",
          "--normal-celltype  --normal-clusters  --cnv-refine",
          "--cnv-method {infercnv|bins}  --cnv-external-reference†"]),
        ("5. Malignant call", ["Per CNV clone (cnv_leiden):", "malignant or normal"],
         ["--malignant-call {auto|reference|gap|median}",
          "--malignant-n-sd 3.0  --min-clone-gap 0.25"]),
        ("6. Annotation, lineage, DE", ["Metacell annotation (anno_label), lineage (Dollo), DE"],
         ["--lineage  --lineage-level  --lineage-events",
          "--no-metacell-annotation (anno_label is on by default; Sec. 6.4)",
          "--sample-kind† {auto|tissue|cell-line} (whether to run DE)"]),
        ("7. Output / plots", ["metacell_obs.csv etc. (Table 3), scanpy/ (UMAP)", "figures and report: metacellcnv_visualize.py"],
         ["--results-dir  --engine {plotly|matplotlib|both}",
          "--no-scanpy  --scanpy-dir (UMAP / QC tables go to <out-dir>/scanpy)"]),
    ],
}
T = {
    "ja": dict(
        t1="パイプラインの流れと、各段階で効くオプション",
        t2="サンプルに応じた「参照」の選び方と、得られる出力の意味",
        yes="はい", no="いいえ",
        q=[["Q1  サンプル内に、マーカーで同定できる", "正常細胞(免疫・内皮など)が十分にある?"],
           ["Q2  マーカーでは同定できないが、", "正常と分かっているクラスタがある?"],
           ["Q3  同じ実験・同じ条件で得た", "別サンプルを複数用意できる?"],
           ["Q4  (どれもない:腫瘍のみ・単一株の細胞株など)"]],
        a=[["A. 既定:--cnv-reference celltype",
            "正常参照 metacell の cnv_score で閾値を決める(reference 法)",
            "→ 絶対的な CNV(参照に対する増減)として読める"],
           ["B. --normal-clusters 3,7 / --normal-celltype Fibroblast",
            "指定したクラスタを参照にする",
            "→ A と同じ読み方だが「正常である」ことは利用者の責任"],
           ["C. 注釈付きの代替:別サンプルから参照を作る",
            "--cnv-reference external† --cnv-external-reference† <CSV>",
            "→ 精度は中程度。サンプル間のずれが混入するため注釈必須"],
           ["D. --cnv-reference none --sample-kind† cell-line",
            "基準は全 metacell の平均 → 均一なクローンの CNV は消える",
            "→ サブクローン間の相対差だけが見える。絶対的な CNV は出ない"]],
        note=["共通の注意",
              "・公開データだけを参照にする場合は C と同じ扱いだが、精度はさらに低い(図3・図4)",
              "・C・D の結果は「判定不明の領域」を含む。不明な領域を「CNV なし」と読んではいけない(§7)",
              "・A・B でも、参照が少ない(例:24 個)と閾値が不安定。参照の数と cnv_score の分布を必ず確認する(§6)",
              "・「正常」と判定された metacell は「参照と区別できない」の意味で、正常細胞であることの確認ではない",
              "† は公開リポジトリ(2026-09-17 時点)にまだ入っていない、作業版のオプション"]),
    "en": dict(
        t1="Pipeline flow and the options that matter at each stage",
        t2="Choosing a reference for your sample, and what the output then means",
        yes="yes", no="no",
        q=[["Q1  Does the sample contain enough normal cells", "that markers can identify (immune, endothelial, ...)?"],
           ["Q2  Are there clusters known to be normal", "even though markers cannot identify them?"],
           ["Q3  Can you prepare several other samples", "from the same experiment and conditions?"],
           ["Q4  None of these (tumor-only, a single cell line, ...)"]],
        a=[["A. Default: --cnv-reference celltype",
            "Threshold from cnv_score of the normal reference metacells (reference method)",
            "-> Read as absolute CNV (gain/loss relative to the reference)"],
           ["B. --normal-clusters 3,7 / --normal-celltype Fibroblast",
            "Use the clusters you name as the reference",
            "-> Read as in A, but 'normal' is your responsibility"],
           ["C. Annotated substitute: build a reference from other samples",
            "--cnv-reference external† --cnv-external-reference† <CSV>",
            "-> Medium accuracy; between-sample shifts leak in, so annotate"],
           ["D. --cnv-reference none --sample-kind† cell-line",
            "Baseline = mean of all metacells; uniform clonal CNV cancels out",
            "-> Only relative differences between subclones; no absolute CNV"]],
        note=["Common cautions",
              "- Using only public data as the reference is treated like C, with even lower accuracy (Figs. 3, 4)",
              "- C and D contain 'undetermined' regions. Never read undetermined as 'no CNV' (Sec. 7)",
              "- Even in A and B, a small reference (e.g. 24) makes the threshold unstable. Check the reference count and cnv_score distribution (Sec. 6)",
              "- A metacell called 'normal' means 'indistinguishable from the reference', not 'confirmed normal'",
              "† = working-version option, not yet in the public repository (as of 2026-09-17)"]),
}


def svg_open(w, h):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
            f'font-family="{FONT}" font-size="13">\n<rect width="{w}" height="{h}" fill="#ffffff"/>\n'
            '<defs><marker id="ar" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto">'
            '<path d="M0,0 L8,4 L0,8 z" fill="#444"/></marker></defs>\n')


def box(x, y, w, h, lines, fill="#f7f7f7", stroke="#444", bold_first=True, size=13):
    s = f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" fill="{fill}" stroke="{stroke}" stroke-width="1.2"/>\n'
    top = y + h / 2 - (len(lines) - 1) * 9 + 4
    for i, t in enumerate(lines):
        weight = ' font-weight="bold"' if (i == 0 and bold_first) else ''
        s += (f'<text x="{x + w/2}" y="{top + i*18}" text-anchor="middle" '
              f'font-size="{size}"{weight}>{esc(t)}</text>\n')
    return s


def arrow(x1, y1, x2, y2, label=None, lx=None, ly=None):
    s = (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#444" stroke-width="1.4" '
         'marker-end="url(#ar)"/>\n')
    if label:
        s += (f'<text x="{lx if lx is not None else (x1+x2)/2+6}" '
              f'y="{ly if ly is not None else (y1+y2)/2}" font-size="12" fill="#b22">{esc(label)}</text>\n')
    return s


def fig1(lang, out):
    stages = STAGES[lang]
    s = svg_open(1000, 40 + len(stages) * 100)
    s += f'<text x="20" y="26" font-size="15" font-weight="bold">{esc(T[lang]["t1"])}</text>\n'
    y = 44
    for i, (title, desc, opts) in enumerate(stages):
        s += box(20, y, 340, 76, [title] + desc, fill="#eef4fb")
        oy = y + 38 - (len(opts) - 1) * 9
        s += f'<line x1="360" y1="{y+38}" x2="388" y2="{y+38}" stroke="#999" stroke-dasharray="3,3"/>\n'
        for j, o in enumerate(opts):
            s += (f'<text x="394" y="{oy + j*18 + 4}" font-family="monospace" font-size="12" '
                  f'fill="#222">{esc(o)}</text>\n')
        if i < len(stages) - 1:
            s += arrow(190, y + 76, 190, y + 100)
        y += 100
    s += '</svg>\n'
    open(os.path.join(out, "fig1_pipeline_flow.svg"), "w", encoding="utf-8").write(s)


def fig2(lang, out):
    t = T[lang]
    s = svg_open(1000, 780)
    s += f'<text x="20" y="26" font-size="15" font-weight="bold">{esc(t["t2"])}</text>\n'
    qx, qw, qh, ex, ew = 20, 380, 74, 470, 510
    ys = [50, 190, 330, 470]
    fills = ["#e8f5e9", "#fff8e1", "#fff3e0", "#fdecea"]
    for i in range(4):
        yy = ys[i]
        s += box(qx, yy, qw, qh, t["q"][i], fill="#f7f7f7", bold_first=False)
        s += arrow(qx + qw, yy + qh / 2, ex, yy + qh / 2, t["yes"] if i < 3 else "", 408, yy + qh / 2 - 6)
        s += box(ex, yy, ew, qh, t["a"][i], fill=fills[i])
        if i < 3:
            s += arrow(qx + qw / 2, yy + qh, qx + qw / 2, ys[i + 1], t["no"],
                       qx + qw / 2 + 8, (yy + qh + ys[i + 1]) / 2 + 4)
    s += box(20, 580, 960, 170, t["note"], fill="#ffffff", bold_first=True, size=13)
    s += '</svg>\n'
    open(os.path.join(out, "fig2_reference_choice.svg"), "w", encoding="utf-8").write(s)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lang", choices=["ja", "en"], default="ja")
    a = ap.parse_args()
    out = os.path.join(HERE, "img") if a.lang == "ja" else os.path.join(HERE, "img", "en")
    os.makedirs(out, exist_ok=True)
    fig1(a.lang, out)
    fig2(a.lang, out)
    print("wrote figures to", out)
