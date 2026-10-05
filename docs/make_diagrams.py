#!/usr/bin/env python3
"""ガイドの図 1(パイプライン概観)と図 2(参照の選び方)を SVG で出力する。"""
import os
from xml.sax.saxutils import escape as esc
HERE = os.path.dirname(os.path.abspath(__file__))
FONT = "'Noto Sans CJK JP','Hiragino Sans','Yu Gothic','Meiryo',sans-serif"

def svg_open(w, h):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
            f'font-family="{FONT}" font-size="13">\n<rect width="{w}" height="{h}" fill="#ffffff"/>\n'
            '<defs><marker id="ar" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto">'
            '<path d="M0,0 L8,4 L0,8 z" fill="#444"/></marker></defs>\n')

def box(x, y, w, h, lines, fill="#f7f7f7", stroke="#444", bold_first=True, size=13):
    s = f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" fill="{fill}" stroke="{stroke}" stroke-width="1.2"/>\n'
    n = len(lines)
    top = y + h / 2 - (n - 1) * 9 + 4
    for i, t in enumerate(lines):
        weight = ' font-weight="bold"' if (i == 0 and bold_first) else ''
        s += f'<text x="{x + w/2}" y="{top + i*18}" text-anchor="middle" font-size="{size}"{weight}>{esc(t)}</text>\n'
    return s

def arrow(x1, y1, x2, y2, label=None, lx=None, ly=None):
    s = f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#444" stroke-width="1.4" marker-end="url(#ar)"/>\n'
    if label:
        s += f'<text x="{lx if lx is not None else (x1+x2)/2+6}" y="{ly if ly is not None else (y1+y2)/2}" font-size="12" fill="#b22">{esc(label)}</text>\n'
    return s

# ---------------- 図 1 ----------------
stages = [
    ("① 入力", ["Cell Ranger の filtered_feature_bc_matrix", "+ GTF + 種別のマーカー表"],
     ["--cellranger-dir  --sample-id  --gtf  --gtf-gene-id auto", "--chromosome-map  --mito-chromosome  --markers {dog|human|mouse|CSV}", "--preflight-only(入力だけ検査して終了)"]),
    ("② 前処理・QC", ["HVG/PCA・粗クラスタ・細胞型注釈", "mtDNA 構成検査・MAD QC・doublet"],
     ["--n-hvg  --n-pcs  --leiden-resolution", "--marker-rule  --marker-min-margin  --marker-min-within-z", "--no-pctmt-filter  --expected-doublet-rate  --skip-doublet-rescue"]),
    ("③ metacell 構築", ["近い細胞を束ねてノイズを減らす", "(SEACells 非依存の自前実装)"],
     ["--cells-per-metacell  --metacell-seed", "--seacell-assignments auto(既存の割り当てを再利用)"]),
    ("④ CNV 推定", ["参照との差を窓移動平均で CNV にする", "参照の選び方が結果の意味を決める(図2)"],
     ["--cnv-reference {celltype|none|external†}", "--normal-celltype  --normal-clusters  --cnv-refine", "--cnv-method {infercnv|bins}  --cnv-external-reference†"]),
    ("⑤ 悪性の判定", ["CNV クローン(cnv_leiden)単位で", "悪性 / 正常を決める"],
     ["--malignant-call {auto|reference|gap|median}", "--malignant-n-sd 3.0  --min-clone-gap 0.25"]),
    ("⑥ 任意の解析", ["系譜(Dollo)・metacell 注釈・DE"],
     ["--lineage  --lineage-level  --lineage-events", "--metacell-annotation", "--sample-kind† {auto|tissue|cell-line}(DE を行うか)"]),
    ("⑦ 出力と可視化", ["metacell_obs.csv ほか(表3)", "metacellcnv_visualize.py で図とレポート"],
     ["--results-dir  --prep-dir  --engine {plotly|matplotlib|both}"]),
]
H = 40 + len(stages) * 100
s = svg_open(1000, H)
s += '<text x="20" y="26" font-size="15" font-weight="bold">パイプラインの流れと、各段階で効くオプション</text>\n'
y = 44
for i, (t, desc, opts) in enumerate(stages):
    s += box(20, y, 340, 76, [t] + desc, fill="#eef4fb")
    oy = y + 38 - (len(opts) - 1) * 9
    s += f'<line x1="360" y1="{y+38}" x2="388" y2="{y+38}" stroke="#999" stroke-dasharray="3,3"/>\n'
    for j, o in enumerate(opts):
        s += f'<text x="394" y="{oy + j*18 + 4}" font-family="monospace" font-size="12" fill="#222">{esc(o)}</text>\n'
    if i < len(stages) - 1:
        s += arrow(190, y + 76, 190, y + 100)
    y += 100
s += '</svg>\n'
open(os.path.join(HERE, "img", "fig1_pipeline_flow.svg"), "w", encoding="utf-8").write(s)

# ---------------- 図 2 ----------------
s = svg_open(1000, 780)
s += '<text x="20" y="26" font-size="15" font-weight="bold">サンプルに応じた「参照」の選び方と、得られる出力の意味</text>\n'
# 質問(左列)と結論(右列)
qx, qw, qh = 20, 380, 74
ex, ew = 470, 510
ys = [50, 190, 330, 470]
qs = [
    ["Q1  サンプル内に、マーカーで同定できる", "正常細胞(免疫・内皮など)が十分にある?"],
    ["Q2  マーカーでは同定できないが、", "正常と分かっているクラスタがある?"],
    ["Q3  同じ実験・同じ条件で得た", "別サンプルを複数用意できる?"],
    ["Q4  (どれもない:腫瘍のみ・単一株の細胞株など)"],
]
ends = [
    (["A. 既定:--cnv-reference celltype", "正常参照 metacell の cnv_score で閾値を決める(reference 法)", "→ 絶対的な CNV(参照に対する増減)として読める"], "#e8f5e9"),
    (["B. --normal-clusters 3,7 / --normal-celltype Fibroblast", "指定したクラスタを参照にする", "→ A と同じ読み方だが「正常である」ことは利用者の責任"], "#fff8e1"),
    (["C. 注釈付きの代替:別サンプルから参照を作る", "--cnv-reference external† --cnv-external-reference† <CSV>", "→ 精度は中程度。サンプル間のずれが混入するため注釈必須"], "#fff3e0"),
    (["D. --cnv-reference none --sample-kind† cell-line", "基準は全 metacell の平均 → 均一なクローンの CNV は消える", "→ サブクローン間の相対差だけが見える。絶対的な CNV は出ない"], "#fdecea"),
]
for i, (q, (e, col)) in enumerate(zip(qs, ends)):
    yy = ys[i]
    s += box(qx, yy, qw, qh, q, fill="#f7f7f7", bold_first=False)
    s += arrow(qx + qw, yy + qh / 2, ex, yy + qh / 2, "はい" if i < 3 else "", 408, yy + qh / 2 - 6)
    s += box(ex, yy, ew, qh, e, fill=col)
    if i < 3:
        s += arrow(qx + qw / 2, yy + qh, qx + qw / 2, ys[i + 1], "いいえ", qx + qw / 2 + 8, (yy + qh + ys[i + 1]) / 2 + 4)
s += box(20, 580, 960, 170, [
    "共通の注意",
    "・公開データだけを参照にする場合は C と同じ扱いだが、精度はさらに低い(図3・図4)",
    "・C・D の結果は「判定不明の領域」を含む。不明な領域を「CNV なし」と読んではいけない(§7)",
    "・A・B でも、参照が少ない(例:24 個)と閾値が不安定。参照の数と cnv_score の分布を必ず確認する(§6)",
    "・「正常」と判定された metacell は「参照と区別できない」の意味で、正常細胞であることの確認ではない",
    "† は公開リポジトリ(2026-09-17 時点)にまだ入っていない、作業版のオプション"],
    fill="#ffffff", bold_first=True, size=13)
s += '</svg>\n'
open(os.path.join(HERE, "img", "fig2_reference_choice.svg"), "w", encoding="utf-8").write(s)
print("done")
