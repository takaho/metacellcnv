#!/usr/bin/env python3
"""pipeline_guide.md から、白背景のシンプルな index.html を生成する(GitHub Pages 用)。
実行: python docs/build_html.py   (リポジトリのルートで)
"""
import os
import re
import markdown

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "pipeline_guide.md"), encoding="utf-8").read()

md = markdown.Markdown(extensions=["tables", "fenced_code", "toc", "sane_lists"],
                       extension_configs={"toc": {"toc_depth": "2-3", "permalink": False}})
body = md.convert(src)

# 表をスクロール可能なラッパで包む(狭い画面で横にはみ出さない)
# 図は原寸の画像へのリンクにする(クリックで拡大)
body = re.sub(r'(<img src="([^"]+)"[^>]*>)', r'<a href="\2">\1</a>', body)
body = re.sub(r"(<table>.*?</table>)", r'<div class="tablewrap">\1</div>', body, flags=re.S)

def toc_html(tokens):
    out = ["<ul>"]
    for t in tokens:
        if t["level"] == 1:
            continue
        out.append(f'<li class="l{t["level"]}"><a href="#{t["id"]}">{t["name"]}</a></li>')
    out.append("</ul>")
    return "\n".join(out)

flat = []
def walk(ts):
    for t in ts:
        flat.append(t)
        walk(t.get("children", []))
walk(md.toc_tokens)
toc = toc_html(flat)

CSS = """
:root { color-scheme: light; }
html { background: #ffffff; }
body { margin: 0; background: #ffffff; color: #222222;
  font-family: -apple-system, BlinkMacSystemFont, "Hiragino Sans", "Yu Gothic", "Noto Sans CJK JP", "Segoe UI", sans-serif;
  line-height: 1.75; font-size: 16px; }
.page { max-width: 900px; margin: 0 auto; padding: 24px 16px 64px; }
h1 { font-size: 1.7rem; margin: 0.4em 0 0.6em; }
h2 { font-size: 1.3rem; margin: 2.2em 0 0.6em; padding-bottom: 0.25em; border-bottom: 1px solid #cccccc; }
h3 { font-size: 1.1rem; margin: 1.6em 0 0.4em; }
a { color: #0b5cad; }
nav.toc { border: 1px solid #dddddd; padding: 8px 16px; margin: 16px 0 28px; background: #fafafa; font-size: 0.92rem; }
nav.toc ul { list-style: none; margin: 0; padding: 0; columns: 2; }
nav.toc li { margin: 2px 0; break-inside: avoid; }
nav.toc li.l3 { padding-left: 1.2em; font-size: 0.9em; }
.tablewrap { overflow-x: auto; margin: 12px 0 20px; }
table { border-collapse: collapse; font-size: 0.9rem; min-width: 100%; }
th, td { border: 1px solid #cccccc; padding: 6px 10px; vertical-align: top; text-align: left; }
th { background: #f3f3f3; }
code { background: #f4f4f4; padding: 1px 4px; border-radius: 3px; font-size: 0.88em; }
pre { background: #f4f4f4; padding: 10px 12px; overflow-x: auto; border-radius: 4px; }
pre code { background: none; padding: 0; }
figure { margin: 20px 0 28px; padding: 0; }
figure img { max-width: 100%; height: auto; display: block; margin: 0 auto; border: 1px solid #eeeeee; background: #ffffff; }
figcaption { font-size: 0.86rem; color: #444444; margin-top: 8px; line-height: 1.6; }
blockquote { border-left: 3px solid #cccccc; margin: 12px 0; padding: 2px 14px; color: #444444; }
footer { margin-top: 48px; padding-top: 12px; border-top: 1px solid #dddddd; font-size: 0.82rem; color: #666666; }
@media (max-width: 640px) { nav.toc ul { columns: 1; } body { font-size: 15px; } }
"""

NAV = '</h1>\n<nav class="toc">\n' + toc + '\n</nav>'
BODY = body.replace("</h1>", NAV, 1)
html = f"""<!doctype html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>metacellcnv 実行ガイド</title>
<style>{CSS}</style>
</head>
<body>
<div class="page">
{BODY}
<footer>metacellcnv 実行ガイド。Markdown 版: <a href="pipeline_guide.md">pipeline_guide.md</a>(同じ内容)。<code>docs/build_html.py</code> で生成。</footer>
</div>
</body>
</html>
"""
open(os.path.join(HERE, "index.html"), "w", encoding="utf-8").write(html)
print("wrote index.html", len(html))
