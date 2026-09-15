#!/usr/bin/env python3
"""metacellcnv_plotly.py -- Plotly-based figure generation (all in-figure text is English).

See README.md for background, rationale, and usage.

Design goals (carried over from the matplotlib version):
- CNV values are polar (gain / neutral / loss): use a diverging scale
  (two colours + neutral grey), not a rainbow-type scale.
- Magnitude data such as cnv_score uses a single-hue sequential scale.
- Identity data relies on small multiples and direct labels, not colour alone.
- Avoid dual-axis plots with different scales on left and right.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

import plotly.graph_objs as go
from plotly.subplots import make_subplots

try:
    from scrna_common import log, warn, MITO_CHROMOSOMES
except Exception:  # pragma: no cover
    def log(m: str) -> None:
        print(f"[plot] {m}", flush=True)

    def warn(m: str) -> None:
        print(f"[plot][warn] {m}", flush=True)

    MITO_CHROMOSOMES = ("NC_002008.4", "chrM", "chrMT", "MT", "M")

#: Default chromosome-name variants dropped as mitochondrial by fig_cnv_heatmap
#: when the caller doesn't resolve --mito-chromosome itself (e.g. metacellcnv_
#: visualize.py's run_plotly() always does; this is the fallback for direct
#: API use). Mitochondrial DNA isn't diploid, so its CNV amplitude otherwise
#: dwarfs real nuclear CNV and masks it in the heatmap's color scale.
DEFAULT_MITO_EXCLUDE = {"chrMito"} | {str(c) for c in MITO_CHROMOSOMES}


__version__ = "1.0"

# Polar data (blue=loss / grey=neutral / red=gain). Colourblind-friendly blue-red.
DIVERGING = [[0.0, "#2166ac"], [0.5, "#f2f2f2"], [1.0, "#b2182b"]]
SEQUENTIAL = "Cividis"
# Discrete colours for identity data. Grey is reserved for "other / unclassified".
QUALITATIVE = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3",
               "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd"]
GRAY = "#b0b0b0"
LAYOUT = dict(plot_bgcolor="white", paper_bgcolor="white",
              font=dict(family="Helvetica, Arial, sans-serif", size=12),
              margin=dict(l=60, r=30, t=60, b=55))
AXIS = dict(linewidth=1, linecolor="black", mirror="ticks",
            showgrid=False, zeroline=False, ticks="outside")


# ---------------------------------------------------------------------------
# Label resolution
# ---------------------------------------------------------------------------

def resolve_group_labels(obs: pd.DataFrame,
                         prefer: Sequence[str] = ("anno_label", "cell_type")
                         ) -> tuple[pd.Series, str]:
    """Determine which column to use as the group label for figures.

    Since v3.0, metacell_obs.csv includes `anno_label` (assigned via a
    validated panel and shuffle calibration). The older `cell_type` column
    can retain unsupported names (e.g. `Myeloid_0`) that don't hold up under
    validation, so `anno_label` is always preferred when present to avoid
    mislabeling in figures.
    """
    for col in prefer:
        if col in obs.columns and obs[col].notna().any():
            s = obs[col].astype(str)
            if col != prefer[0]:
                warn(f"'{prefer[0]}' not found; using '{col}' as the group label"
                     " for figures. Note: legacy cell_type values may include"
                     " unsupported names")
            return s, col
    raise KeyError(f"No group label column found (looked for: {list(prefer)})")


def load_metacell_obs(results_dir: str | Path) -> pd.DataFrame:
    """Load metacell_obs.csv, merging in metacell_annotation.csv results if present.

    If Step 4a-4 was run afterward, or annotation was added to pre-v3.0
    output, the validated labels live in metacell_annotation.csv instead;
    merge them in here as `anno_label` so figures use them.
    """
    d = Path(results_dir)
    obs = pd.read_csv(d / "metacell_obs.csv", index_col=0)
    if "anno_label" not in obs.columns:
        for cand in (d / "metacell_annotation.csv",
                     d.parent / "anno_case1_real" / "metacell_annotation.csv"):
            if cand.exists():
                ann = pd.read_csv(cand, index_col=0)
                col = "label" if "label" in ann.columns else None
                if col:
                    obs["anno_label"] = ann[col].reindex(obs.index)
                    log(f"Merged annotation from: {cand.name}")
                break
    return obs


def prettify_labels(s: pd.Series) -> pd.Series:
    """Convert internal label strings to display-friendly English text."""
    mapping = {
        "Tumor:Unclassified": "Tumour (unclassified type)",
        "Unclassified:low-depth": "Unclassified (low depth)",
        "Unclassified:no-dominant-type": "Unclassified (no dominant type)",
        "Unclassified:multiple-types": "Unclassified (multiple types)",
        "Mixed:tumor+normal": "Mixed (tumour + normal)",
        "putative_malignant": "Putative malignant",
        "malignant": "Malignant",
        "normal": "Normal",
    }
    out = s.astype(str).map(lambda v: mapping.get(v, v.replace("Tumor:", "Tumour ")
                                                  .replace("_", " ")))
    return out


def _palette(levels: Iterable[str]) -> dict[str, str]:
    pal, i = {}, 0
    for lv in levels:
        if "nclassified" in lv or lv.lower() in ("other", "nan", "unknown"):
            pal[lv] = GRAY
        else:
            pal[lv] = QUALITATIVE[i % len(QUALITATIVE)]
            i += 1
    return pal


# ---------------------------------------------------------------------------
# 1. Scatter panels (from metacellcnv_scanpy output)
# ---------------------------------------------------------------------------

def load_prep_outputs(prep_dir: str | Path) -> dict:
    """Load metacellcnv_scanpy outputs (umap3d / clusters)."""
    d = Path(prep_dir)
    out = {}
    for key, name in (("umap", "umap3d.tsv.gz"), ("clusters", "clusters.tsv.gz"),
                      ("qc", "qc_metrics.tsv.gz")):
        p = d / name
        if p.exists():
            out[key] = pd.read_csv(p, sep="\t", index_col=0)
    if "umap" not in out:
        raise FileNotFoundError(f"{d}/umap3d.tsv.gz not found")
    return out


def fig_scatter_panels(prep_dir: str | Path, results_dir: str | Path, *,
                       three_d: bool = False, max_points: int = 60000,
                       seed: int = 0) -> go.Figure:
    """Plot Malignancy / Cell type / Cluster panels on the UMAP embedding.

    Joins `metacellcnv_scanpy`'s umap3d.tsv.gz and clusters.tsv.gz with the
    pipeline's cell_to_metacell.csv and metacell_obs.csv; each cell inherits
    its label from its metacell.
    """
    prep = load_prep_outputs(prep_dir)
    res = Path(results_dir)
    umap = prep["umap"]
    obs = load_metacell_obs(res)
    c2m = pd.read_csv(res / "cell_to_metacell.csv", index_col=0)
    mc_col = next((c for c in ("SEACell", "metacell") if c in c2m.columns),
                  c2m.columns[0])
    labels, used = resolve_group_labels(obs)
    labels = prettify_labels(labels)
    mal = (prettify_labels(obs["putative_malignant"])
           if "putative_malignant" in obs else None)

    cell_mc = c2m[mc_col].astype(str).reindex(umap.index)
    cell_lab = cell_mc.map(labels).fillna("Unassigned")
    cell_mal = cell_mc.map(mal).fillna("Unassigned") if mal is not None else None
    clus = (prep["clusters"].iloc[:, 0].astype(str).reindex(umap.index)
            if "clusters" in prep else None)

    idx = umap.index
    if len(idx) > max_points:
        idx = pd.Index(np.random.default_rng(seed)
                       .choice(idx, max_points, replace=False))
        log(f"Subsampling to {max_points:,} points for plotting (of {len(umap):,} total)")
    U = umap.loc[idx]
    x, y = U.iloc[:, 0].values, U.iloc[:, 1].values
    z = U.iloc[:, 2].values if U.shape[1] > 2 else None

    panels = [("Malignancy", cell_mal), (f"Cell type ({used})", cell_lab)]
    if clus is not None:
        panels.append(("Cluster", clus))
    panels = [(t, s) for t, s in panels if s is not None]

    kind = "scatter3d" if (three_d and z is not None) else "scattergl"
    fig = make_subplots(rows=1, cols=len(panels),
                        subplot_titles=[t for t, _ in panels],
                        specs=[[{"type": kind}] * len(panels)])
    for ci, (title, series) in enumerate(panels, start=1):
        vals = series.reindex(idx).astype(str)
        levels = sorted(vals.unique())
        pal = _palette(levels)
        for lv in levels:
            m = (vals == lv).values
            if not m.any():
                continue
            common = dict(mode="markers", name=lv, legendgroup=lv,
                          showlegend=(ci == 1 or len(panels) == 1),
                          text=[f"{b}<br>{lv}" for b in idx[m]],
                          hoverinfo="text",
                          marker=dict(size=3 if kind == "scattergl" else 2,
                                      opacity=0.55, color=pal[lv]))
            tr = (go.Scatter3d(x=x[m], y=y[m], z=z[m], **common)
                  if kind == "scatter3d"
                  else go.Scattergl(x=x[m], y=y[m], **common))
            fig.add_trace(tr, row=1, col=ci)
    fig.update_layout(title="Cell embedding coloured by annotation",
                      height=430, width=420 * len(panels), **LAYOUT)
    if kind == "scattergl":
        fig.update_xaxes(title_text="UMAP 1", **AXIS)
        fig.update_yaxes(title_text="UMAP 2", **AXIS)
    return fig


# ---------------------------------------------------------------------------
# 2. Chromosome-level CNV heatmap
# ---------------------------------------------------------------------------

def fig_cnv_heatmap(mat: pd.DataFrame, obs: pd.DataFrame | None = None, *,
                    title: str = "Chromosome-level CNV",
                    mito_exclude: set[str] | None = None) -> go.Figure:
    """CNV per metacell x chromosome; rows ordered by hierarchical clustering.

    mito_exclude: chromosome-name variants (case-insensitive) to drop as
    mitochondrial before rendering. Defaults to DEFAULT_MITO_EXCLUDE. Applied
    here (not just upstream) so this figure is safe from mito columns
    regardless of which caller built `mat` -- a pre-computed
    cnv_chromosome_means.csv, or a table derived on the fly.
    """
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

    # Keep only numeric columns (cnv_chromosome_means.csv can include annotation columns)
    M = mat.apply(pd.to_numeric, errors="coerce")
    drop = [c for c in M.columns if M[c].isna().all()]
    if drop:
        log(f"Dropped non-numeric columns: {drop}")
        M = M.drop(columns=drop)
    mito_lower = {str(m).lower() for m in (mito_exclude if mito_exclude is not None else DEFAULT_MITO_EXCLUDE)}
    drop_mito = [c for c in M.columns if str(c).lower() in mito_lower]
    if drop_mito:
        log(f"Dropped mitochondrial chromosome column(s): {drop_mito}")
        M = M.drop(columns=drop_mito)
    M = M.fillna(0.0)
    if M.shape[1] == 0:
        raise ValueError("No numeric chromosome columns found")
    if len(M) > 2:
        order = leaves_list(linkage(pdist(M.values), method="ward",
                                    optimal_ordering=len(M) <= 1500))
        M = M.iloc[order]
    lim = float(np.nanpercentile(np.abs(M.values), 99)) or 1.0
    hover = None
    if obs is not None:
        lab, _ = resolve_group_labels(obs)
        lab = prettify_labels(lab).reindex(M.index)
        hover = [[f"{mc}<br>{lab.get(mc, '')}<br>{c}: {v:+.3f}"
                  for c, v in zip(M.columns, row)] for mc, row in
                 zip(M.index, M.values)]
    fig = go.Figure(go.Heatmap(
        z=M.values, x=list(M.columns), y=list(M.index),
        colorscale=DIVERGING, zmid=0.0, zmin=-lim, zmax=lim,
        text=hover, hoverinfo="text" if hover else None,
        colorbar=dict(title="log2 ratio")))
    fig.update_layout(title=title, height=max(420, min(900, 2 + 2 * len(M))),
                      **LAYOUT)
    fig.update_xaxes(title_text="Chromosome", **AXIS)
    fig.update_yaxes(title_text="Metacell (hierarchically ordered)",
                     showticklabels=len(M) <= 60, **AXIS)
    return fig


# ---------------------------------------------------------------------------
# 3. DEG swarm plot
# ---------------------------------------------------------------------------

def _beeswarm_offsets(y: np.ndarray, width: float = 0.34,
                      n_bins: int = 24) -> np.ndarray:
    """Spread points at the same y-value horizontally (no seaborn dependency)."""
    if y.size == 0:
        return y
    order = np.argsort(y)
    bins = np.linspace(np.nanmin(y), np.nanmax(y) + 1e-12, n_bins + 1)
    which = np.clip(np.digitize(y, bins) - 1, 0, n_bins - 1)
    off = np.zeros_like(y, dtype=float)
    for b in range(n_bins):
        idx = order[which[order] == b]
        k = idx.size
        if k <= 1:
            continue
        pos = (np.arange(k) - (k - 1) / 2.0) / max((k - 1) / 2.0, 1.0)
        off[idx] = pos * width
    return off


def fig_deg_swarm(expr: pd.DataFrame, genes: Sequence[str],
                  obs: pd.DataFrame, *, max_cols: int = 5,
                  stats: pd.DataFrame | None = None) -> go.Figure:
    """Plot per-metacell expression of top DEGs, grouped by label.

    Group labels come from `resolve_group_labels`, which prefers
    `anno_label` over the legacy `cell_type` column, so unsupported names
    like `Myeloid_0` won't appear directly in the figure.

    `stats`, if given, is the DE table (indexed by gene, with
    `log2FoldChange` and `padj` columns) used to annotate each subplot title.
    """
    labels, used = resolve_group_labels(obs)
    labels = prettify_labels(labels).reindex(expr.index)
    keep = labels.notna()
    expr, labels = expr.loc[keep], labels[keep]
    genes = [g for g in genes if g in expr.columns]
    if not genes:
        raise ValueError("No genes available to plot")
    levels = sorted(labels.unique())
    pal = _palette(levels)

    def _title(g: str) -> str:
        if stats is None or g not in stats.index:
            return g
        lfc = float(stats.loc[g, "log2FoldChange"])
        padj = float(stats.loc[g, "padj"])
        ptxt = "padj<1e-300" if padj < 1e-300 else f"padj={padj:.1e}"
        return f"{g}  log2FC={lfc:+.2f}  {ptxt}"

    ncol = min(max_cols, len(genes))
    nrow = int(np.ceil(len(genes) / ncol))
    fig = make_subplots(rows=nrow, cols=ncol,
                        subplot_titles=[_title(g) for g in genes],
                        shared_yaxes=False, vertical_spacing=0.12,
                        horizontal_spacing=0.06)
    for gi, g in enumerate(genes):
        r, c = gi // ncol + 1, gi % ncol + 1
        for li, lv in enumerate(levels):
            m = (labels == lv).values
            yv = expr.loc[m, g].values.astype(float)
            if yv.size == 0:
                continue
            xv = li + _beeswarm_offsets(yv)
            fig.add_trace(go.Scattergl(
                x=xv, y=yv, mode="markers", name=lv, legendgroup=lv,
                showlegend=(gi == 0),
                text=[f"{mc}<br>{lv}<br>{g} = {v:.3f}"
                      for mc, v in zip(expr.index[m], yv)],
                hoverinfo="text",
                marker=dict(size=5, opacity=0.7, color=pal[lv],
                            line=dict(width=0))), row=r, col=c)
        fig.update_xaxes(tickvals=list(range(len(levels))),
                         ticktext=levels, tickangle=-35, row=r, col=c, **AXIS)
        fig.update_yaxes(title_text="log-normalised expression" if c == 1 else None,
                         row=r, col=c, **AXIS)
    fig.update_layout(title=f"Top differential genes by {used}",
                      height=300 * nrow, width=290 * ncol, **LAYOUT)
    return fig


# ---------------------------------------------------------------------------
# 4. Lineage tree (Dollo)
# ---------------------------------------------------------------------------

def fig_lineage_tree(root, *, has_structure: bool = True,
                     title: str = "Tumour lineage (Dollo parsimony)") -> go.Figure:
    """Draw the cnv_lineage.Node tree as a standalone lineage figure (not attached to a heatmap).

    x is the cumulative event count from the root (regions lost along that
    lineage); y is leaf order. When the structure isn't statistically
    supported, a warning is drawn directly into the figure so it isn't
    mistaken for a validated genealogy.
    """
    leaves: list[tuple[str, object]] = []

    def collect(n):
        if not n.children and not n.members:
            leaves.append((n.name, n))
        for m in n.members:
            leaves.append((m, n))
        for c in n.children:
            collect(c)

    collect(root)
    ypos: dict[int, float] = {}
    for i, (_, node) in enumerate(leaves):
        ypos.setdefault(id(node), [])
    ys: dict[int, list[float]] = {}
    for i, (_, node) in enumerate(leaves):
        ys.setdefault(id(node), []).append(float(i))

    xs: dict[int, float] = {}

    def depth(n, d=0.0):
        xs[id(n)] = d
        for c in n.children:
            depth(c, d + max(len(c.events), 1))

    depth(root)

    def ycoord(n) -> float:
        if id(n) in ys and ys[id(n)]:
            vals = list(ys[id(n)])
        else:
            vals = []
        for c in n.children:
            vals.append(ycoord(c))
        return float(np.mean(vals)) if vals else 0.0

    yc: dict[int, float] = {}

    def fill(n):
        for c in n.children:
            fill(c)
        yc[id(n)] = ycoord(n)

    fill(root)

    ex, ey = [], []

    def edges(n):
        # += inside a nested function would rebind ex as a local variable, so use extend instead
        for c in n.children:
            ex.extend([xs[id(n)], xs[id(n)], xs[id(c)], None])
            ey.extend([yc[id(n)], yc[id(c)], yc[id(c)], None])
            edges(c)

    edges(root)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=ex, y=ey, mode="lines", hoverinfo="skip",
                             line=dict(color="#444", width=1.4),
                             showlegend=False))
    nx_, ny_, ntext = [], [], []
    for n in _walk(root):
        nx_.append(xs[id(n)]); ny_.append(yc[id(n)])
        ntext.append(f"{n.name}<br>events: {', '.join(n.events) or '(root)'}"
                     f"<br>metacells here: {len(n.members)}")
    fig.add_trace(go.Scatter(x=nx_, y=ny_, mode="markers",
                             marker=dict(size=8, color="#2166ac"),
                             text=ntext, hoverinfo="text",
                             name="Node", showlegend=False))
    fig.update_layout(title=title + ("" if has_structure else
                                     "  —  NOT SUPPORTED BY THE PERMUTATION TEST"),
                      height=max(420, 14 * max(len(leaves), 10)),
                      **LAYOUT)
    fig.update_xaxes(title_text="Cumulative losses from diploid root", **AXIS)
    fig.update_yaxes(title_text="", showticklabels=False, **AXIS)
    if not has_structure:
        fig.add_annotation(
            xref="paper", yref="paper", x=0.5, y=1.0, showarrow=False,
            text="<b>Lineage structure is not supported</b> — "
                 "do not read this tree as a genealogy",
            font=dict(color="#b2182b", size=13))
    return fig


def _walk(n):
    yield n
    for c in n.children:
        yield from _walk(c)


# ---------------------------------------------------------------------------
# 5. QC
# ---------------------------------------------------------------------------

def fig_qc(qc: pd.DataFrame) -> go.Figure:
    """Per-cell QC (depth, genes detected, mitochondrial fraction)."""
    cols = [c for c in ("total_counts", "n_genes_by_counts", "pct_counts_mt")
            if c in qc.columns]
    if not cols:
        raise ValueError("No QC metric columns found")
    names = {"total_counts": "Total counts per cell",
             "n_genes_by_counts": "Genes detected per cell",
             "pct_counts_mt": "Mitochondrial fraction (%)"}
    fig = make_subplots(rows=1, cols=len(cols),
                        subplot_titles=[names[c] for c in cols])
    for i, c in enumerate(cols, start=1):
        v = qc[c].astype(float).dropna().values
        fig.add_trace(go.Histogram(x=v, nbinsx=60, marker_color="#4c72b0",
                                   showlegend=False), row=1, col=i)
        fig.add_vline(x=float(np.median(v)), line=dict(color="#b2182b", width=1.5,
                                                       dash="dash"), row=1, col=i)
        fig.update_xaxes(row=1, col=i, **AXIS)
        fig.update_yaxes(title_text="Cells" if i == 1 else None, row=1, col=i, **AXIS)
    fig.update_layout(title="Per-cell quality control (dashed line = median)",
                      height=360, width=380 * len(cols), **LAYOUT)
    return fig


def fig_metacell_quality(obs: pd.DataFrame) -> go.Figure:
    """Metacell size, compactness, and separation."""
    panels = [(c, t) for c, t in (("n_cells", "Cells per metacell"),
                                  ("compactness", "Compactness (lower is better)"),
                                  ("separation", "Separation (higher is better)"),
                                  ("cnv_score", "CNV score"))
              if c in obs.columns]
    if not panels:
        raise ValueError("No metacell quality metric columns found")
    fig = make_subplots(rows=1, cols=len(panels),
                        subplot_titles=[t for _, t in panels])
    for i, (c, _) in enumerate(panels, start=1):
        v = obs[c].astype(float).dropna()
        fig.add_trace(go.Box(y=v.values, name="", boxpoints="all", jitter=0.4,
                             pointpos=0, marker=dict(size=3, opacity=0.5,
                                                     color="#4c72b0"),
                             line=dict(color="#2166ac"), showlegend=False,
                             text=list(v.index), hoverinfo="text+y"),
                      row=1, col=i)
        fig.update_xaxes(showticklabels=False, row=1, col=i, **AXIS)
        fig.update_yaxes(row=1, col=i, **AXIS)
    fig.update_layout(title="Metacell quality", height=380,
                      width=300 * len(panels), **LAYOUT)
    return fig


# ---------------------------------------------------------------------------
# 6. Combine into a single HTML report
# ---------------------------------------------------------------------------

def write_report(figures: dict[str, go.Figure], path: str | Path, *,
                 title: str = "metacellcnv report",
                 png_dir: str | Path | None = None,
                 pdf_dir: str | Path | None = None) -> Path:
    """Combine figures into a single HTML report; also writes PNGs/PDFs if those dirs are given.

    Static export (png_dir, pdf_dir) uses go.Figure.write_image (Kaleido) and is
    best-effort: a missing/broken Kaleido install skips it with a warning rather
    than failing the report, since the HTML report is already complete without it.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    parts = [
        "<!doctype html><meta charset='utf-8'>",
        f"<title>{title}</title>",
        "<style>body{font-family:Helvetica,Arial,sans-serif;margin:24px;"
        "background:#fff;color:#222}h1{font-size:20px}h2{font-size:15px;"
        "margin-top:28px;border-bottom:1px solid #ddd;padding-bottom:4px}</style>",
        f"<h1>{title}</h1>",
    ]
    first = True
    for name, fig in figures.items():
        parts.append(f"<h2>{name}</h2>")
        parts.append(fig.to_html(full_html=False,
                                 include_plotlyjs="cdn" if first else False))
        first = False
    out.write_text("\n".join(parts), encoding="utf-8")
    log(f"HTML report: {out} ({len(figures)} figures)")

    for ext, target_dir in (("png", png_dir), ("pdf", pdf_dir)):
        if not target_dir:
            continue
        d = Path(target_dir); d.mkdir(parents=True, exist_ok=True)
        kwargs = {"scale": 2} if ext == "png" else {}
        for i, (name, fig) in enumerate(figures.items(), start=1):
            safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in name)
            try:
                fig.write_image(str(d / f"{i:02d}_{safe}.{ext}"), **kwargs)
            except Exception as exc:
                warn(f"Skipping {ext.upper()} export ({name}): {str(exc).strip()[:80]}."
                     " To enable static image export, run `pip install kaleido`"
                     " (and, on some systems, `plotly_get_chrome`)."
                     f" The HTML report is complete without {ext.upper()}s.")
                break
        else:
            log(f"{ext.upper()}: {d}")
    return out
