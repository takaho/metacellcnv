#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""metacellcnv_de.py -- malignant vs normal DE, overall and within each subtype.

Why this module exists
======================
The DE table of the main run compares ALL malignant metacells with ALL normal
metacells. Malignant and normal metacells differ in cell type as well (tumor
epithelium vs immune cells), so that table mostly reflects subtype
differences. To describe the tumor itself, the same comparison is also made
INSIDE each subtype: malignant vs normal metacells of the same type.

What it does
============
1. Subtype = `anno_label` with the "Tumor:" prefix removed. Labels that name no
   single type (Unclassified:*, Mixed:*) are not subtypes.
2. Comparisons = "Overall" + every subtype that has enough data in BOTH groups
   (at least --min-metacells metacells and --min-cells cells; the cells are the
   sum of `n_cells` of the metacells). The others are listed with the reason.
3. Each comparison is a pyDESeq2 test at the metacell level (same function and
   settings as the main run's DE): malignant vs normal, covariate sample_id.
4. Outputs (under <out-dir>):
     de_malignant_vs_normal_by_subtype.xlsx   one sheet per comparison + Summary
     de_by_subtype/<name>.csv                 the same tables as CSV
     figures/de/de_<name>.html (+ index.html) volcano plot + top-gene expression
   (figures are drawn by `--visualize`, or here unless --no-figures)

Reading the results
===================
* Metacells of one sample are pseudo-replicates, so p-values are not
  independent evidence; use them to rank genes, not as proof.
* The malignant / normal groups come from the CNV call, and the subtype from
  the annotation, so the groups are defined partly from the same data
  (circularity). A subtype whose "normal" group is small or low in CNV
  signal can be unreliable.
* padj is corrected within each comparison, not across comparisons.

Usage
-----
    python metacellcnv_de.py --results-dir results/<sample>
    python metacellcnv.py --de --results-dir results/<sample>      # same thing
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scrna_common import log, warn
except Exception:  # pragma: no cover
    def log(m: str) -> None:
        print(f"[de] {m}", flush=True)

    def warn(m: str) -> None:
        print(f"[de][warn] {m}", flush=True)

__version__ = "1.0"

GROUP_KEY = "putative_malignant"
TEST, REF = "malignant", "normal"
# Thresholds used to call a gene "differential" (same as the DEG swarm plot)
PADJ_MAX, LFC_MIN, BASEMEAN_MIN = 0.05, 0.5, 5.0
XLSX_NAME = "de_malignant_vs_normal_by_subtype.xlsx"
CSV_DIR = "de_by_subtype"


# ---------------------------------------------------------------------------
# Subtypes and comparisons
# ---------------------------------------------------------------------------

def subtype_labels(obs: pd.DataFrame, key: str = "anno_label") -> tuple[pd.Series, str]:
    """Subtype of each metacell (NaN when no single type is named).

    Uses `key` (default anno_label). If it is missing, falls back to
    `cell_type` with a warning (its names are weaker evidence).
    """
    col = key
    if col not in obs.columns or obs[col].isna().all():
        col = "cell_type"
        warn(f"'{key}' not found; using 'cell_type' as the subtype. These names "
             "have weaker evidence (re-run with the metacell annotation to get anno_label)")
    s = obs[col].astype(str).str.replace(r"^Tumor:", "", regex=True)
    bad = s.str.startswith(("Unclassified", "Mixed")) | s.isin(["nan", "None", ""])
    return s.mask(bad), col


def plan_comparisons(obs: pd.DataFrame, sub: pd.Series, *, min_metacells: int = 3,
                     min_cells: int = 100) -> pd.DataFrame:
    """Table with one row per candidate comparison, with `status` and `reason`."""
    n_cells = obs["n_cells"].astype(float) if "n_cells" in obs.columns else pd.Series(1.0, index=obs.index)
    g = obs[GROUP_KEY].astype(str)
    rows = []

    def one(name: str, mask: pd.Series) -> dict:
        r = {"comparison": name}
        for grp in (TEST, REF):
            m = mask & (g == grp)
            r[f"n_metacells_{grp}"] = int(m.sum())
            r[f"n_cells_{grp}"] = int(n_cells[m].sum())
        why = []
        for grp in (TEST, REF):
            if r[f"n_metacells_{grp}"] < min_metacells:
                why.append(f"{grp}: {r[f'n_metacells_{grp}']} metacells (< {min_metacells})")
            elif r[f"n_cells_{grp}"] < min_cells:
                why.append(f"{grp}: {r[f'n_cells_{grp}']} cells (< {min_cells})")
        r["status"] = "run" if not why else "skipped"
        r["reason"] = "; ".join(why)
        return r

    rows.append(one("Overall", pd.Series(True, index=obs.index)))
    order = (n_cells.groupby(sub).sum().sort_values(ascending=False).index
             if sub.notna().any() else [])
    for st in order:
        rows.append(one(str(st), sub == st))
    return pd.DataFrame(rows)


def safe_name(name: str, maxlen: int = 31) -> str:
    """File- and Excel-sheet-safe name."""
    s = re.sub(r"[\[\]:*?/\\]", "_", name)
    s = re.sub(r"[^A-Za-z0-9._ +-]", "_", s).strip()
    return s[:maxlen] or "x"


# ---------------------------------------------------------------------------
# Computation
# ---------------------------------------------------------------------------

def _default_de_fn():
    import importlib
    return importlib.import_module("metacellcnv").run_de_branch


def compute_tables(adata, obs: pd.DataFrame, plan: pd.DataFrame, sub: pd.Series, *,
                   de_fn=None, overall: pd.DataFrame | None = None,
                   covariate_key: str | None = "sample_id", min_total_counts: int = 10,
                   n_cpus: int = 4) -> dict[str, pd.DataFrame]:
    """Run the DE of every comparison marked `run`. Returns {comparison: DE table}."""
    de_fn = de_fn or _default_de_fn()
    tables: dict[str, pd.DataFrame] = {}
    for _, row in plan[plan.status == "run"].iterrows():
        name = row["comparison"]
        if name == "Overall" and overall is not None:
            tables[name] = overall.copy()
            log(f"DE [{name}]: reusing the main run's table")
            continue
        mask = (pd.Series(True, index=obs.index) if name == "Overall" else (sub == name))
        mask &= obs[GROUP_KEY].astype(str).isin([TEST, REF])
        sel = adata[mask.reindex(adata.obs_names).fillna(False).values].copy()
        sel.obs[GROUP_KEY] = obs[GROUP_KEY].astype(str).reindex(sel.obs_names).values
        log(f"DE [{name}]: {int((sel.obs[GROUP_KEY] == TEST).sum())} malignant vs "
            f"{int((sel.obs[GROUP_KEY] == REF).sum())} normal metacells")
        try:
            _, res = de_fn(sel, condition_key=GROUP_KEY, condition_test=TEST,
                           condition_ref=REF, covariate_key=covariate_key,
                           min_total_counts=min_total_counts, n_cpus=n_cpus)
            tables[name] = res
        except Exception as exc:
            warn(f"DE [{name}] failed: {exc}")
            plan.loc[plan.comparison == name, ["status", "reason"]] = ["failed", str(exc)[:200]]
    return tables


def _log_cpm(adata) -> pd.DataFrame:
    """log1p(CPM) of the raw metacell counts (genes as columns)."""
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    X = X.tocsr() if hasattr(X, "tocsr") else np.asarray(X)
    tot = np.asarray(X.sum(axis=1)).ravel().astype(float)
    tot[tot == 0] = 1.0
    return X, tot


def _expr_for(adata, genes: list[str]) -> pd.DataFrame:
    X, tot = _log_cpm(adata)
    idx = [adata.var_names.get_loc(g) for g in genes if g in adata.var_names]
    cols = [g for g in genes if g in adata.var_names]
    sub = X[:, idx]
    sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
    return pd.DataFrame(np.log1p(sub / tot[:, None] * 1e6), index=adata.obs_names, columns=cols)


def annotate_table(res: pd.DataFrame, adata, obs: pd.DataFrame, mask: pd.Series) -> pd.DataFrame:
    """Add direction and per-group mean expression columns to a DE table."""
    out = res.copy()
    ok = out["padj"].notna() & (out["padj"] < PADJ_MAX) & (out["log2FoldChange"].abs() >= LFC_MIN) \
        & (out["baseMean"] >= BASEMEAN_MIN)
    out.insert(0, "direction", np.where(ok & (out.log2FoldChange > 0), "up in malignant",
                                        np.where(ok & (out.log2FoldChange < 0), "down in malignant", "ns")))
    g = obs[GROUP_KEY].astype(str).reindex(adata.obs_names)
    m = mask.reindex(adata.obs_names).fillna(False).values
    X, tot = _log_cpm(adata)
    for grp, col in ((TEST, "mean_logCPM_malignant"), (REF, "mean_logCPM_normal")):
        rows = np.where(m & (g == grp).values)[0]
        if rows.size == 0:
            out[col] = np.nan
            continue
        sub = X[rows]
        sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
        v = pd.Series(np.log1p(sub / tot[rows][:, None] * 1e6).mean(axis=0), index=adata.var_names)
        out[col] = v.reindex(out.index).values
    out.index.name = "gene"
    return out.sort_values(["padj", "pvalue"], na_position="last")


# ---------------------------------------------------------------------------
# Excel / CSV
# ---------------------------------------------------------------------------

def write_outputs(out_dir: Path, tables: dict[str, pd.DataFrame], plan: pd.DataFrame,
                  params: dict) -> Path:
    """Write the CSVs and the Excel workbook (Summary + one sheet per comparison)."""
    out_dir = Path(out_dir)
    csv_dir = out_dir / CSV_DIR
    csv_dir.mkdir(parents=True, exist_ok=True)
    summ = plan.copy()
    n_up, n_down = [], []
    for name in summ.comparison:
        t = tables.get(name)
        n_up.append(int((t.direction == "up in malignant").sum()) if t is not None else np.nan)
        n_down.append(int((t.direction == "down in malignant").sum()) if t is not None else np.nan)
        if t is not None:
            t.to_csv(csv_dir / f"{safe_name(name, 80)}.csv")
    summ["n_up_in_malignant"], summ["n_down_in_malignant"] = n_up, n_down
    summ.insert(1, "sheet", [safe_name(n) if n in tables else "" for n in summ.comparison])

    notes = pd.DataFrame({"item": [
        "Comparison", "Subtype", "Test", "Differential gene", "Multiple testing",
        "Caution 1", "Caution 2", "Parameters", "Columns"], "description": [
        "malignant vs normal metacells (log2FoldChange > 0 = higher in malignant)",
        "anno_label with the 'Tumor:' prefix removed; Unclassified:* and Mixed:* are not subtypes",
        "pyDESeq2 on raw metacell counts; covariate sample_id (dropped if one level)",
        f"padj < {PADJ_MAX}, |log2FoldChange| >= {LFC_MIN}, baseMean >= {BASEMEAN_MIN}",
        "padj is corrected within each comparison, not across comparisons",
        "Metacells of one sample are pseudo-replicates: use the p-values to rank genes, not as proof",
        "Malignant/normal comes from the CNV call and the subtype from the annotation, so the groups "
        "are partly defined by the same data (circularity)",
        "; ".join(f"{k}={v}" for k, v in params.items()),
        "direction, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj (pyDESeq2); "
        "mean_logCPM_* = mean log1p(CPM) of the metacells in each group"]})

    path = out_dir / XLSX_NAME
    tmp = out_dir / (XLSX_NAME + ".tmp")   # write aside, then rename: no half-written workbook
    with pd.ExcelWriter(tmp, engine="openpyxl") as xw:
        summ.to_excel(xw, sheet_name="Summary", index=False)
        notes.to_excel(xw, sheet_name="Notes", index=False)
        used = {"Summary", "Notes"}
        for name, t in tables.items():
            sn, k = safe_name(name), 2
            while sn in used:
                sn = f"{safe_name(name, 28)}_{k}"
                k += 1
            used.add(sn)
            t.reset_index().to_excel(xw, sheet_name=sn, index=False)
        _format_workbook(xw.book)
    tmp.replace(path)
    log(f"Wrote: {path} ({len(tables)} comparison sheet(s) + Summary, Notes)")
    return path


def _format_workbook(wb) -> None:
    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter
    fill = PatternFill("solid", fgColor="EEEEEE")
    for ws in wb.worksheets:
        ws.freeze_panes = "A2"
        for c in ws[1]:
            c.font = Font(bold=True)
            c.fill = fill
        for i, col in enumerate(ws.columns, start=1):
            head = str(col[0].value or "")
            width = max([len(head)] + [len(str(c.value)) for c in col[1:60] if c.value is not None])
            ws.column_dimensions[get_column_letter(i)].width = min(max(10, width + 2), 70 if ws.title == "Notes" else 40)
            if head in ("pvalue", "padj"):
                for c in col[1:]:
                    c.number_format = "0.00E+00"
            elif head in ("baseMean", "log2FoldChange", "lfcSE", "stat", "mean_logCPM_malignant",
                          "mean_logCPM_normal"):
                for c in col[1:]:
                    c.number_format = "0.000"
        if ws.title == "Notes":
            from openpyxl.styles import Alignment
            for row in ws.iter_rows(min_row=2):
                for c in row:
                    c.alignment = Alignment(wrap_text=True, vertical="top")
        elif ws.title != "Summary" and ws.max_row > 1:
            ws.auto_filter.ref = ws.dimensions


# ---------------------------------------------------------------------------
# Figures (plotly)
# ---------------------------------------------------------------------------

def pick_genes(t: pd.DataFrame, top_n: int) -> list[str]:
    """Top-N differential genes, half up / half down (the rest filled from the other side)."""
    sig = t[t.direction != "ns"]
    up = sig[sig.direction == "up in malignant"].sort_values(["padj", "pvalue"])
    dn = sig[sig.direction == "down in malignant"].sort_values(["padj", "pvalue"])
    k = top_n // 2
    take_up, take_dn = list(up.index[:k]), list(dn.index[:top_n - k])
    short = top_n - len(take_up) - len(take_dn)
    if short > 0:
        take_up += [g for g in up.index[k:] if g not in take_up][:short]
        take_dn += [g for g in dn.index[top_n - k:] if g not in take_dn][:max(0, short - (len(take_up) - k))]
    return (take_up + take_dn)[:top_n]


def fig_volcano(t: pd.DataFrame, title: str, n_label: int = 12):
    import plotly.graph_objs as go
    import metacellcnv_plotly as PL
    d = t.dropna(subset=["padj", "log2FoldChange"]).copy()
    d["nlp"] = -np.log10(d.padj.clip(lower=1e-300))
    fig = go.Figure()
    for lab, col in (("ns", "#b0b0b0"), ("down in malignant", "#2166ac"), ("up in malignant", "#b2182b")):
        s = d[d.direction == lab]
        if s.empty:
            continue
        fig.add_trace(go.Scattergl(
            x=s.log2FoldChange, y=s.nlp, mode="markers", name=f"{lab} ({len(s)})",
            text=[f"{g}<br>log2FC={a:+.2f}<br>padj={b:.1e}" for g, a, b in zip(s.index, s.log2FoldChange, s.padj)],
            hoverinfo="text", marker=dict(size=4, opacity=0.7, color=col)))
    lab = d[d.direction != "ns"].sort_values("padj").head(n_label)
    if len(lab):
        fig.add_trace(go.Scatter(x=lab.log2FoldChange, y=lab.nlp, mode="text", text=list(lab.index),
                                 textposition="top center", showlegend=False, hoverinfo="skip",
                                 textfont=dict(size=10)))
    fig.add_hline(y=-np.log10(PADJ_MAX), line=dict(dash="dash", color="#888", width=1))
    for x in (-LFC_MIN, LFC_MIN):
        fig.add_vline(x=x, line=dict(dash="dash", color="#888", width=1))
    fig.update_layout(title=title, height=520, width=760, **PL.LAYOUT)
    fig.update_xaxes(title_text="log2 fold change (malignant / normal)", **PL.AXIS)
    fig.update_yaxes(title_text="-log10 padj", **PL.AXIS)
    return fig


def fig_top_genes(expr: pd.DataFrame, groups: pd.Series, genes: list[str], t: pd.DataFrame, title: str,
                  max_cols: int = 4):
    import plotly.graph_objs as go
    from plotly.subplots import make_subplots
    import metacellcnv_plotly as PL
    genes = [g for g in genes if g in expr.columns]
    if not genes:
        raise ValueError("no genes to plot")
    ncol = min(max_cols, len(genes))
    nrow = int(np.ceil(len(genes) / ncol))

    def ttl(g):
        return f"<b>{g}</b><br>log2FC={t.loc[g, 'log2FoldChange']:+.2f}  padj={t.loc[g, 'padj']:.1e}"

    fig = make_subplots(rows=nrow, cols=ncol, subplot_titles=[ttl(g) for g in genes],
                        vertical_spacing=0.16, horizontal_spacing=0.07)
    for a in fig.layout.annotations:
        a.font.size = 11
    cols = {REF: "#4c72b0", TEST: "#c44e52"}
    for gi, g in enumerate(genes):
        r, c = gi // ncol + 1, gi % ncol + 1
        for li, lv in enumerate((REF, TEST)):
            m = (groups == lv).values
            y = expr.loc[m, g].values.astype(float)
            if y.size == 0:
                continue
            fig.add_trace(go.Scattergl(
                x=li + PL._beeswarm_offsets(y), y=y, mode="markers", name=lv, legendgroup=lv,
                showlegend=(gi == 0),
                text=[f"{mc}<br>{lv}<br>{g} = {v:.2f}" for mc, v in zip(expr.index[m], y)],
                hoverinfo="text", marker=dict(size=5, opacity=0.7, color=cols[lv])), row=r, col=c)
            fig.add_trace(go.Scatter(x=[li - 0.3, li + 0.3], y=[np.median(y)] * 2, mode="lines",
                                     line=dict(color="black", width=2), showlegend=False,
                                     hoverinfo="skip"), row=r, col=c)
        fig.update_xaxes(tickvals=[0, 1], ticktext=["normal", "malignant"], range=[-0.6, 1.6],
                         row=r, col=c, **PL.AXIS)
        fig.update_yaxes(title_text="log1p CPM" if c == 1 else None, row=r, col=c, **PL.AXIS)
    fig.update_layout(title=title, height=300 * nrow + 90, width=250 * ncol + 100, **{**PL.LAYOUT, 'margin': dict(l=60, r=30, t=120, b=55)})
    return fig


def draw_figures(fig_root: Path, adata, obs: pd.DataFrame, sub: pd.Series, plan: pd.DataFrame,
                 tables: dict[str, pd.DataFrame], *, top_n: int = 12) -> list[Path]:
    """One HTML per comparison (volcano + top-gene expression) and an index page, in <fig_root>/de."""
    import metacellcnv_plotly as PL
    fig_dir = Path(fig_root) / "de"
    fig_dir.mkdir(parents=True, exist_ok=True)
    made = []
    for _, row in plan.iterrows():
        name = row["comparison"]
        t = tables.get(name)
        if t is None:
            continue
        mask = (pd.Series(True, index=obs.index) if name == "Overall" else (sub == name))
        mask &= obs[GROUP_KEY].astype(str).isin([TEST, REF])
        keep = mask.reindex(adata.obs_names).fillna(False).values
        a = adata[keep]
        groups = obs[GROUP_KEY].astype(str).reindex(a.obs_names)
        head = (f"{name}: malignant ({row.n_metacells_malignant} metacells, {row.n_cells_malignant} cells) vs "
                f"normal ({row.n_metacells_normal} metacells, {row.n_cells_normal} cells)")
        figs = {}
        try:
            figs["Volcano"] = fig_volcano(t, head)
            genes = pick_genes(t, top_n)
            if genes:
                figs[f"Top {len(genes)} genes (log1p CPM per metacell)"] = fig_top_genes(
                    _expr_for(a, genes), groups, genes, t, head)
            else:
                warn(f"DE [{name}]: no gene passes padj<{PADJ_MAX}, |log2FC|>={LFC_MIN}; volcano only")
        except Exception as exc:
            warn(f"Figure for [{name}] failed: {exc}")
            continue
        p = fig_dir / f"de_{safe_name(name, 80)}.html"
        PL.write_report(figs, p, title=f"DE {name} (malignant vs normal)", png_dir=None,
                        pdf_dir=None)
        made.append(p)
    idx = ["<!doctype html><meta charset='utf-8'><title>DE by subtype</title>",
           "<style>body{font-family:Helvetica,Arial,sans-serif;margin:24px;background:#fff;color:#222}"
           "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:4px 10px;text-align:left}"
           "th{background:#f3f3f3}</style>", "<h1>Malignant vs normal: overall and within each subtype</h1>",
           f"<p>Statistics (one sheet per comparison): <code>{XLSX_NAME}</code> in the results directory</p>",
           "<table><tr><th>comparison</th><th>malignant (metacells / cells)</th><th>normal (metacells / cells)</th>"
           "<th>up</th><th>down</th><th>status</th></tr>"]
    for _, r in plan.iterrows():
        t = tables.get(r.comparison)
        link = (f"<a href='de_{safe_name(r.comparison, 80)}.html'>{r.comparison}</a>"
                if t is not None else r.comparison)
        up = int((t.direction == "up in malignant").sum()) if t is not None else ""
        dn = int((t.direction == "down in malignant").sum()) if t is not None else ""
        idx.append(f"<tr><td>{link}</td><td>{r.n_metacells_malignant} / {r.n_cells_malignant}</td>"
                   f"<td>{r.n_metacells_normal} / {r.n_cells_normal}</td><td>{up}</td><td>{dn}</td>"
                   f"<td>{r.status}{(': ' + r.reason) if r.reason else ''}</td></tr>")
    idx.append("</table>")
    (fig_dir / "index.html").write_text("\n".join(idx), encoding="utf-8")
    log(f"DE figures: {fig_dir}/index.html ({len(made)} comparison(s))")
    return made


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def load_results(res_dir: Path):
    """metacells.h5ad (raw counts) and metacell_obs.csv (+ annotation if separate)."""
    import anndata as ad
    res_dir = Path(res_dir)
    adata = ad.read_h5ad(res_dir / "metacells.h5ad")
    obs = pd.read_csv(res_dir / "metacell_obs.csv", index_col=0)
    if "anno_label" not in obs.columns:
        ann = res_dir / "metacell_annotation.csv"
        if ann.exists():
            a = pd.read_csv(ann, index_col=0)
            if "label" in a.columns:
                obs["anno_label"] = a["label"].reindex(obs.index)
                log(f"Merged annotation from {ann.name}")
    return adata, obs


def run(res_dir, *, out_dir=None, fig_dir=None, adata=None, obs=None, de_fn=None, overall=None,
        group_key: str = "anno_label", min_metacells: int = 3, min_cells: int = 100,
        top_n: int = 12, n_cpus: int = 4, min_total_counts: int = 10, figures: bool = True,
        recompute: bool = False) -> dict:
    """Compute (or load) the tables, write Excel/CSV, and draw figures."""
    res_dir = Path(res_dir)
    out_dir = Path(out_dir) if out_dir else res_dir
    if adata is None or obs is None:
        adata, obs = load_results(res_dir)
    for c in (GROUP_KEY, "n_cells"):
        if c not in obs.columns and c in adata.obs.columns:
            obs[c] = adata.obs[c].reindex(obs.index)
    if GROUP_KEY not in obs.columns:
        raise KeyError(f"{GROUP_KEY} is not in metacell_obs.csv")
    if "anno_label" not in obs.columns and "anno_label" in adata.obs.columns:
        obs["anno_label"] = adata.obs["anno_label"].reindex(obs.index)
    sub, used = subtype_labels(obs, group_key)
    plan = plan_comparisons(obs, sub, min_metacells=min_metacells, min_cells=min_cells)
    log("DE comparisons (subtype from '%s'):\n%s" % (used, plan.drop(columns="reason").to_string(index=False)))

    xlsx = out_dir / XLSX_NAME
    tables: dict[str, pd.DataFrame] = {}
    if xlsx.exists() and not recompute and adata is not None:
        try:
            sheets = pd.read_excel(xlsx, sheet_name=None, index_col=None)
            for name in plan[plan.status == "run"].comparison:
                sn = safe_name(name)
                if sn in sheets:
                    tables[name] = sheets[sn].set_index("gene")
            log(f"Using {len(tables)} existing table(s) from {xlsx.name} (use --recompute to redo)")
        except Exception as exc:
            warn(f"Could not read {xlsx.name} ({exc}); recomputing")
    todo = plan[(plan.status == "run") & (~plan.comparison.isin(tables))]
    if len(todo):
        raw = compute_tables(adata, obs, plan, sub, de_fn=de_fn, overall=overall,
                             n_cpus=n_cpus, min_total_counts=min_total_counts)
        for name, t in raw.items():
            mask = (pd.Series(True, index=obs.index) if name == "Overall" else (sub == name))
            tables[name] = annotate_table(t, adata, obs, mask)
        params = dict(subtype_key=used, min_metacells=min_metacells, min_cells=min_cells,
                      padj_max=PADJ_MAX, lfc_min=LFC_MIN, version=__version__)
        write_outputs(out_dir, tables, plan, params)
    if figures and tables:
        draw_figures(Path(fig_dir) if fig_dir else out_dir / "figures", adata, obs, sub, plan, tables,
                     top_n=top_n)
    return {"plan": plan, "tables": tables}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Malignant vs normal DE, overall and within each subtype "
                                f"(metacellcnv_de v{__version__})")
    p.add_argument("--results-dir", required=True, help="Output directory of metacellcnv.py")
    p.add_argument("--out-dir", default=None, help="Where to write the tables (default: --results-dir)")
    p.add_argument("--fig-dir", default=None, help="Figure directory (default: <out-dir>/figures; files go in its de/)")
    p.add_argument("--group-key", default="anno_label",
                   help="Column that names the subtype (default anno_label; the 'Tumor:' prefix is removed)")
    p.add_argument("--de-min-metacells", type=int, default=3,
                   help="Minimum metacells in EACH of malignant and normal for a subtype comparison (default 3)")
    p.add_argument("--de-min-cells", type=int, default=100,
                   help="Minimum cells in EACH group for a subtype comparison (default 100)")
    p.add_argument("--de-top-n", type=int, default=12, help="Genes shown in the expression panel (default 12)")
    p.add_argument("--n-cpus", type=int, default=4)
    p.add_argument("--min-total-counts", type=int, default=10)
    p.add_argument("--no-figures", action="store_true", help="Tables only")
    p.add_argument("--recompute", action="store_true", help="Redo the DE even if the workbook exists")
    a = p.parse_args(argv)
    run(a.results_dir, out_dir=a.out_dir, fig_dir=a.fig_dir, group_key=a.group_key, min_metacells=a.de_min_metacells,
        min_cells=a.de_min_cells, top_n=a.de_top_n, n_cpus=a.n_cpus,
        min_total_counts=a.min_total_counts, figures=not a.no_figures, recompute=a.recompute)
    return 0


if __name__ == "__main__":
    sys.exit(main())
