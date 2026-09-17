#!/usr/bin/env python3
"""
Export the cell-barcode -> metacell membership table from singlecells_qc.h5ad
===============================================================================

Pipeline v1.9+ writes results/cell_to_metacell.csv automatically at run time,
but older outputs only carry this mapping in singlecells_qc.h5ad's obs. Since
the h5ad can be hundreds of MB, this reads obs only and never touches the
expression matrix.

Usage:
    python export_metacell_membership.py results_c1
    python export_metacell_membership.py results_c1 --metacell MC-3 MC-4

Column names accept both the new "metacell" name and the legacy "SEACell"
name (from an h5ad written before the rename).

See README.md for details.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from scrna_common import resolve_metacell_key


def read_obs_only(h5ad_path: Path) -> pd.DataFrame:
    """Read only obs via h5py (skips X, so even a multi-hundred-MB file is fast)."""
    import h5py

    with h5py.File(h5ad_path, "r") as f:
        obs = f["obs"]
        index_key = obs.attrs["_index"]
        barcodes = [
            b.decode() if isinstance(b, bytes) else str(b) for b in obs[index_key][()]
        ]
        data: dict[str, object] = {}
        for key in obs.keys():
            if key == index_key:
                continue
            node = obs[key]
            if isinstance(node, h5py.Group):  # categorical
                if "categories" not in node or "codes" not in node:
                    continue
                cats = [
                    c.decode() if isinstance(c, bytes) else c
                    for c in node["categories"][()]
                ]
                codes = node["codes"][()]
                data[key] = [cats[c] if c >= 0 else None for c in codes]
            else:
                arr = node[()]
                if arr.dtype.kind == "S":
                    arr = np.array([x.decode() for x in arr])
                data[key] = arr
    return pd.DataFrame(data, index=pd.Index(barcodes, name="barcode"))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Export the cell-barcode -> metacell membership table")
    ap.add_argument("results_dir", help="Pipeline output directory")
    ap.add_argument(
        "--h5ad", default="singlecells_qc.h5ad", help="Single-cell h5ad file name"
    )
    ap.add_argument("--out", default="cell_to_metacell.csv", help="Output CSV name")
    ap.add_argument(
        "--metacell",
        nargs="*",
        default=None,
        help="Print the member barcodes of the given metacell(s) to stdout",
    )
    args = ap.parse_args(argv)

    results = Path(args.results_dir)
    h5ad = results / args.h5ad
    if not h5ad.exists():
        print(f"Not found: {h5ad}", file=sys.stderr)
        return 1

    obs = read_obs_only(h5ad)
    mc_key = resolve_metacell_key(obs.columns)
    if mc_key not in obs:
        print(
            "obs has no metacell column ('metacell' / legacy 'SEACell'). Pass an"
            " h5ad that went through build_metacells().",
            file=sys.stderr,
        )
        return 1

    keep = [
        c
        for c in (mc_key, "sample_id", "coarse_cluster", "cell_type", "total_counts",
                  "n_genes_by_counts", "pct_counts_mt", "doublet_score")
        if c in obs
    ]
    table = obs[keep]
    out_path = results / args.out
    table.to_csv(out_path)
    print(f"Saved: {out_path} ({len(table):,} cells x {len(keep)} columns)")

    sizes = table[mc_key].value_counts()
    print(
        f"metacells: {len(sizes)} / member cells: median {int(sizes.median())}"
        f" range {int(sizes.min())}-{int(sizes.max())}"
    )

    if args.metacell:
        for mc in args.metacell:
            bcs = table.index[table[mc_key].astype(str) == str(mc)].tolist()
            print(f"\n=== {mc}: {len(bcs)} cells ===")
            for b in bcs:
                print(b)
    return 0


if __name__ == "__main__":
    sys.exit(main())
