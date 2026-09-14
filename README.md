# metacellcnv

**Marker-free identification of malignant cells in tumor scRNA-seq, by copy-number variation (CNV) over metacells.**

Most single-cell tumor pipelines call malignant cells with marker genes, but marker expression is not stable across tumor origins, and outside human tissue (dog, mouse, ...) markers are often only an analogy borrowed from human biology. `metacellcnv` instead aggregates cells into metacells, estimates copy-number variation per metacell directly from expression counts, and calls malignancy from chromosomal-level CNV evidence — the same logic pathologists use (aneuploidy), applied per metacell rather than per cell to make the signal statistically tractable.

The pipeline takes CellRanger output (`filtered_feature_bc_matrix/`) and a gene-annotation GTF as input, and needs no R/Seurat, no SEACells, and no infercnvpy: metacell construction and CNV inference are both native reimplementations (see `metacells_native.py`, `cnv_native.py`), validated against those reference tools during development.

## What it does

1. **Preprocess** (`metacellcnv_scanpy.py`) — cell/gene QC, normalization, HVGs, PCA, UMAP, clustering. Scanpy-only; usable standalone.
2. **Metacell construction** (`metacells_native.py`) — archetypal-analysis metacells (same objective as SEACells), with a factored-kernel implementation that avoids materializing an n×n cell-cell matrix.
3. **CNV inference** (`cnv_native.py`) — per-metacell copy-number signal from binned gene expression, using an **empirical (per-bin) null** rather than a normal approximation — see "On z-scores" below.
4. **Malignancy calling & lineage** (`metacellcnv.py`, `metacell_annotation.py`, `cnv_lineage.py`) — CNV-based malignant/normal calls, cell-type annotation from a marker table with an explicit evidence trail (typed / measured-mixture / unclassified — never a bare guess), and a Dollo-parsimony tumor lineage tree (chromosome loss is treated as irreversible) with output in Newick format.
5. **Differential expression** (`metacellcnv.py`, via pyDESeq2) — malignant vs. normal metacells, with sample as a covariate.
6. **Visualization** (`metacellcnv_visualize.py` / `metacellcnv_plotly.py`) — an interactive, self-contained HTML report (Plotly): embeddings, QC, CNV heatmaps, DEG swarm plots, the lineage tree.

Run the whole thing with `python metacellcnv.py`, or invoke the preprocessing/visualization stages standalone with `python metacellcnv.py --scanpy ...` / `python metacellcnv.py --visualize ...` (equivalently, run `metacellcnv_scanpy.py` / `metacellcnv_visualize.py` directly).

## On z-scores

CNV bin statistics are **not** normally distributed — cell-to-cell variation in expression programs (not sampling noise) makes the tails heavier than a normal approximation predicts, so a fixed z-score cutoff runs at several times its nominal false-positive rate. `cnv_native.py` uses an empirical, per-bin null (built by shuffling gene-to-bin assignment) instead of a parametric normal/binomial/beta-binomial quantile — see the module docstring in `cnv_native.py` for the short version, and `smoke_test.py` (`########## v3.1: cnv_native ##########`) for the calibration check.

## Species support

Cell-type marker panels are external CSV tables (`markers/markers_{dog,human,mouse}.csv`), not hardcoded — pass `--markers <species-or-csv-path>` to `metacellcnv.py`. Dog panels are validated against real depth-matched data; the human and mouse panels were produced by orthology conversion from the dog panel and are **not yet independently validated** — the pipeline prints a warning at runtime whenever an unvalidated panel is used for labeling. See `markers/build_marker_tables.py` to regenerate or extend the tables, and the CSV schema (`cell_type,gene,role,use_as_normal_reference,use_for_labels,validation`) to add a new species or tissue panel.

## Files

| File | Role |
|---|---|
| `scrna_common.py` | Shared building blocks: GTF parsing, chromosome-name handling, loading/normalization, mtDNA identification, QC, doublet detection. No CNV/DE dependency. |
| `metacellcnv_scanpy.py` | Unified preprocessing: cell selection → normalization → HVG → PCA → UMAP (3D) → clustering → mtx/TSV export. Depends only on `scrna_common.py`. |
| `metacellcnv.py` | Main pipeline: QC → metacells → CNV → malignancy calls → DE. Also dispatches `--scanpy`/`--visualize` to the two scripts below. |
| `metacellcnv_visualize.py` / `metacellcnv_plotly.py` | Turns pipeline output into an interactive HTML report. |
| `metacell_annotation.py` | Classifies metacells as typed / measured-mixture / unclassified, from validated marker panels + karyotype projection + shuffle-calibrated thresholds — never a silent default label. |
| `metacells_native.py` | Metacell construction (archetypal analysis). No SEACells dependency. Factored kernel, exact RSS tracking, size-degeneracy detection. |
| `cnv_native.py` | CNV estimation. No infercnvpy dependency. Pyramid-window running mean (matches infercnvpy's definition) plus a bin-composition beta-binomial/negative-binomial test with an empirical null. |
| `cnv_lineage.py` | Dollo-parsimony tumor lineage tree from CNV loss events, with a permutation test (`independent` and degree-preserving `curveball` nulls) for whether the data actually supports tree structure at all. |
| `subtype_within_compartment.py` | Splits metacells into normal/malignant compartments by CNV, then classifies subtypes within each — with tests for whether CNV itself is confounding the subtype call. |
| `export_metacell_membership.py` | Writes out which barcodes belong to which metacell. |
| `blas_check.py` | Standalone diagnostic: BLAS thread count and effective dense-matmul throughput. |
| `smoke_test.py` | Synthetic-data regression test suite. Run this after any change. |
| `markers/` | `markers_{dog,human,mouse}.csv` marker tables and `build_marker_tables.py` to regenerate them. |
| `mito_gene_profile_case1_reference.csv` | Reference mtDNA gene-composition profile, used by the mitochondrial-composition QC check. |

### Module dependency graph

```
scrna_common.py                    (no CNV/DE dependency)
   ├── metacellcnv_scanpy.py       preprocessing only; usable standalone
   ├── metacells_native.py         metacell construction (no SEACells)
   ├── cnv_native.py               CNV inference (no infercnvpy)
   ├── metacell_annotation.py      metacell typing
   ├── cnv_lineage.py              tumor lineage tree
   └── metacellcnv.py              metacells -> CNV -> malignancy -> DE
```

`metacellcnv_scanpy.py` does not import `metacellcnv.py` (or pydeseq2), so it can be used for preprocessing alone, independent of the CNV/DE analysis. `smoke_test.py` checks this independence, and that the pipeline never imports `SEACells` / `infercnvpy` / `palantir`, by walking the AST of the pipeline module rather than trusting a comment.

## Installation and quick start

See [INSTALL.md](INSTALL.md).

```bash
# preprocessing
python metacellcnv.py --scanpy --cellranger-dir <sample>/outs/filtered_feature_bc_matrix \
    --sample-id <sample> --gtf genes.gtf --mito-chromosome <mito_contig> --out-dir prep/<sample>

# main pipeline
python metacellcnv.py --cellranger-dir <sample>/outs/filtered_feature_bc_matrix \
    --gtf genes.gtf --sample-id <sample> --markers dog --out-dir results/<sample>

# interactive HTML report
python metacellcnv.py --visualize --results-dir results/<sample> --prep-dir prep/<sample>
```

Run `python metacellcnv.py --help` for the full option list.

## Testing

```bash
python smoke_test.py
```

builds small synthetic CellRanger-format datasets and exercises every pipeline stage end to end (metacell construction, CNV inference and its empirical-null calibration, malignancy calling, Dollo lineage reconstruction, the plotly report, and the marker-table loader). It should print `ALL SMOKE TESTS PASSED` and exit 0 with no further errors; run it after any change before trusting the result.

## Status

This is a research pipeline, not a validated clinical tool. In particular: the empirical evaluation to date rests on a single tumor specimen, so treat performance numbers as a proof of concept rather than a generalization claim; the human and mouse marker panels are orthology conversions, not independently validated (the pipeline warns about this at runtime); and CNV inference detects chromosome/arm-level events well but has limited power for small focal events. A manuscript describing the method and its validation is in preparation.

## License

No license file is included yet — until one is added, all rights are reserved by the author. Contact the author before reuse.
