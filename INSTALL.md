# Installation

## Requirements

- Python >= 3.10
- A gene-annotation GTF matching your reference genome (for mapping genes to genomic coordinates)
- CellRanger output: one `filtered_feature_bc_matrix/` directory per sample

## Core dependencies

Required for preprocessing, metacell construction, and CNV inference:

```bash
pip install numpy pandas scipy scikit-learn anndata scanpy
```

Required for the full pipeline (malignancy calling already covered above; these two are needed for differential expression and doublet detection):

```bash
pip install pydeseq2 scrublet
```

Required for the default (Plotly) visualization report:

```bash
pip install plotly
```

## Optional dependencies

| Package | Needed for |
|---|---|
| `matplotlib` | Legacy static-figure report (`metacellcnv.py --visualize --engine matplotlib`) |
| `infercnvpy` | One figure in the legacy matplotlib report only (`cnv.pl.chromosome_heatmap`); **not** needed by the CNV pipeline itself, which does not depend on infercnvpy at all |
| `kaleido` | Exporting the Plotly report's figures to static PNG in addition to the HTML report. Kaleido needs a Chrome/Chromium install; if it's missing, the HTML report is still written in full, just without PNGs |
| `threadpoolctl` | Only for the standalone `blas_check.py` diagnostic |

Nothing else in this repository imports `SEACells` or `infercnvpy` for its own computation — metacell construction and CNV inference are both native reimplementations. (`decoupler` is checked for at pipeline startup and reported if present, but is not actually required by anything here — safe to ignore.)

## Suggested setup

```bash
conda create -n metacellcnv python=3.11 -y
conda activate metacellcnv
pip install numpy pandas scipy scikit-learn anndata scanpy pydeseq2 scrublet plotly
```

Add `matplotlib` and/or `infercnvpy` only if you specifically want the legacy matplotlib report path.

## Verifying the install

```bash
python smoke_test.py
```

This builds small synthetic CellRanger-format datasets on the fly (no external data needed) and runs every pipeline stage against them. A successful run prints `ALL SMOKE TESTS PASSED` partway through and continues to completion with no errors — if any dependency is missing or misconfigured, this is where it will surface.

## Running on your own data

```bash
# 1. Preprocess (per sample)
python metacellcnv.py --scanpy \
    --cellranger-dir <sample>/outs/filtered_feature_bc_matrix \
    --sample-id <sample> \
    --gtf genes.gtf --gtf-gene-id auto \
    --mito-chromosome <mito_contig_or_MT> \
    --out-dir prep/<sample>

# 2. Run the pipeline (metacells -> CNV -> malignancy -> DE)
python metacellcnv.py \
    --cellranger-dir <sample>/outs/filtered_feature_bc_matrix \
    --gtf genes.gtf \
    --sample-id <sample> \
    --markers dog \
    --out-dir results/<sample>

# 3. Build the interactive report
python metacellcnv.py --visualize \
    --results-dir results/<sample> --prep-dir prep/<sample>
```

`--markers` accepts `dog`, `human`, or `mouse` (built-in tables in `markers/`), or a path to your own CSV following the same schema (`cell_type,gene,role,use_as_normal_reference,use_for_labels,validation`) — see `markers/build_marker_tables.py`.

Run `python metacellcnv.py --help`, `python metacellcnv.py --scanpy --help`, and `python metacellcnv.py --visualize --help` for the full option lists, including CNV window/step size, chromosome exclusion, lineage-reconstruction parameters, and report engine selection.
