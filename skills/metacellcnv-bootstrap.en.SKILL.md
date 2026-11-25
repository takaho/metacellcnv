---
name: metacellcnv-bootstrap
description: >
  Given only a CellRanger output directory, fetch the metacellcnv pipeline
  (https://github.com/takaho/metacellcnv) from GitHub, build a matching
  isolated Python environment with all dependencies, auto-resolve the gene
  annotation (GTF) and chromosome map from a genome name, and run
  preprocessing -> metacell construction -> CNV inference -> malignancy
  calling -> DE -> the interactive HTML report, end to end. The user never
  has to install anything themselves. Trigger on requests to run metacellcnv,
  get CNV results from CellRanger output, or identify malignant cells from
  scRNA-seq via copy-number variation.
---

# metacellcnv bootstrap & run skill

You (the agent) are given nothing but a CellRanger output and have to make
`metacellcnv` work end to end yourself. Never tell the user "please install
X first" — whatever needs installing, you install it.

## 0. Parameters you take from the request

Read these from the user's request. Anything not given gets the stated
default; state the default you're using in one line and proceed — don't
stall on a clarifying question.

| Parameter | Required | Meaning | Default |
|---|---|---|---|
| `cellranger_dir` | yes | CellRanger output directory containing `filtered_feature_bc_matrix/` (one per sample, repeatable) | — |
| `sample_id` | — | Sample name for each `cellranger_dir` | inferred from the directory path |
| `genome` | yes, unless `gtf` given | Reference genome name, e.g. `GRCh38`, `GRCm39`, `CanFam3.1`, `CanFam6`, `Dog10K_Boxer_Tasha`, or any NCBI assembly name/accession | — |
| `gtf` | — | Path to an existing GTF. If given, skip auto-fetch entirely | none (auto-fetch) |
| `chromosome_map` | — | NCBI assembly report, or a 2-column chromosome-name table | auto-determined from where the GTF came from |
| `mito_chromosome` | — | Sequence ID of the mitochondrial genome | auto-detected from the GTF + known IDs |
| `markers` | — | Cell-type marker set: `dog` / `human` / `mouse`, or a CSV path | inferred from `genome` (ask the user if it can't be inferred) |
| `cnv_method` | — | `infercnv` (default, windowed moving average) or `bins` | `infercnv` |
| `cnv_reference` | — | `celltype` (tissue with normal cells mixed in) / `none` (pure tumor cell line) / `external` | inferred from sample composition (default to `celltype` if unclear) |
| `sample_kind` | — | `tissue` / `cell-line` | `auto` (pipeline's own default) |
| `n_cpus` | — | parallelism | number of logical cores on the machine |
| `out_dir` | — | where outputs go | `./metacellcnv_run/<timestamp>/` |
| `repo_url` | — | metacellcnv's GitHub URL | `https://github.com/takaho/metacellcnv.git` |
| `repo_ref` | — | branch/tag/commit | repo's default branch (usually `main`) |
| `python_version` | — | Python version to use | `3.11` (per the repo's `INSTALL.md`; `>=3.10` also works) |

## 1. Diagnose the environment (never touch the user's existing setup)

1. Check `python3 --version`, and whether `pyenv`, `uv`, `conda`/`mamba`,
   `asdf` are present via `command -v`.
2. Priority order for getting the required Python version (default 3.11,
   `>=3.10` acceptable): `uv` (`uv venv --python 3.11`) → `conda`/`mamba`
   (`conda create -n metacellcnv python=3.11`) → `pyenv` (install the tool
   itself with the user's OK if missing) → the system `python3` directly if
   it's already `>=3.10`, via `venv`.
3. Either way, create a **dedicated virtual environment** (e.g. named
   `metacellcnv`) and never modify the system Python or its packages. Reuse
   an existing environment of the same name rather than rebuilding it.
4. Report the chosen method and the environment's path in one line.

## 2. Fetch the repository

1. Set up a working directory for the source (e.g. `./_metacellcnv_src/`),
   kept separate from `cellranger_dir` and the output directory.
2. If already cloned there, `git fetch && git checkout <repo_ref>` to
   update; otherwise `git clone --depth 1 -b <repo_ref> <repo_url> .`
   (a normal shallow clone if no branch is given).
3. If `git` isn't usable or the network is closed to it, fall back to
   GitHub's tarball API
   (`https://github.com/<owner>/<repo>/archive/refs/heads/<ref>.tar.gz`)
   and extract it.
4. Once fetched, skim `README.md` and `INSTALL.md` — the CLI options may
   have moved on from what this skill describes (the upstream repo keeps
   changing), so sanity-check against them.

## 3. Install dependencies and self-verify

1. Inside the virtual environment, install per the repo's `INSTALL.md`
   (install the full production set, not the bare minimum, unless told
   otherwise):
   ```
   pip install numpy pandas scipy scikit-learn anndata scanpy \
               pydeseq2 scrublet plotly
   ```
   Add `plotly` unless the user explicitly wants to skip visualization.
   Only add `infercnvpy` / `matplotlib` if `--cnv-method bins` or the
   legacy matplotlib report is explicitly requested.
2. **Always self-verify**: run `python smoke_test.py` and confirm
   `ALL SMOKE TESTS PASSED` before touching real data. If it fails, stop —
   report the error as-is to the user rather than improvising a workaround.
3. Report a one-line summary once setup and self-verification are done
   (e.g. "environment ready (uv, python 3.11, metacellcnv HEAD, smoke test
   PASSED)").

## 4. Auto-resolve the genome annotation (GTF)

Skip this entirely and use the given `gtf` if one was provided.

1. **Check the cache first**: keep a dedicated cache directory (e.g.
   `reference/<genome>/`) in the working area, and reuse a previously
   fetched GTF (+ assembly report) for the same `genome` name rather than
   re-downloading (speeds up every run after the first).
2. **Human/mouse go through GENCODE first** (chromosome names already
   `chr1...`, so no `--chromosome-map` needed):
   - `genome` matching `GRCh38`/`hg38`/`human` →
     list `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/`
     and fetch the `gencode.v*.annotation.gtf.gz` found there.
   - `genome` matching `GRCm39`/`mm39`/`mouse` →
     same pattern against
     `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/latest_release/`.
   - For an explicitly named older assembly (`GRCh37`/`hg19`,
     `GRCm38`/`mm10`), look for the matching archived GENCODE release
     directory instead.
   - After fetching, `--gtf-gene-id` should be `gene_name`
     (`--gtf-gene-id auto` should detect this on its own — confirm with
     `--preflight-only`).
3. **Everything else (non-human/mouse organisms such as dog, or an
   explicit accession) goes through NCBI RefSeq/GenBank**:
   - If `genome` is already an NCBI accession (`GCF_*`/`GCA_*`), use it
     directly.
   - If it's an assembly name (e.g. `CanFam6`, `Dog10K_Boxer_Tasha`,
     `CanFam3.1`), resolve the current accession live via the NCBI
     Datasets API/CLI or `esearch -db assembly`. **Do not hardcode a
     name-to-accession table** — always resolve it at run time, since
     assemblies get superseded. This project has previously used
     `CanFam3.1` = `GCF_000002285.5` (useful as a reference point, but
     defer to whatever the live API returns if a newer assembly is
     wanted).
   - From the resolved accession, fetch the `*_genomic.gtf.gz` and
     `*_assembly_report.txt` under NCBI's genomes path convention
     (`.../genomes/all/GCF/NNN/NNN/NNN/GCF_XXXXXXXXX.X_<AssemblyName>/`),
     or use `datasets download genome accession <ACC> --include gtf,seq-report`
     if available to get both in one shot.
   - Pass the fetched `*_assembly_report.txt` directly as
     `--chromosome-map` (it maps RefSeq `NC_...` to `chr1...`). In this
     branch `--gtf-gene-id` is `gene` (RefSeq's attribute name, unlike
     GENCODE's `gene_name`).
4. For the mitochondrial sequence ID, use a known mapping if one applies
   (e.g. dog CanFam assemblies use `NC_002008.4`); otherwise rely on the
   pipeline's built-in auto-detection from the GTF.
5. Record the fetched files' paths and provenance (URL/accession) in the
   output directory's manifest, for reproducibility.

## 5. Validate the input (preflight)

Before the real run, always do one preflight pass per `cellranger_dir`:

```
python metacellcnv.py --cellranger-dir <dir> --gtf <gtf> --gtf-gene-id auto \
    [--chromosome-map <map>] --markers <markers> --preflight-only
```

Check the gene-name/GTF match rate, presence of mtDNA genes, presence of
marker genes, and chromosome naming. If the match rate is low or markers
barely hit, suspect a wrong `--gtf-gene-id` or a mismatched genome — find
and fix the cause before proceeding (never proceed past an unexplained
warning).

## 6. The real run (three stages per sample)

```
# 6a. Preprocessing
python metacellcnv.py --scanpy \
    --cellranger-dir <cellranger_dir> --sample-id <sample_id> \
    --gtf <gtf> --gtf-gene-id auto [--chromosome-map <map>] \
    --mito-chromosome <mito> \
    --out-dir <out_dir>/prep/<sample_id>

# 6b. metacells -> CNV -> malignancy calls -> DE
python metacellcnv.py \
    --cellranger-dir <cellranger_dir> --sample-id <sample_id> \
    --gtf <gtf> --gtf-gene-id auto [--chromosome-map <map>] \
    --mito-chromosome <mito> \
    --markers <markers> \
    --cnv-method <cnv_method> --cnv-reference <cnv_reference> \
    [--cnv-external-reference <path>] [--sample-kind <sample_kind>] \
    --n-cpus <n_cpus> \
    --out-dir <out_dir>/results/<sample_id>

# 6c. Interactive report
python metacellcnv.py --visualize \
    --results-dir <out_dir>/results/<sample_id> \
    --prep-dir <out_dir>/prep/<sample_id>
```

With multiple samples where a matched normal-tissue sample should serve as
an external CNV reference (`--cnv-reference external`), build that pooled
profile (the equivalent of `build_external_cnv_reference.py`) before this
step.

## 7. Progress reporting and ETA

At the start of each stage, print one line saying which stage is running
(e.g. `[2/3 samples] Case1: building metacells (seacells-fw-iters 50),
started 14:32`). For the ETA, don't fake confidence from the start:

1. Timestamp the start and end of every stage yourself.
2. Use the first light stage to finish (`--scanpy` preprocessing) as a
   throughput baseline: elapsed time per cell, on this machine, right now.
3. For the heavy stages (metacell construction, CNV inference),
   extrapolate linearly from that baseline for the same cell/gene counts,
   and label it explicitly as "a rough estimate from the previous stage's
   timing."
4. For stages whose stdout already carries progress (iteration counts,
   barcodes processed), read it and refresh the ETA on a cadence
   (every 30s-1min, or at each natural checkpoint). For a stage with no
   progress signal at all, say so plainly — "progress unknown, showing
   elapsed time only" — rather than inventing one.
5. With multiple samples, use the actual time taken by completed samples
   to estimate the next one (the first sample's ETA is the least certain;
   accuracy improves from the second sample on).
6. If a stage is running far longer than expected (roughly 2x or more),
   check whether it has actually stalled (CPU usage, whether the log is
   still advancing) and say so before continuing to report.

## 8. On failure

- `smoke_test.py` fails → do not proceed to real data; report the error
  as-is and stop.
- `--preflight-only` shows an extremely low match rate or almost no
  markers hitting → suspect `--gtf-gene-id` or a mismatched genome,
  find the actual cause, and don't proceed to the real run until it's
  resolved.
- The network is closed to repository/GTF fetching → say so explicitly,
  name the exact URL(s)/accession(s) needed, and ask the user to place
  the files manually (don't silently give up).

## 9. When it's done

- List, per sample, the paths to `prep/<sample_id>/`, `results/<sample_id>/`,
  and the self-contained HTML report.
- Write the actual parameters used (genome, accession, GTF provenance,
  every CLI option) to a manifest under `out_dir`, so a later run's diff
  is obvious.
- Leave the virtual environment, cloned repository, and downloaded genome
  files in place. On a later invocation that only changes
  `cellranger_dir` (a new sample), skip environment setup, repository
  fetch, and GTF fetch entirely and go straight to the real run — tell the
  user this makes repeat runs much faster.
