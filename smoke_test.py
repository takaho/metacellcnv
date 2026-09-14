"""Synthetic-data regression tests for the metacellcnv pipeline; run with `python smoke_test.py`."""
from __future__ import annotations

import gzip
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.io import mmwrite

import metacellcnv as P
import scrna_common as C

# Work dir is created next to this script (no hard-coded path per environment).
# Can also be overridden explicitly via SMOKE_TMP.
TMP = Path(os.environ.get("SMOKE_TMP") or (Path(__file__).resolve().parent / "_smoke"))
if TMP.exists():
    shutil.rmtree(TMP)
TMP.mkdir(parents=True)
# scanpy's read_10x_mtx(cache=True) keys its cache only on the path name, so
# rebuilding the synthetic data would otherwise load a stale matrix. Clear it every run.
if Path("cache").exists():
    shutil.rmtree("cache")

rng = np.random.default_rng(0)

# To actually exercise the v2.4 8-set panel, include markers for the 3 immune
# lineages plus endothelial and fibroblast. Fibroblast is there to test the
# "gets a type but is never used as a normal reference" path.
MK_TNK = ["CD3D", "CD3E", "CD3G", "CD2", "LCK", "NKG7", "GNLY", "KLRD1"]
MK_B = ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "BANK1"]
MK_MYE = ["LYZ", "CD14", "ITGAM", "FCGR3A", "CSF1R", "CD68", "C1QA", "AIF1"]
MK_ENDO = ["PECAM1", "CDH5", "VWF", "KDR", "CLDN5", "FLT1"]
MK_FIB = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "POSTN"]
MK_BLOCKS = [MK_TNK, MK_B, MK_MYE, MK_ENDO, MK_FIB]
GENES = (
    MK_TNK + MK_B + MK_MYE + MK_ENDO + MK_FIB
    + ["PTPRC", "ND1", "COX1", "CYTB"]
    + [f"LOC{100000 + i}" for i in range(262)]
)
assert len(GENES) == 300, len(GENES)
# Start position of each block within GENES (make_sample boosts each group accordingly)
_off, MK_SPANS = 0, []
for blk in MK_BLOCKS:
    MK_SPANS.append((_off, _off + len(blk)))
    _off += len(blk)


def make_sample(name: str, n_cells: int) -> Path:
    d = TMP / name / "outs" / "filtered_feature_bc_matrix"
    d.mkdir(parents=True)
    X = rng.negative_binomial(5, 0.3, size=(n_cells, len(GENES)))
    # Build a 5-group structure (T/NK, B, Myeloid, Endothelial, Fibroblast each boosted
    # in its own group). The relative-z-threshold defect only shows up when "all groups
    # are the same type", so here each group is given a distinct type instead, to see
    # whether the margin rule can assign a type to all 5 groups.
    k = n_cells // len(MK_SPANS)
    for gi, (lo, hi) in enumerate(MK_SPANS):
        beg = gi * k
        end = (gi + 1) * k if gi < len(MK_SPANS) - 1 else n_cells
        X[beg:end, lo:hi] += 60
    mmwrite(str(d / "matrix.mtx"), sparse.csr_matrix(X.T.astype(int)), field="integer")
    with open(d / "matrix.mtx", "rb") as fin, gzip.open(d / "matrix.mtx.gz", "wb") as fout:
        shutil.copyfileobj(fin, fout)
    (d / "matrix.mtx").unlink()
    with gzip.open(d / "features.tsv.gz", "wt") as f:
        for g in GENES:
            f.write(f"{g}\t{g}\tGene Expression\n")
    with gzip.open(d / "barcodes.tsv.gz", "wt") as f:
        for i in range(n_cells):
            f.write(f"{name}_BC{i:05d}-1\n")
    return d


dirs = [str(make_sample("CF6", 420)), str(make_sample("CF7", 390))]

# NCBI/RefSeq-style GTF (no gene_name, no chr prefix)
gtf = TMP / "genes.gtf"
with open(gtf, "wt") as f:
    for i, g in enumerate(GENES):
        chrom = f"NC_0065{83 + (i % 5)}.4"
        start = 1000 + i * 5000
        f.write(
            f'{chrom}\tGnomon\tgene\t{start}\t{start + 2000}\t.\t+\t.\t'
            f'gene_id "{g}"; gene "{g}"; gene_biotype "protein_coding";\n'
        )

print("\n########## preflight ##########")
info = P.preflight(dirs, str(gtf), None)
assert info["gene_id_attr"] == "gene", info["gene_id_attr"]
assert info["has_chr_prefix"] is False
assert info["exclude_chromosomes"] == ()

print("\n########## chromosome map from assembly report ##########")
report = TMP / "GCF_000002285.5_assembly_report.txt"
rows = ["# Assembly name: Dog10K_Boxer_Tasha", "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn"]
mols = ["1", "2", "3", "X", "MT"]
for i, mol in enumerate(mols):
    rows.append(
        f"{mol}\tassembled-molecule\t{mol}\tChromosome\tCM0000{i}.1\t=\tNC_0065{83 + i}.4"
    )
rows.append("UNPL1\tunplaced-scaffold\tna\tna\tAAEX0001.1\t=\tNW_00312345.1")
report.write_text("\n".join(rows) + "\n")
cmap = P.load_chromosome_map(str(report))
assert cmap["NC_006583.4"] == "chr1" and cmap["NC_006586.4"] == "chrX", cmap
assert cmap["NC_006587.4"] == "chrM"
assert "NW_00312345.1" not in cmap
norm = P.normalize_chromosome_names(pd.DataFrame({"chromosome": ["NC_006583.4", "NC_006586.4", "NW_9.1", np.nan]}), cmap)
assert list(norm["chromosome"][:3]) == ["chr1", "chrX", "chrNW_9.1"], list(norm["chromosome"])
assert P.check_exclude_chromosomes(("chrX", "chrY"), set(cmap.values())) == ("chrX",)

print("\n########## preflight (with assembly report) ##########")
info2 = P.preflight(dirs, str(gtf), None, str(report))
assert "chr1" in info2["normalized_chromosomes"], info2["normalized_chromosomes"]

print("\n########## GFF3 detection ##########")
gff3 = TMP / "bad.gff3"
gff3.write_text("NC_1\tX\tgene\t1\t100\t.\t+\t.\tID=gene-A;Name=A\n")
try:
    P.inspect_gtf(str(gff3), None)
    raise AssertionError("failed to detect GFF3")
except ValueError as e:
    print("OK:", str(e)[:60])

print("\n########## Step 1: load_and_preprocess ##########")
adata = P.load_and_preprocess(dirs, n_hvg=200, n_pcs=20)
P.ensure_sample_id(adata)
assert adata.n_obs == 810
assert list(adata.obs["sample_id"].cat.categories) == ["CF6", "CF7"]
assert "counts" in adata.layers
assert adata.obs_names.is_unique

print("\n########## Step 1a: leiden + annotate ##########")
sc.pp.neighbors(adata, use_rep="X_pca")
sc.tl.leiden(adata, resolution=0.5, key_added="coarse_cluster", flavor="igraph", n_iterations=2, directed=False)
adata.obs["cell_type"] = P.annotate_coarse_celltype(adata, cluster_key="coarse_cluster")
print(adata.obs["cell_type"].value_counts().to_dict())
assert not [c for c in adata.obs.columns if c.startswith("_score_")], "leftover temp columns"

print("\n########## v2.4: comparing decision rules (margin vs relative-z) ##########")
# leiden(resolution 0.5) collapses the 5 groups down to 3, so to actually exercise
# the decision rule itself we use the "planted" 5 groups as the cluster labels
# (recovered from the barcode index).
_sizes = {"CF6": 420, "CF7": 390}
def _planted(bc: str) -> str:
    _name, _rest = str(bc).split("_BC")
    _i = int(_rest.split("-")[0])
    _n = _sizes[_name]; _k = _n // len(MK_SPANS)
    _g = min(_i // _k, len(MK_SPANS) - 1)
    return f"g{_g}"
adata.obs["planted_group"] = pd.Categorical([_planted(b) for b in adata.obs_names])
print("  planted groups:", adata.obs["planted_group"].value_counts().sort_index().to_dict())
_res = {}
for _rule in ("margin", "relative-z"):
    _lab = P.annotate_coarse_celltype(adata.copy(), cluster_key="planted_group",
                                      rule=_rule, per_cluster_labels=False)
    _base = {v.split("_")[0] for v in set(_lab)} - {"Other"}
    _res[_rule] = {"types": _base, "n_ref": int(P.is_known_normal(_lab).sum()),
                   "labels": pd.Series(_lab).value_counts().to_dict()}
    print(f"  {_rule:<11}: types {sorted(_base)} / normal ref {_res[_rule]['n_ref']} cells"
          f" / {_res[_rule]['labels']}")
# When all 5 groups are distinct types (the ideal case) both rules should pass.
# Check that margin is not worse than relative-z.
assert len(_res["margin"]["types"]) == 5, _res["margin"]["types"]
assert _res["margin"]["n_ref"] >= _res["relative-z"]["n_ref"], _res
# Endothelial belongs in the normal reference; Fibroblast gets typed but stays out of it.
assert "Endothelial" in _res["margin"]["types"], _res["margin"]["types"]
assert "Fibroblast" in _res["margin"]["types"], "Fibroblast did not get a type"
assert "Fibroblast" not in P.KNOWN_NORMAL_CELLTYPES
assert not P.is_known_normal(pd.Series(["Fibroblast"]))[0], "Fibroblast is in the normal reference"
# The reference should be only the 4 groups other than Fibroblast
_n_per = 810 // 5
assert _res["margin"]["n_ref"] == 4 * _n_per, (_res["margin"]["n_ref"], 4 * _n_per)
print(f"  OK: all 5 types got typed, and the reference is the 4 groups excluding "
      f"Fibroblast ({_res['margin']['n_ref']} cells)")

print("\n########## v2.4: adding reference types via --normal-celltype ##########")
_before = list(P.KNOWN_NORMAL_CELLTYPES)
P.set_normal_celltypes(["Fibroblast"])
assert P.is_known_normal(pd.Series(["Fibroblast", "Fibroblast_3"])).all()
P.KNOWN_NORMAL_CELLTYPES[:] = _before          # restore for later tests
assert not P.is_known_normal(pd.Series(["Fibroblast"]))[0]
print("  OK: can add and restore")

print("\n########## v2.4: fail-fast on an invalid rule ##########")
try:
    P.annotate_coarse_celltype(adata.copy(), cluster_key="planted_group", rule="nonsense")
    raise AssertionError("no ValueError raised")
except ValueError as e:
    print("  OK:", str(e)[:60])
try:
    P.annotate_coarse_celltype(adata.copy(), cluster_key="planted_group", rule="absolute")
    raise AssertionError("no ValueError raised (absolute without min_score)")
except ValueError as e:
    print("  OK:", str(e)[:60])

print("\n########## v2.4: how many cells relative-z drops when all groups are one type ##########")
# Reproduces what actually happened on Case1: every cluster is Myeloid, and scores
# are uniformly high. The relative-z threshold only ever lets the top ~20% through,
# so most cells fall to Other.
_rng2 = np.random.default_rng(7)
_n2, _k2 = 700, 7
_M = _rng2.negative_binomial(5, 0.3, size=(_n2, len(GENES))).astype(float)
_grp = np.repeat(np.arange(_k2), _n2 // _k2)
_lo, _hi = MK_SPANS[2]                      # Myeloid block
for _g in range(_k2):
    # All groups are Myeloid; only the level varies slightly (matches Case1's 0.193-0.294).
    _M[_grp == _g, _lo:_hi] += 50 + 6 * _g
_uni = ad.AnnData(X=sparse.csr_matrix(_M))
_uni.obs_names = [f"U{i:05d}" for i in range(_n2)]
_uni.var_names = GENES
_uni.layers["counts"] = _uni.X.copy()
_uni.obs["planted_group"] = pd.Categorical([f"g{g}" for g in _grp])
sc.pp.normalize_total(_uni, target_sum=1e4); sc.pp.log1p(_uni)
_m = P.annotate_coarse_celltype(_uni.copy(), cluster_key="planted_group", rule="margin",
                                per_cluster_labels=False)
_z = P.annotate_coarse_celltype(_uni.copy(), cluster_key="planted_group", rule="relative-z",
                                per_cluster_labels=False)
_nm, _nz = int(P.is_known_normal(_m).sum()), int(P.is_known_normal(_z).sum())
_cm = pd.Series(_m).value_counts().to_dict()
_cz = pd.Series(_z).value_counts().to_dict()
print(f"  All {_k2} groups set to Myeloid (total {_n2} cells):")
print(f"    margin     : normal ref {_nm} cells {_cm}")
print(f"    relative-z : normal ref {_nz} cells {_cz}")
assert _nm > _nz, (_nm, _nz)
assert _nm == _n2, f"margin only picked up {_nm}/{_n2}"
assert _nz <= 0.35 * _n2, f"relative-z let through more than expected: {_nz}/{_n2}"
print("  OK: the margin rule fixes the relative-threshold defect (same-type groups mostly falling to Other)")

print("\n########## fallback markers ##########")
sub = adata[:, [g for g in adata.var_names if g.startswith("LOC")] + ["PTPRC"]].copy()
lab = P.annotate_coarse_celltype(sub, cluster_key="coarse_cluster")
assert "Leukocyte" in set(lab) or set(lab) == {"Other"}, set(lab)
print("fallback result:", pd.Series(lab).value_counts().to_dict())

print("\n########## fail-fast when all markers are missing ##########")
sub2 = adata[:, [g for g in adata.var_names if g.startswith("LOC")]].copy()
try:
    P.annotate_coarse_celltype(sub2, cluster_key="coarse_cluster")
    raise AssertionError("no ValueError raised")
except ValueError as e:
    print("OK:", str(e)[:50])

print("\n########## Step 1b: parametric_qc_filter ##########")
qc = P.parametric_qc_filter(adata, sample_key="sample_id", coarse_cluster_key="coarse_cluster")
assert not (qc["qc_fail"] & qc["outlier_high_count"] & ~qc["outlier_low_count"]
            & ~qc["outlier_low_feature"] & ~qc["outlier_high_pctmt"]).any(), \
    "high_count leaked into qc_fail (violates spec section 5-1)"
adata_qc = adata[~qc["qc_fail"].values].copy()

print("\n########## Step 1c: detect_doublets_cluster_aware (scrublet mocked) ##########")
def fake_scrublet(a, batch_key=None, expected_doublet_rate=0.05, random_state=0, **kw):
    r = np.random.default_rng(1)
    a.obs["doublet_score"] = r.random(a.n_obs)
    a.obs["predicted_doublet"] = a.obs["doublet_score"] > 0.9
    if batch_key is not None:
        assert batch_key in a.obs, "batch_key missing from work.obs"
        assert a.obs[batch_key].value_counts().min() >= 30, "small batches were not merged"
    return None

orig = sc.pp.scrublet
sc.pp.scrublet = fake_scrublet
try:
    calls = P.detect_doublets_cluster_aware(adata_qc, coarse_cluster_key="coarse_cluster")
finally:
    sc.pp.scrublet = orig
assert calls.index.equals(adata_qc.obs_names)
adata_qc.obs["doublet_score"] = calls["doublet_score"].values
adata_qc.obs["predicted_doublet"] = calls["predicted_doublet"].values

print("\n########## Step 1d: rescue_doublets_by_cnv_consistency ##########")
n = adata_qc.n_obs
clone = pd.Series(np.where(np.arange(n) % 2 == 0, "0", "1"), index=adata_qc.obs_names)
base = rng.normal(size=(2, 40))
cnv = np.vstack([base[int(c)] + rng.normal(scale=0.05, size=40) for c in clone])
rescued = P.rescue_doublets_by_cnv_consistency(
    adata_qc.obs[["doublet_score", "predicted_doublet"]], sparse.csr_matrix(cnv), clone
)
assert rescued.sum() <= adata_qc.obs["predicted_doublet"].sum()
print("doublets after rescue:", int(rescued.sum()), "/ original:", int(adata_qc.obs['predicted_doublet'].sum()))

# Row-count mismatch should raise
try:
    P.rescue_doublets_by_cnv_consistency(
        adata_qc.obs[["doublet_score", "predicted_doublet"]], cnv[:-1], clone
    )
    raise AssertionError("failed to detect row-count mismatch")
except ValueError as e:
    print("OK:", str(e)[:40])

# A single clone skips rescue
one = pd.Series(["0"] * n, index=adata_qc.obs_names)
r1 = P.rescue_doublets_by_cnv_consistency(
    adata_qc.obs[["doublet_score", "predicted_doublet"]], cnv, one
)
assert r1.equals(adata_qc.obs["predicted_doublet"].astype(bool))

adata_qc.obs["final_doublet"] = rescued.values
adata_f = adata_qc[~adata_qc.obs["final_doublet"].values].copy()

print("\n########## Step 3: aggregate_to_metacells ##########")
adata_f.obs["SEACell"] = [f"SEACell-{i % 12}" for i in range(adata_f.n_obs)]
mc = P.aggregate_to_metacells(adata_f, sample_key="sample_id", celltype_key="cell_type")
assert mc.n_obs == 12
assert P.looks_like_raw_counts(mc.layers["counts"])
assert mc.layers["counts"].sum() == adata_f.layers["counts"].sum()
print(mc.obs.head(3))

print("\n########## Step 4a: validating reference_key (spec section 4-#9) ##########")
mc_no_ct = mc.copy()
del mc_no_ct.obs["cell_type"]
try:
    P.run_cnv_branch(mc_no_ct, gtf_path=str(gtf), gtf_gene_id="gene", reference_key="cell_type")
    raise AssertionError("no KeyError raised")
except KeyError as e:
    print("OK:", str(e)[:70])

try:
    P.run_cnv_branch(mc, gtf_path=str(gtf), gtf_gene_id="gene",
                     reference_key="cell_type", reference_cat=["NoSuchType"])
    raise AssertionError("no ValueError raised")
except ValueError as e:
    print("OK:", str(e)[:70])

print("\n########## Step 4a: run_cnv_branch, actually run ##########")
cnv_mc = P.run_cnv_branch(
    mc, gtf_path=str(gtf), gtf_gene_id="gene",
    reference_key="cell_type", reference_cat=P.KNOWN_NORMAL_CELLTYPES,
    window_size=20, step=5, exclude_chromosomes=("chrX",), chromosome_map=cmap,
)
assert "chrX" not in set(cnv_mc.uns["cnv"]["chr_pos"]), cnv_mc.uns["cnv"]["chr_pos"]
print("CNV target chromosomes:", sorted(cnv_mc.uns["cnv"]["chr_pos"]))
assert "X_cnv" in cnv_mc.obsm
assert "cnv_score" in cnv_mc.obs and np.isfinite(cnv_mc.obs["cnv_score"]).all()
assert "cnv_leiden" in cnv_mc.obs
print("cnv_score:", cnv_mc.obs["cnv_score"].round(3).tolist()[:5])

print("\n########## single-cell-level CNV (equivalent to Step 1d) ##########")
prelim = P.run_cnv_branch(
    adata_f, gtf_path=str(gtf), gtf_gene_id="gene",
    reference_key=None, reference_cat=None, window_size=20, step=5,
    exclude_chromosomes=("chrX", "chrY"),   # doesn't match anything -> warns and continues
)
assert prelim.n_obs == adata_f.n_obs
r2 = P.rescue_doublets_by_cnv_consistency(
    adata_f.obs[["doublet_score", "predicted_doublet"]],
    prelim.obsm["X_cnv"], prelim.obs["cnv_leiden"],
)
print("doublets after rescue with real CNV:", int(r2.sum()))

print("\n########## GTF coordinate parser ##########")
pos = P.parse_gtf_gene_positions(str(gtf), gene_id_type="gene")
assert set(pos.columns) == {"chromosome", "start", "end"}
assert len(pos) == 300
print(pos.head(2))
P.check_exclude_chromosomes(("chrX", "chrY"), set(pos["chromosome"]))

print("\n########## Step 4b: run_de_branch ##########")
mc.obs["putative_malignant"] = np.where(np.arange(mc.n_obs) % 2 == 0, "malignant", "normal")
dds, res = P.run_de_branch(
    mc, "putative_malignant", "malignant", "normal",
    covariate_key="sample_id", min_total_counts=5, n_cpus=2,
)
print(res.head(3))
assert "padj" in res.columns

print("\n########## DE with a single sample (covariate has one level) ##########")
mc1 = mc.copy()
mc1.obs["sample_id"] = "CF6"
dds1, res1 = P.run_de_branch(
    mc1, "putative_malignant", "malignant", "normal",
    covariate_key="sample_id", min_total_counts=5, n_cpus=2,
)
assert "padj" in res1.columns

print("\n########## fail-fast on undersized groups ##########")
mc2 = mc.copy()
mc2.obs["putative_malignant"] = ["malignant"] + ["normal"] * (mc2.n_obs - 1)
try:
    P.run_de_branch(mc2, "putative_malignant", "malignant", "normal", covariate_key=None)
    raise AssertionError("no ValueError raised")
except ValueError as e:
    print("OK:", str(e)[:60])

print("\n########## environment.lock.txt ##########")
lock = P.write_environment_lock(TMP / "results")
assert lock and lock.exists() and "scanpy" in lock.read_text()

print("\n########## h5ad save ##########")
P._write_h5ad(mc, TMP / "results" / "metacells.h5ad")
assert (TMP / "results" / "metacells.h5ad").exists()
assert ad.read_h5ad(TMP / "results" / "metacells.h5ad").n_obs == 12

print("\n########## looks_like_raw_counts ##########")
assert P.looks_like_raw_counts(sparse.csr_matrix(np.array([[0, 3], [2, 0]])))
assert not P.looks_like_raw_counts(np.array([[0.0, 1.7], [2.3, 0.0]]))

print("\n########## v3.1: metacells_native (no SEACells dependency) ##########")
import sys as _sys, types as _types, time as _time
import metacells_native as MC
import scipy.sparse as _sp

_rng = np.random.default_rng(0)
_X = np.vstack([_rng.normal(c, 1.0, (300, 20)) for c in (0, 4, 8, 12)])
_ad4 = ad.AnnData(np.zeros((1200, 5), dtype=np.float32))
_ad4.obsm["X_pca"] = _X
_ad4.obs_names = [f"c{i}" for i in range(1200)]
_ad4.obs["truth"] = np.repeat(["a", "b", "c", "d"], 300)
_M = MC.adaptive_rbf_kernel(_ad4, k=15)   # leaves obsp['distances'] on _ad4

# --- (1) the kernel is symmetric with a diagonal of 1 (self as a neighbor) ---
assert abs(_M - _M.T).max() == 0.0, "kernel is not symmetric"
assert np.allclose(_M.diagonal(), 1.0), _M.diagonal()[:5]
print(f"  OK: kernel is symmetric with diagonal 1 (nnz {_M.nnz:,})")

# --- (2) the bandwidth matches SEACells' kth_neighbor_distance ---
# The original implementation picks "the element at descending rank num_nonzero -
# kth", which is equivalent to the kth-smallest nonzero in ascending order. Confirm
# the equivalence against real data shapes.
_dist = _ad4.obsp["distances"].tocsr()
def _seacells_kth(distances, kth, i):
    row = distances[i, :].toarray().ravel()
    nz = int(np.sum(row > 0))
    mask = np.argsort(np.argsort(-row)) == nz - kth
    return float(np.linalg.norm(row[mask]))
_own = MC._kth_smallest_nonzero(_dist, 15 // 2)
_ref = np.array([_seacells_kth(_dist, 15 // 2, i) for i in range(0, 1200, 37)])
assert np.allclose(_own[::37], _ref), np.abs(_own[::37] - _ref).max()
print(f"  OK: adaptive bandwidth matches SEACells' definition (max diff {np.abs(_own[::37]-_ref).max():.2e})")

# --- (3) the factorized K@B matches the explicit K@B ---
_K = (_M @ _M.T).toarray()
_B = _rng.random((1200, 30)); _B /= _B.sum(0)
_d = np.abs(MC._KX(_M, _B) - _K @ _B).max()
assert _d < 1e-10, _d
print(f"  OK: K@B = M@(M@B) matches (max diff {_d:.2e}, never materializes n×n)")

# --- (4) RSS matches the naive dense computation (never materializes n×n) ---
_A = _rng.random((30, 1200)); _A /= _A.sum(0)
_naive = float(np.linalg.norm(_M.toarray() - (_M.toarray() @ _B) @ _A))
_fast = MC._rss(_M, _A, _B, float(_M.multiply(_M).sum()))
assert abs(_naive - _fast) / _naive < 1e-10, (_naive, _fast)
print(f"  OK: RSS is exact (naive {_naive:.6f} / fast {_fast:.6f})")

# --- (5) different seeds converge to the same solution (SEACells has no seed) ---
_labs = []
for _sd in (0, 1):
    _a = _ad4.copy()
    _m = MC.MetacellModel(n_metacells=30, seed=_sd, verbose=False).fit(
        _a, max_iter=20, min_iter=5)
    _labs.append(_a.obs["SEACell"].astype(str).values)
from sklearn.metrics import adjusted_rand_score as _ARI
_ari = _ARI(_labs[0], _labs[1])
assert _ari > 0.8, f"solution changes a lot with seed: ARI {_ari}"
print(f"  OK: ARI between seed 0 and 1 is {_ari:.3f} (deterministic init)")

# --- (6) degenerate size distributions can be detected ---
_even = np.repeat([f"m{i}" for i in range(30)], 40)
_degen = np.array([f"m{i}" for i in range(29)] + ["m29"] * (1200 - 29))
_be, _bd = MC.size_balance(_even), MC.size_balance(_degen)
assert _be["max_share"] < 0.05 and _be["gini"] < 0.05, _be
assert _bd["max_share"] > 0.9 and _bd["n_singleton"] == 29, _bd
print(f"  OK: degeneracy detection (even: max share {_be['max_share']:.1%} /"
      f" degenerate: {_bd['max_share']:.1%}, {_bd['n_singleton']} singleton metacells)")

# --- (7) line search works but warns (disabled by default) ---
assert MC.MetacellModel(n_metacells=5).line_search is False, "default has line search on"
print("  OK: default is fixed step size (line search lowers RSS but wrecks size balance)")

# --- (8) diffusion components drop the λ=1 eigenvalue (diverges on a disconnected graph) ---
_dc = MC.diffusion_components(_X)
assert np.isfinite(_dc).all() and _dc.var(0).max() < 1e3, _dc.var(0)
print(f"  OK: diffusion components are finite (max variance {_dc.var(0).max():.4g}, λ=1 removed)")

# --- (9) evaluation metrics work without palantir ---
_a = _ad4.copy()
MC.MetacellModel(n_metacells=30, seed=0, verbose=False).fit(_a, max_iter=20, min_iter=5)
_c = MC.compactness(_a); _s = MC.separation(_a); _pu = MC.celltype_purity(_a, "truth")
assert len(_c) == len(_s) == 30 and _pu["truth_purity"].median() > 0.9
print(f"  OK: compactness/separation/purity (median purity {_pu['truth_purity'].median():.3f})")

# --- (10) the pipeline does not import SEACells ---
# Parse only the import statements via AST so we don't pick up explanatory text in
# strings or comments.
import ast as _ast
_tree = _ast.parse(Path("metacellcnv.py").read_text(encoding="utf-8"))
_imported = set()
for _n in _ast.walk(_tree):
    if isinstance(_n, _ast.Import):
        _imported |= {a.name.split(".")[0] for a in _n.names}
        _imported |= {(a.asname or a.name).split(".")[0] for a in _n.names}
    elif isinstance(_n, _ast.ImportFrom) and _n.module:
        _imported.add(_n.module.split(".")[0])
for _bad in ("SEACells", "infercnvpy", "palantir"):
    assert _bad not in _imported, f"pipeline imports {_bad}"
for _mod in ("metacells_native", "cnv_native"):
    assert _mod in _imported, f"pipeline does not import {_mod}"
print("  OK: the pipeline does not import SEACells / infercnvpy")


# --- (11) the pipeline's wrapper functions work (names and obs keys stay compatible) ---
_a = _ad4.copy()
_mdl = P.build_seacells(_a, cells_per_metacell=40, min_iter=3, max_iter=10,
                        verbose_iterations=False)
assert _a.obs["SEACell"].nunique() == 30, _a.obs["SEACell"].nunique()
_met = P.evaluate_metacells(_a, _mdl, celltype_key="truth")
assert {"SEACell", "compactness", "separation"} <= set(_met.columns), list(_met.columns)
assert len(_met) == 30
# Passing line_search should still get disabled with a warning
_a2 = _ad4.copy()
_m2 = P.build_seacells(_a2, cells_per_metacell=40, min_iter=3, max_iter=6,
                       line_search=True, verbose_iterations=False)
assert _m2.line_search is False, "line search was not disabled"
print(f"  OK: build_seacells / evaluate_metacells work, and line_search=True is "
      f"disabled ({_a.obs['SEACell'].nunique()} metacells)")


print("\n########## v3.1: cnv_native (no infercnvpy dependency) ##########")
import cnv_native as CV
from scipy.stats import norm as _norm
def stats_sf(x):
    return _norm.sf(x)

# --- (1) the pyramid-window moving average matches infercnvpy's definition ---
_x = _rng.normal(size=(6, 500))
def _ref_running_mean(x, n=100, step=10):
    r = np.arange(1, n + 1)
    pyr = np.minimum(r, r[::-1])
    sm = np.apply_along_axis(lambda row: np.convolve(row, pyr, mode="valid"), 1, x) / pyr.sum()
    return sm[:, np.arange(0, sm.shape[1], step)]
_d = np.abs(CV.running_mean_pyramid(_x, 100, 10) - _ref_running_mean(_x)).max()
assert _d < 1e-12, _d
print(f"  OK: pyramid-window moving average matches (max diff {_d:.2e})")

# --- (2) bounded mode: values inside the reference range become 0 ---
_v = np.array([[0.0, 1.0, 2.0, 3.0]])
_r2 = np.array([[1.0, 1.0, 1.0, 1.0], [2.0, 2.0, 2.0, 2.0]])
_bc = CV._bounded_center(_v, _r2)
assert np.allclose(_bc, [[-1.0, 0.0, 0.0, 1.0]]), _bc
assert np.allclose(CV._bounded_center(_v, _r2[:1]), _v - 1.0)
print("  OK: bounded mode (only values outside the reference min/max become logFC)")

# --- (3) synthetic AnnData with coordinates ---
_ng, _nc = 900, 200
_gp = pd.DataFrame({
    "chromosome": np.repeat(["chr1", "chr2", "chr3"], _ng // 3),
    "start": np.tile(np.arange(_ng // 3) * 1000, 3).astype(float),
})
_gp["end"] = _gp["start"] + 500
_gp.index = [f"g{i}" for i in range(_ng)]
_base = _rng.lognormal(0.0, 1.0, _ng)
_depth = _rng.integers(3000, 12000, _nc)
_gain = np.zeros(_nc, dtype=bool); _gain[:60] = True     # double all of chr2
_prob = np.tile(_base, (_nc, 1))
_prob[np.ix_(_gain, np.where(_gp["chromosome"].values == "chr2")[0])] *= 2.0
_prob /= _prob.sum(1, keepdims=True)
_cnt = np.vstack([_rng.multinomial(int(_depth[i]), _prob[i]) for i in range(_nc)])
_cad = ad.AnnData(sparse.csr_matrix(_cnt.astype(np.float32)))
_cad.var_names = _gp.index; _cad.obs_names = [f"c{i}" for i in range(_nc)]
_cad.var = _cad.var.join(_gp)
_cad.layers["counts"] = _cad.X.copy()
_cad.obs["grp"] = np.where(_gain, "malignant", "normal")

# --- (4) infercnv-equivalent output is the same regardless of block size
#     (infercnvpy's result changes with chunksize) ---
_out = []
for _blk in (_nc, 50, 17):
    _c2 = _cad.copy()
    _c2.X = _c2.X.copy()
    import scanpy as _sc
    _sc.pp.normalize_total(_c2, target_sum=1e4); _sc.pp.log1p(_c2)
    CV.infercnv_scores(_c2, reference_key="grp", reference_cat=["normal"],
                       window_size=50, step=5, block=_blk,
                       exclude_chromosomes=())
    _out.append(np.asarray(_c2.obsm["X_cnv"].todense()))
_dd = max(float(np.abs(_out[i] - _out[0]).max()) for i in range(len(_out)))
assert _dd == 0.0, f"result changes with block: {_dd}"
print(f"  OK: result is identical across block sizes {_nc}/50/17 (diff {_dd})")

# --- (5) count-balanced bins have less variance in weight ---
_w = np.asarray(_cad.layers["counts"].sum(0)).ravel().astype(float)
_gb_g, _if_g = CV.genomic_bins(_cad.var, 100, gene_weight=_w, bin_by="genes",
                               exclude=())
_gb_c, _if_c = CV.genomic_bins(_cad.var, 100, gene_weight=_w, bin_by="counts",
                               exclude=())
_cv_g = _if_g["weight"].std() / _if_g["weight"].mean()
_cv_c = _if_c["weight"].std() / _if_c["weight"].mean()
assert _cv_c < _cv_g, (_cv_g, _cv_c)
print(f"  OK: bin-weight CV goes from {_cv_g:.3f} (gene-balanced) to {_cv_c:.3f} (count-balanced)")

# --- (6) distribution model: binomial overstates z on overdispersed data ---
_n = np.full(400, 5000.0)
_p = 0.02
_k_bb = _rng.binomial(5000, _rng.beta(2, 98, 400))       # overdispersed
_W = _k_bb.reshape(-1, 1).astype(float)
_exp = np.outer(_n, np.array([_p]))
_z_bi, _ = CV._dispersion_z(_W, _exp, _n, np.array([_p]), "binom")
_z_bb, _ = CV._dispersion_z(_W, _exp, _n, np.array([_p]), "betabinom")
_z_nb, _ = CV._dispersion_z(_W, _exp, _n, np.array([_p]), "nb")
assert _z_bi.std() > 3 * _z_bb.std(), (_z_bi.std(), _z_bb.std())
print(f"  OK: SD(z) on overdispersed data — binomial {_z_bi.std():.2f} /"
      f" beta-binomial {_z_bb.std():.2f} / negative binomial {_z_nb.std():.2f}")

# --- (7) negative-binomial alpha can be recovered by the method of moments ---
_mu = np.full(2000, 50.0)
_alpha_true = 0.5
_kk = _rng.negative_binomial(1.0 / _alpha_true, 1.0 / (1.0 + _alpha_true * 50.0), 2000).astype(float)
_ah = CV._nb_alpha(_kk, _mu, trim=1.0)
assert 0.2 < _ah < 1.2, _ah
print(f"  OK: recovered negative-binomial alpha (true {_alpha_true} -> estimate {_ah:.3f})")

# --- (8) supplying a reference raises detection power (no longer estimating
#     overdispersion together with the signal) ---
_pw = {}
for _tag, _rm in (("no reference", None), ("normal reference", ~_gain)):
    _r = CV.bin_composition_cnv(_cad.copy(), whole_chromosome=True, bin_by="genes",
                                dist="betabinom", layer="counts",
                                reference_mask=_rm, n_control=0, n_null_shuffle=3,
                                exclude_chromosomes=(), seed=0, inplace=False)
    _info = _r["bins"]
    _j = _info.index.get_loc(_info.index[_info["chromosome"] == "chr2"][0])
    _pw[_tag] = float((_r["q"][_gain, _j] < 0.05).mean())
assert _pw["normal reference"] >= _pw["no reference"], _pw
print(f"  OK: power to detect the planted chr2 event — no reference {_pw['no reference']:.3f}"
      f" -> normal reference {_pw['normal reference']:.3f}")

# --- (9) an expression-program control is returned ---
_r = CV.bin_composition_cnv(_cad.copy(), bin_genes=100, reference_mask=~_gain,
                            layer="counts", n_control=5, n_null_shuffle=3,
                            exclude_chromosomes=(), seed=0, inplace=False)
assert _r["control_amplitude"] > 0 and _r["amplitude_ratio"] > 1.0, _r["amplitude_ratio"]
print(f"  OK: amplitude {_r['amplitude']:.4f} / control {_r['control_amplitude']:.4f}"
      f" -> ratio {_r['amplitude_ratio']:.2f}x")

# --- (10) empirical p-values are the default and are better calibrated in the tail
#     than the normal approximation ---
# Build overdispersed (null) counts and compare the actual false-positive rate at a
# nominal alpha.
_nc2, _nb2 = 300, 40
_dep = _rng.integers(20000, 60000, _nc2)
_pb2 = _rng.dirichlet(np.full(_nb2, 3.0))
# Add "program differences" per metacell = perturb the bin proportions themselves
_W2 = np.vstack([_rng.multinomial(int(_dep[i]),
                                  _rng.dirichlet(_pb2 * 60.0)) for i in range(_nc2)])
_n2 = _W2.sum(1).astype(float)
_p2 = _W2.sum(0) / _W2.sum()
_z2, _ = CV._dispersion_z(_W2.astype(float), np.outer(_n2, _p2), _n2, _p2, "betabinom")
_z2 = _z2 / max(_z2.std(), 1e-9)
_fpr_norm = float(np.mean(2 * stats_sf(np.abs(_z2)) < 0.001))
assert _fpr_norm > 0.003, f"normal approximation does not break down on this synthetic data: {_fpr_norm}"
print(f"  OK: on overdispersed null data, the normal approximation's alpha=0.001 false-positive"
      f" rate is {_fpr_norm:.4f} ({_fpr_norm/0.001:.0f}x) — the reason empirical p is the default")
assert CV.DEFAULT_P_METHOD == "empirical", CV.DEFAULT_P_METHOD
_r2 = CV.bin_composition_cnv(_cad.copy(), bin_genes=100, reference_mask=~_gain,
                             layer="counts", n_control=0, n_null_shuffle=30,
                             exclude_chromosomes=(), seed=0, inplace=False)
assert _r2["p_method"] == "empirical" and _r2["resolution"] is not None
assert _r2["p_value"].min() >= _r2["resolution"] - 1e-12, "p goes below the resolution floor"
print(f"  OK: default is empirical p (resolution floor {_r2['resolution']:.2e},"
      f" min p {_r2['p_value'].min():.2e})")

# --- (11) contiguity test: whole-chromosome-only changes should not support
#     subchromosomal structure ---
_ct = CV.contiguity_test(_cad.copy(), bin_genes=100, n_perm=8, layer="counts",
                         exclude_chromosomes=(), seed=0)
assert not _ct["supports_segmental"], _ct
print(f"  OK: whole-chromosome-only data does not support subchromosomal structure"
      f" (ratio {_ct['ratio']:.2f}, p={_ct['p_value']:.3f})")

# --- (12) CNV embedding works with scanpy alone ---
_c3 = _cad.copy()
CV.bin_composition_cnv(_c3, bin_genes=100, reference_mask=~_gain, layer="counts",
                       n_control=0, n_null_shuffle=0, exclude_chromosomes=())
CV.cnv_embedding(_c3, key="cnv", n_comps=10, cluster_key="cnv_leiden")
assert "X_cnv_pca" in _c3.obsm and _c3.obs["cnv_leiden"].nunique() >= 2
print(f"  OK: CNV embedding (Leiden {_c3.obs['cnv_leiden'].nunique()} clusters)")


print("\n########## metacellcnv_scanpy v1.0: unified Seurat-replacement preprocessing ##########")
import json as _json
import gzip as _gz
import scipy.io as _sio
import metacellcnv_scanpy as PS

# --- (1) mtx.gz round-trips (both orientations, integer field) ---
_rgp = np.random.default_rng(7)
_Xp = sparse.csr_matrix(_rgp.poisson(0.4, size=(60, 25)).astype(np.int64))
for _orient, _shape in [("cells-genes", (60, 25)), ("genes-cells", (25, 60))]:
    _pth = TMP / f"m_{_orient}.mtx.gz"
    PS.write_mtx_gz(_Xp, _pth, orient=_orient, integer=True)
    with _gz.open(_pth, "rb") as _fh:
        _Y = _sio.mmread(_fh)
    assert _Y.shape == _shape, (_orient, _Y.shape)
    _back = sparse.csr_matrix(_Y).T if _orient == "genes-cells" else sparse.csr_matrix(_Y)
    assert abs(_back - _Xp).max() == 0, _orient
    assert np.issubdtype(_Y.data.dtype, np.integer), _Y.data.dtype
print("  OK: mtx.gz round-trips in both orientations, and raw counts are saved as integers")

# integer=True on non-integer data should stop (not silently round)
try:
    PS.write_mtx_gz(_Xp.astype(float) * 1.5, TMP / "bad.mtx.gz", integer=True)
    raise AssertionError("was able to write a non-integer array with integer=True")
except ValueError as _e:
    print("  OK: writing non-integer data with integer=True stops:", str(_e)[:40])

# --- (2) the containment metric separates "metacells built within a cluster" from
#     "metacells that ignore cluster boundaries" ---
_bc = [f"C{i:04d}" for i in range(600)]
_clu = pd.Series([str(i % 5) for i in range(600)], index=_bc)
_within, _nxt = np.empty(600, dtype=object), 0
for _c in sorted(_clu.unique()):
    _idx = np.where(_clu.values == _c)[0]
    for _chunk in np.array_split(_idx, 4):
        _within[_chunk] = f"MC{_nxt}"; _nxt += 1
_res_in = PS.metacell_cluster_containment(pd.Series(_within, index=_bc), _clu)
_rand = pd.Series([f"MC{i}" for i in _rgp.integers(0, _nxt, 600)], index=_bc)
_res_rd = PS.metacell_cluster_containment(_rand, _clu)
_a = _res_in["summary"]["containment_weighted_mean"]
_b = _res_rd["summary"]["containment_weighted_mean"]
print(f"  metacells built within a cluster: {_a:.3f} / metacells ignoring clusters: {_b:.3f}")
assert _a > 0.99 and _b < 0.45, (_a, _b)
assert _res_in["summary"]["frac_fully_contained"] == 1.0
assert set(_res_in["per_metacell"].columns) == {
    "n_cells", "containment", "dominant_cluster", "n_clusters_spanned"}
PS.report_containment(_res_rd)
print("  OK: containment clearly separates full containment (1.00) from random (0.4)")

# --- (3) auto-detection of the metacell column, and error handling ---
for _col in ("SEACell", "metacell", "metacell_id"):
    _df = pd.DataFrame({_col: _within}, index=_bc)
    assert PS.metacell_cluster_containment(_df, _clu)["summary"]["n_metacells"] == _nxt, _col
try:
    PS.metacell_cluster_containment(pd.DataFrame({"foo": _within}, index=_bc), _clu)
    raise AssertionError("passed despite a wrong column name")
except ValueError as _e:
    print("  OK: stops and names the expected column names when none match:", str(_e)[:40])
# Should stop if not a single barcode matches (never silently return an empty result)
try:
    PS.metacell_cluster_containment(pd.Series(_within, index=[f"X{i}" for i in range(600)]), _clu)
    raise AssertionError("passed despite mismatched barcodes")
except ValueError as _e:
    print("  OK: stops when barcodes don't match at all:", str(_e)[:40])
# A partial match should warn and compute over the intersection
_part = PS.metacell_cluster_containment(
    pd.Series(_within[:400], index=_bc[:400]), _clu)
assert _part["summary"]["n_cells"] == 400, _part["summary"]
print("  OK: a partial match warns and computes over the intersection")

# --- (4) --seurat-compat argument mapping ---
_a_def = PS.parse_args(["--cellranger-dir", "x"])
assert _a_def.norm_target_sum is None and _a_def.hvg_flavor == "seurat"
assert _a_def.scale_before_pca is False and _a_def.umap_components == 3
assert _a_def.orient == "cells-genes" and _a_def.no_refit is False
_a_sc = PS.parse_args(["--cellranger-dir", "x", "--seurat-compat"])
assert _a_sc.norm_target_sum == 1e4, _a_sc.norm_target_sum
assert _a_sc.hvg_flavor == "seurat_v3" and _a_sc.scale_before_pca is True
assert _a_sc.umap_neighbors == 30 and _a_sc.umap_min_dist == 0.3
assert _a_sc.min_cells_per_gene == 3
assert _a_sc.cluster_algo == "leiden", "no SLM available, so stays Leiden"
print("  OK: defaults match the current pipeline; --seurat-compat switches to Seurat defaults")

# --- (5) end to end: pick up dog-style mtDNA (no MT- prefix) via the GTF sequence ID ---
# GENES already contains ND1/COX1/CYTB as "genes on nuclear chromosomes". This is
# exactly the situation where name-only mtDNA detection (MT- prefix or symbol match)
# gets it wrong — a real problem for dogs. Here only the remaining 10 genes are placed
# on the mtDNA sequence NC_002008.4, to confirm that "sequence-ID-based detection wins
# over name-based detection".
_MITO_ALL = ["ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
             "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB"]
_NUCLEAR_LOOKALIKE = [_g for _g in _MITO_ALL if _g in GENES]
_MITO = [_g for _g in _MITO_ALL if _g not in GENES]
assert len(_NUCLEAR_LOOKALIKE) == 3 and len(_MITO) == 10, (_NUCLEAR_LOOKALIKE, _MITO)
_pd_dir = TMP / "prep_src" / "filtered_feature_bc_matrix"
_pd_dir.mkdir(parents=True)
_ng, _nc = len(GENES) + len(_MITO), 700
_Xd = rng.negative_binomial(5, 0.3, size=(_nc, _ng))
_kk = _nc // len(MK_SPANS)
for _gi, (_lo, _hi) in enumerate(MK_SPANS):
    _Xd[_gi * _kk:(_gi + 1) * _kk, _lo:_hi] += 60
_Xd[:, len(GENES):] += 25                      # make sure mtDNA is expressed
mmwrite(str(_pd_dir / "matrix.mtx"), sparse.csr_matrix(_Xd.T.astype(int)), field="integer")
with open(_pd_dir / "matrix.mtx", "rb") as _fi, _gz.open(_pd_dir / "matrix.mtx.gz", "wb") as _fo:
    shutil.copyfileobj(_fi, _fo)
(_pd_dir / "matrix.mtx").unlink()
with _gz.open(_pd_dir / "features.tsv.gz", "wt") as _f:
    for _g in list(GENES) + _MITO:
        _f.write(f"{_g}\t{_g}\tGene Expression\n")
with _gz.open(_pd_dir / "barcodes.tsv.gz", "wt") as _f:
    for _i in range(_nc):
        _f.write(f"PREP_BC{_i:05d}-1\n")
_gtf_d = TMP / "genes_dog.gtf"
with open(_gtf_d, "wt") as _f:
    for _i, _g in enumerate(GENES):
        _c = f"NC_0065{83 + (_i % 5)}.4"; _s = 1000 + _i * 5000
        _f.write(f'{_c}\tGnomon\tgene\t{_s}\t{_s + 2000}\t.\t+\t.\tgene_id "{_g}"; gene "{_g}";\n')
    for _j, _g in enumerate(_MITO):     # dog mtDNA — names don't carry an MT- prefix
        _s = 100 + _j * 1200
        _f.write(f'NC_002008.4\tGnomon\tgene\t{_s}\t{_s + 900}\t.\t+\t.\tgene_id "{_g}"; gene "{_g}";\n')

if Path("cache").exists():
    shutil.rmtree("cache")
_out_prep = TMP / "prep_out"
_args = PS.parse_args([
    "--cellranger-dir", str(_pd_dir), "--sample-id", "PREP",
    "--gtf", str(_gtf_d), "--gtf-gene-id", "gene",
    "--mito-chromosome", "NC_002008.4",
    "--n-hvg", "120", "--n-pcs", "15", "--umap-components", "3",
    "--out-dir", str(_out_prep), "--no-doublet", "--no-pctmt-filter",
])
_ad_prep, _man = PS.run_prep(_args)

assert _man["matrix_orientation"] == "cells-genes"
assert _man["n_cells"] == _ad_prep.n_obs and _man["n_genes"] == _ad_prep.n_vars
assert _man["versions"]["scrna_common"] == C.__version__, _man["versions"]["scrna_common"]
for _key in ("barcodes", "features", "matrix_raw", "matrix_lognorm",
             "pca", "umap3d", "clusters", "qc_metrics"):
    assert _key in _man["files"], _key
    _fp = _out_prep / _man["files"][_key]["path"]
    assert _fp.exists(), _fp
    assert PS.sha256_file(_fp) == _man["files"][_key]["sha256"], _key
print(f"  OK: all 8 artifacts + manifest are present, and every sha256 matches ({_man['n_cells']} cells)")

_ft = pd.read_csv(_out_prep / "features.tsv.gz", sep="\t", index_col=0)
assert int(_ft["mt"].sum()) == len(_MITO), int(_ft["mt"].sum())
assert not any(str(_g).upper().startswith("MT-") for _g in _ft.index[_ft["mt"]]), \
    "the no-MT--prefix assumption has broken"
assert set(_ft.index[_ft["mt"]]) == set(_MITO)
# The 3 genes that merely look mtDNA-ish by name, but sit on a nuclear chromosome,
# must not be picked up as mt
for _g in _NUCLEAR_LOOKALIKE:
    assert not bool(_ft.loc[_g, "mt"]), f"{_g} was misclassified as mtDNA by name alone"
    assert _ft.loc[_g, "chromosome"] != "NC_002008.4", _g
assert int((_ft["mt"] & _ft["hvg_for_pca"]).sum()) == 0, "mtDNA still present in HVGs"
assert _ft.loc[_MITO[0], "chromosome"] == "NC_002008.4"
print(f"  OK: identified {len(_MITO)} mtDNA genes with no MT- prefix via the GTF sequence ID and"
      f" excluded them from HVGs. The name-only lookalikes {_NUCLEAR_LOOKALIKE} are not misclassified")

_bcp = [_l.strip() for _l in _gz.open(_out_prep / "barcodes.tsv.gz", "rt")]
for _nm, _fn in [("pca", "pca.tsv.gz"), ("umap", "umap3d.tsv.gz"),
                 ("clusters", "clusters.tsv.gz"), ("qc", "qc_metrics.tsv.gz")]:
    _df = pd.read_csv(_out_prep / _fn, sep="\t", index_col=0)
    assert list(_df.index) == _bcp, f"{_nm} barcode order does not match barcodes.tsv.gz"
assert len(_bcp) == _man["n_cells"]
print("  OK: every artifact shares the same barcode set in the same order")

with _gz.open(_out_prep / "matrix_raw.mtx.gz", "rb") as _fh:
    _R = _sio.mmread(_fh)
with _gz.open(_out_prep / "matrix_lognorm.mtx.gz", "rb") as _fh:
    _L = sparse.csr_matrix(_sio.mmread(_fh))
assert _R.shape == (_man["n_cells"], _man["n_genes"]) == _L.shape
assert np.allclose(_R.data, np.rint(_R.data))
_ex = _L.copy(); _ex.data = np.expm1(_ex.data)
_tot = np.asarray(_ex.sum(axis=1)).ravel()
assert _tot.std() / _tot.mean() < 1e-6, f"post-normalization row totals are not constant, CV={_tot.std()/_tot.mean():.2e}"
print(f"  OK: matrix shapes match the manifest, and post-normalization row totals are constant ({np.median(_tot):.0f})")

_umap = pd.read_csv(_out_prep / "umap3d.tsv.gz", sep="\t", index_col=0)
assert list(_umap.columns) == ["UMAP1", "UMAP2", "UMAP3"], list(_umap.columns)
_clup = pd.read_csv(_out_prep / "clusters.tsv.gz", sep="\t", index_col=0)
assert "cluster" in _clup.columns and "cluster_provisional" in _clup.columns
assert _clup["cluster"].nunique() >= 2, _clup["cluster"].nunique()
print(f"  OK: UMAP is 3D / {_clup['cluster'].nunique()} clusters"
      f" ({_clup['cluster_provisional'].nunique()} provisional clusters also recorded)")

# --- (6) containment between an actual metacell assignment and the output clusters
#     can be measured ---
_mc_path = TMP / "mc_for_prep.csv"
_lab = _clup["cluster"].astype(str)
_mcv, _mn = np.empty(len(_lab), dtype=object), 0
for _c in sorted(_lab.unique()):
    _idx = np.where(_lab.values == _c)[0]
    for _chunk in np.array_split(_idx, max(1, len(_idx) // 40)):
        _mcv[_chunk] = f"MC{_mn}"; _mn += 1
pd.DataFrame({"SEACell": _mcv}, index=_lab.index).rename_axis("barcode").to_csv(_mc_path)
_cargs = PS.parse_args(["--containment", str(_mc_path),
                        "--clusters", str(_out_prep / "clusters.tsv.gz")])
_cres = PS.run_containment(_cargs)
assert _cres["summary"]["containment_weighted_mean"] > 0.99
assert not Path("prep").exists(), "wrote output despite --out-dir not being given"
print("  OK: containment can be measured against the saved clusters.tsv.gz, and nothing is written without --out-dir")

# --- (7) --no-refit reuses the pre-QC PCA as is (compatible with the current pipeline) ---
if Path("cache").exists():
    shutil.rmtree("cache")
_out_nr = TMP / "prep_out_norefit"
_ad_nr, _man_nr = PS.run_prep(PS.parse_args([
    "--cellranger-dir", str(_pd_dir), "--sample-id", "PREP",
    "--gtf", str(_gtf_d), "--gtf-gene-id", "gene",
    "--mito-chromosome", "NC_002008.4",
    "--n-hvg", "120", "--n-pcs", "15", "--out-dir", str(_out_nr),
    "--no-doublet", "--no-pctmt-filter", "--no-refit",
]))
assert _man_nr["n_cells"] == _man["n_cells"], "cell set unexpectedly changed"
_p1 = pd.read_csv(_out_prep / "pca.tsv.gz", sep="\t", index_col=0).to_numpy()
_p2 = pd.read_csv(_out_nr / "pca.tsv.gz", sep="\t", index_col=0).to_numpy()
_d12 = float(np.abs(np.abs(_p1) - np.abs(_p2)).max())
assert _d12 > 1e-6, "PCA is identical even with --no-refit (the recompute isn't taking effect)"
print(f"  OK: PCA changes under --no-refit (max diff {_d12:.3g}) = the recompute is actually taking effect")


print("\n########## metacellcnv_scanpy v1.1: MT-based dead-cell detection (adaptive threshold, off by default) ##########")

def _plant_dead(n=2000, dead_frac=0.02, base_pct=3.0, lam=0.25, seed=0):
    """Healthy cells + planted dead cells (nuclear mRNA drops to lam-fold, mtDNA stays).

    Builds a straightforward membrane-rupture signature: mitochondria are retained so mt
    is unchanged, nuclear-derived transcripts drop, and transcript diversity (detected
    gene count) also falls.
    """
    r = np.random.default_rng(seed)
    nuc0 = r.lognormal(np.log(3000), 0.35, n)
    mt0 = r.gamma(8, base_pct / 100 / 8 * nuc0)
    dead = np.zeros(n, dtype=bool)
    if dead_frac > 0:
        dead[r.choice(n, int(n * dead_frac), replace=False)] = True
    nuc = r.poisson(np.where(dead, nuc0 * lam, nuc0)).astype(float)
    mt = r.poisson(mt0).astype(float)
    tot = nuc + mt
    ng = r.poisson(np.maximum(1200 * (1 - np.exp(-tot / 2500))
                              * np.where(dead, 0.72, 1.0), 1)).astype(float)
    return tot, mt, ng, dead

# --- (1) complexity residual ---
# Fitting a quadratic needs enough points (4 points would let the fit absorb the drop).
_tt = np.logspace(3, 4, 60)
_gg = 1200 * (1 - np.exp(-_tt / 2500))
_rr = PS.complexity_residual(_tt, _gg)
assert np.abs(_rr).max() < 0.02, np.abs(_rr).max()   # nothing left once depth is accounted for
_gg2 = _gg.copy(); _gg2[30] *= 0.5                   # halve diversity for a single cell
_rr2 = PS.complexity_residual(_tt, _gg2)
assert _rr2[30] < -0.25, _rr2[30]                    # log10(0.5) = -0.301
assert np.abs(np.delete(_rr2, 30)).max() < 0.05
print(f"  OK: complexity residual factors out what depth explains (max {np.abs(_rr).max():.3f}),"
      f" and catches the half-diversity cell at {_rr2[30]:+.2f}")

# --- (2) without overdispersion, no beta-binomial is fit / trimming the upper tail
#     narrows the null distribution ---
_r0 = np.random.default_rng(5)
_n0 = np.full(3000, 3000.0)
_k0 = _r0.binomial(3000, 0.03, 3000).astype(float)     # pure binomial (no overdispersion)
assert PS._betabinom_mom(_k0, _n0) is None, "fit a beta-binomial despite no overdispersion"
_tt3, _mt3, _ng3, _dd3 = _plant_dead(seed=11)
_ab_full = PS._betabinom_mom(_mt3, _tt3, fit_quantile=None)
_ab_trim = PS._betabinom_mom(_mt3, _tt3, fit_quantile=0.90)
assert _ab_full and _ab_trim
assert sum(_ab_trim) > sum(_ab_full), (sum(_ab_full), sum(_ab_trim))   # larger rho = narrower
print(f"  OK: falls back to binomial with no overdispersion, and trimming the top 10%"
      f" raises rho from {sum(_ab_full):.0f} to {sum(_ab_trim):.0f} (narrows the null distribution)")

# --- (3) the planted 2% dead cells are detected, and robustification is effective ---
_tt, _mt, _ng, _dd = _plant_dead()
_rep_t = PS.mito_dead_cell_report(_tt, _mt, _ng, n_mito_genes=13, fit_quantile=0.90)
_rep_f = PS.mito_dead_cell_report(_tt, _mt, _ng, n_mito_genes=13, fit_quantile=1.0)
_f = _rep_t["per_cell"].mito_flagged.values
_prec = (_f & _dd).sum() / max(_f.sum(), 1); _rec = (_f & _dd).sum() / _dd.sum()
assert _rep_t["recommend"] == "remove", _rep_t["reason"]
assert _prec >= 0.75 and _rec >= 0.75, (_prec, _rec)
assert _f.sum() > _rep_f["per_cell"].mito_flagged.values.sum(), \
    "candidate count did not increase relative to the untrimmed null distribution (robustification not effective)"
for _g in ("1_測定可能性", "2_分離可能性", "3_機構整合性", "4_非交絡"):
    assert _rep_t["gates"][_g]["passed"], (_g, _rep_t["gates"][_g]["detail"])
print(f"  OK: detected the planted dead cells at precision {_prec:.2f} / recall {_rec:.2f},"
      f" and all 4 gates pass so it recommends removal")

# --- (4) with no dead cells, the call is "do not remove" (gate 2) ---
_t0, _m0, _g0, _ = _plant_dead(dead_frac=0.0, seed=1)
_rep0 = PS.mito_dead_cell_report(_t0, _m0, _g0, n_mito_genes=13)
assert _rep0["recommend"] == "do_not_remove"
assert not _rep0["gates"]["2_分離可能性"]["passed"], _rep0["gates"]["2_分離可能性"]
print(f"  OK: with no dead cells present, {_rep0['stats']['n_flagged']} candidate cells fail gate 2")

# --- (5) when the mechanism runs backwards (mtDNA just increases), it stops (gate 3)
#     — the shape of this particular synthetic case ---
_r2 = np.random.default_rng(2); _n2 = 2000
_nuc2 = _r2.poisson(_r2.lognormal(np.log(3000), 0.35, _n2)).astype(float)
_mt2 = _r2.poisson(_nuc2 * 0.03).astype(float)
_mt2[_r2.choice(_n2, 60, replace=False)] *= 30      # cells with excess mitochondria
_tot2 = _nuc2 + _mt2
_ng2 = _r2.poisson(np.maximum(1200 * (1 - np.exp(-_tot2 / 2500)), 1)).astype(float)
_rep2 = PS.mito_dead_cell_report(_tot2, _mt2, _ng2, n_mito_genes=13)
assert _rep2["recommend"] == "do_not_remove"
_g3 = _rep2["gates"]["3_機構整合性"]
assert not _g3["passed"] and _g3["mt_ratio"] > 3.0, _g3
assert "mitochondria-rich" in _g3["detail"]
print(f"  OK: with mtDNA absolute level {_g3['mt_ratio']:.0f}x and no drop in nuclear signal,"
      " gate 3 identifies these as cells with excess mitochondria and stops")

# --- (6) if candidates are confounded with the malignant label, it stops (gate 4) ---
_t4, _m4, _g4v, _d4 = _plant_dead(seed=3)
_mal = _d4 | (np.random.default_rng(4).random(len(_d4)) < 0.05)
_rep4 = PS.mito_dead_cell_report(_t4, _m4, _g4v, n_mito_genes=13,
                                 labels={"malignant": _mal})
assert _rep4["recommend"] == "do_not_remove"
_gg4 = _rep4["gates"]["4_非交絡"]
assert not _gg4["passed"] and _gg4["tests"]["malignant"]["odds_ratio"] > 2
print(f"  OK: when candidates associate with malignant status at OR={_gg4['tests']['malignant']['odds_ratio']:.0f},"
      " gate 4 stops (can't tell whether cells are dead or just skewed)")

# --- (7) gate 1: composition test fails / insufficient counts ---
_rep_c = PS.mito_dead_cell_report(_tt, _mt, _ng, n_mito_genes=13, composition_ok=False)
assert _rep_c["recommend"] == "do_not_remove"
assert not _rep_c["gates"]["1_測定可能性"]["passed"]
assert "composition" in _rep_c["gates"]["1_測定可能性"]["detail"]
_rep_l = PS.mito_dead_cell_report(_tt, np.zeros_like(_mt), _ng, n_mito_genes=13)
assert not _rep_l["gates"]["1_測定可能性"]["passed"]
_rep_n = PS.mito_dead_cell_report(_tt, _mt, _ng, n_mito_genes=2)
assert not _rep_n["gates"]["1_測定可能性"]["passed"]
print("  OK: composition failure / insufficient mtDNA counts / too few genes all fail gate 1")

# --- (8) adaptive_mito_filter behavior by mode ---
def _pseudo_adata(tot, mt, ng, prefix):
    """A pseudo-AnnData with 13 mtDNA genes + 1 nuclear-derived gene.

    adaptive_mito_filter recounts total and mtDNA counts from counts and var['mt'],
    so split mtDNA across 13 columns to also exercise the default path (counting
    var['mt'] entries).
    """
    n = len(tot)
    X = np.zeros((n, 14))
    base = np.floor(mt / 13.0)
    for j in range(13):
        X[:, j] = base
    X[:, 0] += mt - base * 13          # remainder goes to column 1 (keeps the sum exact)
    X[:, 13] = tot - mt
    a = ad.AnnData(sparse.csr_matrix(X))
    a.obs_names = [f"{prefix}{i:05d}" for i in range(n)]
    a.var_names = [f"MT{j}" for j in range(13)] + ["NUCGENE"]
    a.var["mt"] = [True] * 13 + [False]
    a.layers["counts"] = a.X.copy()
    a.obs["n_genes_by_counts"] = ng
    a.obs["cluster_provisional"] = "0"
    return a

_adm = _pseudo_adata(_tt, _mt, _ng, "MC")
assert int(np.asarray(_adm.layers["counts"].sum(axis=1)).ravel()[0]) == int(_tt[0])

_res_rep = PS.adaptive_mito_filter(_adm, mode="report")
assert _res_rep["recommend"] == "remove", _res_rep["reason"]
assert _res_rep["n_marked_dead"] == 0, "report mode marked cells for removal"
assert _res_rep["n_flagged"] > 0
for _c in ("pct_counts_mt_adaptive", "mito_complexity_resid", "mito_qvalue",
           "mito_flagged", "mito_dead"):
    assert _c in _adm.obs, _c
assert not _adm.obs["mito_dead"].any()
print(f"  OK: mode='report' (default) records {_res_rep['n_flagged']} candidate cells"
      " but marks none for removal")

_res_auto = PS.adaptive_mito_filter(_adm, mode="auto")
assert _res_auto["n_marked_dead"] == _res_auto["n_flagged"] > 0
print(f"  OK: mode='auto' marks {_res_auto['n_marked_dead']} cells for removal once gates pass")

# On data that doesn't pass the gates, auto should still mark 0
_ad2 = _pseudo_adata(_tot2, _mt2, _ng2, "HM")
_res_auto2 = PS.adaptive_mito_filter(_ad2, mode="auto")
assert _res_auto2["recommend"] == "do_not_remove"
assert _res_auto2["n_marked_dead"] == 0, "marked cells despite failing gates"
_res_force = PS.adaptive_mito_filter(_ad2, mode="force")
assert _res_force["n_marked_dead"] > 0, "force still marked no cells"
print(f"  OK: on gate failure, auto marks 0 cells while force marks "
      f"{_res_force['n_marked_dead']} cells (with a warning)")

_res_fix = PS.adaptive_mito_filter(_ad2, mode="fixed", fixed_pct=5.0)
_pct2 = np.asarray(_ad2.obs["pct_counts_mt_adaptive"].values)
assert _res_fix["n_marked_dead"] == int((_pct2 > 5.0).sum()) > 0
print(f"  OK: mode='fixed' directly marks the {_res_fix['n_marked_dead']} cells with pctMT>5% for removal")

_res_off = PS.adaptive_mito_filter(_ad2, mode="off")
assert _res_off["recommend"] is None and _res_off["n_flagged"] == 0
print("  OK: mode='off' skips diagnostics entirely")
try:
    PS.adaptive_mito_filter(_ad2, mode="nonsense")
    raise AssertionError("an invalid mode was accepted")
except ValueError as _e:
    print("  OK: an invalid mode stops:", str(_e)[:44])

# --- (9) the CLI default is report (= no removal) ---
_amt = PS.parse_args(["--cellranger-dir", "x"])
assert _amt.mito_filter == "report", _amt.mito_filter
assert _amt.mito_filter_pct == 5.0 and _amt.mito_fit_quantile == 0.90
print("  OK: CLI default is --mito-filter report (diagnostics only, no removal)")

# --- (10) wired into run_prep: a report is written and cell count is unchanged ---
if Path("cache").exists():
    shutil.rmtree("cache")
_out_mt = TMP / "prep_out_mito"
_ad_mt, _man_mt = PS.run_prep(PS.parse_args([
    "--cellranger-dir", str(_pd_dir), "--sample-id", "PREP",
    "--gtf", str(_gtf_d), "--gtf-gene-id", "gene",
    "--mito-chromosome", "NC_002008.4",
    "--n-hvg", "120", "--n-pcs", "15", "--out-dir", str(_out_mt),
    "--no-doublet", "--no-pctmt-filter", "--mito-filter", "report",
]))
_rp = _out_mt / "mito_filter_report.json"
assert _rp.exists(), _rp
_j = _json.loads(_rp.read_text())
assert _j["mode"] == "report" and _j["n_marked_dead"] == 0
assert set(_j["gates"]) == {"1_測定可能性", "2_分離可能性", "3_機構整合性", "4_非交絡"}
assert _j["recommend"] in ("remove", "do_not_remove")
assert (_out_mt / "mito_filter_percell.tsv.gz").exists()
assert _man_mt["n_cells"] == _man["n_cells"], \
    "cell count changed under the default --mito-filter report"
print(f"  OK: run_prep writes mito_filter_report.json (recommend={_j['recommend']}),"
      " and cell count is unchanged by default")


print("\n########## scrna_common v1.0: separating out the shared components ##########")
import importlib as _il
import subprocess as _sp2

# --- (1) metacellcnv_scanpy does not load the CNV/DE pipeline or SEACells ---
_probe = _sp2.run(
    [_sys.executable, "-c",
     "import sys; import metacellcnv_scanpy; "
     "bad=[m for m in ('metacellcnv','SEACells','infercnvpy','pydeseq2') "
     "if m in sys.modules]; print('LOADED:'+','.join(bad) if bad else 'CLEAN')"],
    cwd=str(Path(__file__).resolve().parent), capture_output=True, text=True)
assert "CLEAN" in _probe.stdout, (_probe.stdout, _probe.stderr[-500:])
print("  OK: metacellcnv_scanpy does not import metacellcnv / SEACells / infercnvpy / pydeseq2 at all")

# --- (2) scrna_common itself has no CNV/DE dependency ---
_probe2 = _sp2.run(
    [_sys.executable, "-c",
     "import sys; import scrna_common; "
     "bad=[m for m in ('SEACells','infercnvpy','pydeseq2','metacellcnv') "
     "if m in sys.modules]; print('LOADED:'+','.join(bad) if bad else 'CLEAN')"],
    cwd=str(Path(__file__).resolve().parent), capture_output=True, text=True)
assert "CLEAN" in _probe2.stdout, (_probe2.stdout, _probe2.stderr[-500:])
print("  OK: scrna_common also does not import the CNV/DE modules")

# --- (3) re-exports from the pipeline still work (don't break existing call sites) ---
_REEXPORT = ["log", "warn", "inspect_gtf", "parse_gtf_gene_positions",
             "detect_mito_chromosomes", "load_chromosome_map",
             "normalize_chromosome_names", "check_exclude_chromosomes",
             "validate_dimensionality", "flag_mito_genes", "load_and_preprocess",
             "check_embedding_mito_influence", "ensure_sample_id", "mad_outlier_mask",
             "parametric_qc_filter", "detect_doublets_cluster_aware",
             "load_mito_reference_profile", "check_mito_composition",
             "check_gene_overlap", "looks_like_raw_counts", "_write_h5ad",
             "_leiden_supports_igraph",
             "MITO_CHROMOSOMES", "MITO_PREFIXES", "MITO_SYMBOLS", "RIBO_PREFIXES",
             "NUCLEAR_MITO_PREFIXES", "NUCLEAR_MITO_GENES", "N_HVG", "N_PCS",
             "CELLS_PER_METACELL", "GTF_PATH", "GTF_GENE_ID_ATTR",
             "CHROMOSOME_MAP_PATH", "MITO_REFERENCE_PROFILE"]
for _n in _REEXPORT:
    assert hasattr(P, _n), f"{_n} is missing from metacellcnv"
    assert hasattr(C, _n), f"{_n} is missing from scrna_common"
    assert getattr(P, _n) is getattr(C, _n), f"{_n} is a different object (got copied?)"
print(f"  OK: all {len(_REEXPORT)} names are re-exported from the pipeline and are the same"
      " object as in the shared module (no duplicate implementation)")

# --- (4) no duplicated definitions (nothing left behind where code was moved from) ---
_psrc = open(Path(__file__).resolve().parent / "metacellcnv.py").read()
for _n in ("def log(", "def warn(", "def inspect_gtf(", "def load_and_preprocess(",
           "def parametric_qc_filter(", "def check_mito_composition(",
           "def detect_doublets_cluster_aware("):
    assert _n not in _psrc, f"pipeline still has a definition of {_n} (duplicate)"
_csrc = open(Path(__file__).resolve().parent / "scrna_common.py").read()
for _n in ("def load_and_preprocess(", "def parametric_qc_filter(",
           "def check_mito_composition(", "def inspect_gtf("):
    assert _n in _csrc, f"scrna_common is missing {_n}"
assert "MARKER_SETS" not in _csrc, "cell-type panels are a CNV-analysis concern, keep them out of the shared module"
for _bad in ("import SEACells", "import infercnvpy", "import pydeseq2",
             "from SEACells", "_import_seacells"):
    assert _bad not in _csrc, f"scrna_common still contains '{_bad}'"
# Mentioning SEACells in a comment or warning message is fine (just explaining the
# downstream impact) — only the absence of an actual import is pinned here.
print("  OK: the moved functions leave no residue in the pipeline and exist only in the shared module")

# --- (5) the manifest records the scrna_common version ---
assert _man_mt["versions"]["scrna_common"] == C.__version__
assert "pipeline" not in _man_mt["versions"], "manifest still references the pipeline"
print(f"  OK: the manifest records scrna_common {C.__version__}")


print("\n########## metacell_annotation v1.0: metacell type calling ##########")
import metacell_annotation as MA

# --- (1) panel coverage is reported (missing markers are not silently ignored) ---
_av = MA.panel_gene_availability(["CD68", "AIF1", "TYROBP", "PECAM1"],
                                 {"Macrophage": MA.VALIDATED_PANELS["Macrophage"]})
assert _av.loc["Macrophage", "n_present"] == 3
assert "C1QA" in _av.loc["Macrophage", "missing"]
print(f"  OK: panel coverage is reported (Macrophage 3/10, missing markers listed)")

# --- (2) a threshold is rejected when the distribution isn't bimodal ---
_r = np.random.default_rng(0)
_uni = _r.normal(0, 1, 3000)
_f1 = MA.fit_bimodal_threshold(_uni)
assert not _f1["usable"], _f1
_bi = np.r_[_r.normal(-2, 0.4, 1200), _r.normal(2, 0.4, 1800)]
_f2 = MA.fit_bimodal_threshold(_bi)
assert _f2["usable"], _f2
assert abs(_f2["threshold"]) < 0.6, _f2["threshold"]
assert _f2["separation"] > 4, _f2["separation"]
# a split with a small tail should also be rejected (weight under 10%)
_tail = np.r_[_r.normal(0, 1, 3000), _r.normal(-8, 0.3, 60)]
_f3 = MA.fit_bimodal_threshold(_tail)
assert not _f3["usable"] and "weight" in _f3["reason"], _f3
print(f"  OK: unimodal is rejected / bimodal gives threshold {_f2['threshold']:+.2f}"
      f" (separation {_f2['separation']:.1f} SD) / a small-tail split is also rejected")

# --- (3) the mixture test correctly accounts for metacell size ---
# The same "10% minority" is much stronger evidence at n=150 than at n=20 cells.
def _mixq(k, n, eps=0.02):
    """Return the larger of the two null p-values (all-normal / all-tumor) as a scalar."""
    a, b = MA._mixture_pvalues(np.array([k]), np.array([n]), eps)
    return float(np.maximum(a, b)[0])

_ps, _pl = _mixq(2, 20), _mixq(15, 150)
assert _ps > _pl, (_ps, _pl)
# a pure metacell should never be called mixed
assert _mixq(0, 80) > 0.5, _mixq(0, 80)
assert _mixq(80, 80) > 0.5, _mixq(80, 80)
print(f"  OK: the binomial test reflects metacell size (10% minority at n=20 gives p={_ps:.3f},"
      f" at n=150 gives p={_pl:.1e}) / pure metacells are never called mixed")

# --- (4) end to end: can the 5 planted metacell kinds be recovered? ---
def _plant_metacells(seed=0, n_mc=120, per=60, frac_mixed=0.15):
    """Tumor (epithelial phenotype, chr1/2 gain) + 3 normal types + tumor/normal mixtures."""
    r = np.random.default_rng(seed)
    PAN, REP = MA.VALIDATED_PANELS, MA.REPORTED_PANELS
    pg = list(dict.fromkeys([g for v in PAN.values() for g in v]
                            + [g for v in REP.values() for g in v]))
    allg = pg + [f"BG{i:05d}" for i in range(4000)]
    gi = {g: i for i, g in enumerate(allg)}
    chrom = np.array([f"chr{1 + i % 8}" for i in range(len(allg))])
    base = r.lognormal(0, 1.0, size=len(allg)); base /= base.sum()
    CT = ["Macrophage", "Endothelial", "Fibroblast"]
    kinds, rows, mcs = [], [], []
    for m in range(n_mc):
        u = r.random()
        kind = "tumor" if u < 0.45 else ("mixed" if u < 0.45 + frac_mixed else CT[m % 3])
        kinds.append(kind)
        for _ in range(per):
            ck = ("tumor" if (kind == "mixed" and r.random() < 0.5)
                  else ("Macrophage" if kind == "mixed" else kind))
            p = base.copy()
            for g in (PAN["Epithelial"] if ck == "tumor" else PAN[ck]):
                p[gi[g]] = base.max() * 6
            if ck == "tumor":
                p[(chrom == "chr1") | (chrom == "chr2")] *= 1.5
            p = p / p.sum()
            rows.append(r.multinomial(r.integers(1200, 3000), p).astype(np.float32))
            mcs.append(f"MC{m:03d}")
    X = sparse.csr_matrix(np.array(rows))
    A = ad.AnnData(X)
    A.obs_names = [f"C{i:06d}" for i in range(X.shape[0])]
    A.var_names = allg; A.var["chromosome"] = chrom
    A.layers["counts"] = A.X.copy(); A.obs["metacell"] = mcs
    truth = pd.Series(kinds, index=[f"MC{m:03d}" for m in range(n_mc)])
    A.obs["grp"] = np.where(np.array([truth[m] for m in mcs]) == "tumor",
                            "malignant", "normal")
    return A, truth

_ad_mc, _truth = _plant_metacells()
_det = np.asarray((_ad_mc.layers["counts"] > 0).sum(1)).ravel()
assert np.median(_det) / _ad_mc.n_vars < 0.35, "synthetic data isn't sparse enough (should match real data)"
_ANN = MA.annotate_metacells(_ad_mc, "metacell", group_key="grp",
                             tumor_groups=("malignant", "normal"), n_perm=12, seed=0)
_exp = _truth.reindex(_ANN.index).map(
    lambda k: "Mixed:tumor+normal" if k == "mixed"
    else ("Tumor:Epithelial" if k == "tumor" else k))
_hit = (_ANN.label == _exp)
print(f"  planted composition {dict(_truth.value_counts())}")
assert _hit.all(), pd.crosstab(_truth.reindex(_ANN.index), _ANN.label).to_string()
print(f"  OK: recovered all 5 kinds at {int(_hit.sum())}/{len(_ANN)} = 100%"
      " (tumor, 3 normal types, mixed)")

# --- (5) regression: the gene-count condition should be consistent with the
#     calibrated threshold ---
# There used to be an inconsistency where the median exceeded the threshold, yet the
# test failed on "fewer than 3 genes above 2-fold".
_thr = _ANN.attrs["panel_thresholds"]
_tum = _ANN[_ANN.label == "Tumor:Epithelial"]
assert (_tum.score_Epithelial >= _thr["Epithelial"]["threshold"]).all()
assert (_tum.ngene_over_thr_Epithelial >= 3).all()
assert (_tum.Epithelial_nsig == 0).any(), \
    "want to confirm this passes even when the fixed-2x nsig is 0 (proof it counts via the calibrated threshold)"
print(f"  OK: the gene-count condition is evaluated against the calibrated threshold"
      f" (median genes above the {_thr['Epithelial']['threshold']:+.2f} threshold is "
      f"{int(_tum.ngene_over_thr_Epithelial.median())} / would be 0 at a fixed 2x)")

# --- (6) thresholds always sit above the shuffled null ---
for _p, _t in _thr.items():
    assert _t["threshold"] > _t["null_median"], (_p, _t)
    assert _t["n_null"] > 100, (_p, _t)
print("  OK: every calibrated threshold sits above the null median")

# --- (7) stops when a group is too small / coordinates are missing ---
try:
    MA.chromosome_karyotype(_ad_mc, ["a"] * _ad_mc.n_obs, "a", "b")
    raise AssertionError("passed despite an empty group")
except ValueError as _e:
    assert "size" in str(_e), str(_e)
_noc = _ad_mc.copy(); del _noc.var["chromosome"]
try:
    MA.chromosome_karyotype(_noc, _noc.obs.grp.values, "malignant", "normal")
    raise AssertionError("passed despite missing coordinates")
except ValueError as _e:
    assert "chromosome" in str(_e), str(_e)
print("  OK: stops with a stated reason when a group is too small / chromosome coordinates are missing")

# --- (8) type calling still works without karyotyping (mixture calling is skipped) ---
_A2 = MA.annotate_metacells(_ad_mc, "metacell", n_perm=8, seed=0)
assert (_A2.label_kind == "measured-mixture").sum() == 0, "called mixtures despite no karyotype"
assert (_A2.label_kind == "typed").sum() > 30, _A2.label_kind.value_counts().to_dict()
assert _A2.cnv_mixture_q.isna().all()
print(f"  OK: without karyotyping, mixture calling is skipped, and type calling still works for {int((_A2.label_kind=='typed').sum())} metacells")

# --- (9) the label vocabulary never uses "Admixture" (keeps causes distinguished) ---
_labels = set(_ANN.label) | set(_A2.label)
assert not any("Admixture" in str(l) for l in _labels), _labels
assert {"measured-mixture", "typed"} <= set(_ANN.label_kind)
print("  OK: never uses Admixture; splits into measured-mixture / typed / unclassified")

# --- (10) pipeline Step 4a-4 wiring: passing chromosome coordinates through ---
# Coordinates are written into cnv_adata.var by add_genomic_positions, and the chr-
# prefix normalization happens inside run_cnv_branch too — the single-cell-side
# adata.var never gets them. Forgetting to copy them makes chromosome_karyotype raise,
# which falls into an except clause and silently skips annotation entirely (this is
# how it was discovered on real data).
_src = Path("metacellcnv.py").read_text(encoding="utf-8")
_i = _src.find("Step 4a-4")
assert _i > 0, "could not find the Step 4a-4 block"
_blk = _src[_i:_i + 4000]
_j = _blk.find("annotate_metacells")
assert _j > 0, "could not find the annotate_metacells call"
_before = _blk[:_j]
assert 'cnv_adata.var["chromosome"]' in _before, (
    "Step 4a-4 does not copy cnv_adata.var['chromosome'] to the single-cell side")
assert 'adata.var["chromosome"] =' in _before, (
    "Step 4a-4 does not set adata.var['chromosome']")
assert "chromosome_karyotype" in _before
assert 'if "chromosome" in cnv_adata.var.columns' in _before, (
    "no branch to skip karyotyping gracefully when coordinates are missing (would skip everything via an exception)")
print("  OK: Step 4a-4 copies cnv_adata's coordinates to the single-cell side before computing the karyotype")

print("\n\n===== ALL SMOKE TESTS PASSED =====")

print("\n########## v3.3: externalized marker tables (dog / human / mouse) ##########")
# --- (1) all 3 species load, and their structure matches ---
_tabs = {sp: C.load_marker_table(sp) for sp in C.BUILTIN_MARKER_SPECIES}
_types = {sp: set(t.df.cell_type) for sp, t in _tabs.items()}
assert _types["dog"] == _types["human"] == _types["mouse"], _types
assert len(_types["dog"]) == 11, sorted(_types["dog"])
print(f"  OK: all {len(C.BUILTIN_MARKER_SPECIES)} species share the same {len(_types['dog'])} types")

# --- (2) mouse naming convention is converted, and orthologless genes are dropped ---
_mg = set(_tabs["mouse"].df.gene)
_hg = set(_tabs["human"].df.gene)
assert "Cd3d" in _mg and "CD3D" not in _mg, sorted(_mg)[:5]
assert "GNLY" in _hg and not any(g.lower() == "gnly" for g in _mg), "GNLY is still present"
assert "Lyz2" in _mg and "Fcgr3" in _mg, "the LYZ/FCGR3A special-case conversion isn't working"
assert len(_mg) < len(_hg), (len(_mg), len(_hg))
print(f"  OK: mouse is Title-cased with orthologless genes excluded (human {len(_hg)} -> mouse {len(_mg)} genes)")

# --- (3) only the 4 validated types are used for labels; types that can turn malignant
#     never enter the normal reference ---
for sp, t in _tabs.items():
    assert set(t.label_panels) == {"Macrophage", "Endothelial", "Fibroblast", "Epithelial"}, \
        (sp, list(t.label_panels))
    assert "Fibroblast" not in t.normal_celltypes and "Epithelial" not in t.normal_celltypes, \
        (sp, t.normal_celltypes)
print("  OK: 4 label types, with Fibroblast/Epithelial excluded from the normal reference")

# --- (4) the pipeline default can be swapped ---
_before = list(P.MARKER_SETS["T/NK"])
P.apply_marker_table("mouse")
assert P.MARKER_SETS["T/NK"][0] == "Cd3d", P.MARKER_SETS["T/NK"][:3]
P.apply_marker_table("dog")
assert P.MARKER_SETS["T/NK"] == _before, "did not restore correctly"
print("  OK: --markers can swap out MARKER_SETS / KNOWN_NORMAL_CELLTYPES")

# --- (5) an unknown species stops and lists the valid options ---
try:
    C.load_marker_table("chicken")
    raise AssertionError("an unknown species was accepted")
except FileNotFoundError as _e:
    assert "dog" in str(_e) and "human" in str(_e), str(_e)
print("  OK: an unknown species stops and lists the built-in options")

# --- (6) a CSV with missing columns stops ---
_bad = TMP / "bad_markers.csv"
_bad.write_text("cell_type,gene\nT/NK,CD3D\n", encoding="utf-8")
try:
    C.load_marker_table(str(_bad))
    raise AssertionError("a CSV with missing columns was accepted")
except ValueError as _e:
    assert "missing required columns" in str(_e), str(_e)
print("  OK: a CSV missing required columns stops and lists them")

print("\n########## v3.5: CNV-based tumor lineage (Dollo parsimony) ##########")
import cnv_lineage as LIN

def _make_lineage_profiles(planted: bool, seed: int = 0):
    """planted=True: a known lineage root->A->B and root->C. False: events scattered at random."""
    r = np.random.default_rng(seed)
    n_bins, n_norm = 40, 60
    rows, names, norm = [], [], []
    for i in range(n_norm):
        rows.append(np.zeros(n_bins)); names.append(f"N{i}"); norm.append(True)
    if planted:
        for cl, ev in {"A": [0, 1, 2], "B": [0, 1, 2, 3, 4], "C": [5, 6]}.items():
            for i in range(25):
                v = np.zeros(n_bins); v[ev] = -1.0
                rows.append(v); names.append(f"{cl}{i}"); norm.append(False)
    else:
        for i in range(75):
            v = np.zeros(n_bins); v[r.choice(n_bins, 3, replace=False)] = -1.0
            rows.append(v); names.append(f"R{i}"); norm.append(False)
    X = np.vstack(rows) + r.normal(0, 0.12, (len(rows), n_bins))
    cols = [f"chr{j//4+1}:{j%4}" for j in range(n_bins)]
    return pd.DataFrame(X, index=names, columns=cols), np.array(norm)

# --- (1) a planted lineage is detected, and randomly scattered data is not ---
_pf, _nm = _make_lineage_profiles(True, 0)
_r1 = LIN.build_lineage(_pf, _nm, n_perm=80, seed=0)
_pf0, _nm0 = _make_lineage_profiles(False, 1)
_r0 = LIN.build_lineage(_pf0, _nm0, n_perm=80, seed=0)
assert _r1["has_structure"], _r1["null_test"]
assert not _r0["has_structure"], _r0["null_test"]
print(f"  OK: the planted lineage is detected at p={_r1['null_test']['p_value']:.4f},"
      f" and random data is not detected at p={_r0['null_test']['p_value']:.4f}")

# --- (2) the minimum support count is set from the noise floor (a fixed value alone
#     lets spurious events through) ---
_B, _ei = LIN.cnv_events(_pf, _nm, min_support=3)
assert _ei["support_threshold"] > 3, _ei
assert _B.shape[1] <= 12, f"too many spurious events remain: {_B.shape[1]}"
print(f"  OK: the support floor is set from noise at {_ei['support_threshold']}"
      f" (not the fixed 3) / {_B.shape[1]} events adopted")

# --- (3) contaminating the normal reference with tumor cells doesn't change which
#     events get called (order-statistic threshold) ---
# Turn 2 of the 60 normal cells (3.3%) into "actually tumor" and set them to -1 in a
# real deletion region. A quantile-based threshold would get dragged by the
# contamination and miss the real event; an order statistic tolerates contamination
# up to max_normal_rate by construction.
_pf3 = _pf.copy()
_pf3.iloc[np.where(_nm)[0][:2], 0:3] = -1.0
_Bc, _ = LIN.cnv_events(_pf, _nm)
_B3, _ei3 = LIN.cnv_events(_pf3, _nm)
assert list(_B3.columns) == list(_Bc.columns), (list(_Bc.columns), list(_B3.columns))
assert _ei3["max_normal_carriers"] >= 2, _ei3
print(f"  OK: called events stay the same even with 3.3% tumor contamination in the reference"
      f" (contamination tolerance {_ei3['max_normal_carriers']})")

# --- (4) the four-gamete condition holds for a perfectly nested case ---
_perfect = pd.DataFrame({
    "loss:a": [1, 1, 1, 1, 0, 0],
    "loss:b": [1, 1, 0, 0, 0, 0],
    "loss:c": [0, 0, 0, 0, 1, 1],
}, index=[f"m{i}" for i in range(6)]).astype(bool)
assert LIN.four_gamete_violations(_perfect, min_count=1)["n_violating_pairs"] == 0
_conflict = _perfect.copy(); _conflict["loss:d"] = [1, 0, 1, 0, 1, 0]
assert LIN.four_gamete_violations(_conflict, min_count=1)["n_violating_pairs"] > 0
print("  OK: 0 violations for nested traits; adding a crossing trait produces violations")

# --- (5) the parsimony score never improves by dropping events (no loophole) ---
_root, _st = LIN.dollo_tree(_perfect, min_count=1, verbose=False)
assert _st["n_events_dropped"] == 0
_rootc, _stc = LIN.dollo_tree(_conflict, min_count=1, verbose=False)
assert _stc["n_events_dropped"] >= 1
assert _stc["parsimony_score"] > _st["parsimony_score"], (_stc, _st)
print(f"  OK: dropping events makes the tree's score worse"
      f" ({_st['parsimony_score']} -> {_stc['parsimony_score']})")

# --- (6) Newick output is well-formed, and the branch table comes out ---
_nwk = LIN.to_newick(_root)
assert _nwk.endswith(";") and _nwk.count("(") == _nwk.count(")"), _nwk
assert all(m in _nwk for m in _perfect.index), _nwk
_tab = LIN.branch_table(_root)
assert {"node", "parent", "events", "n_metacells_subtree"} <= set(_tab.columns)
print(f"  OK: Newick has balanced parentheses and includes every metacell / branch table has {len(_tab)} rows")

# --- (7) the comparison neighbor-joining method also works ---
_njr, _nji = LIN.nj_tree(_perfect)
assert _nji["n_leaves"] == 6 and LIN.to_newick(_njr).endswith(";")
print("  OK: the comparison neighbor-joining method also returns a tree")

print("\n########## v3.6: speedups (incremental updates, vectorization) don't change results ##########")
# --- (1) the incremental update to _update_B matches a naive recomputation exactly ---
_rng6 = np.random.default_rng(3)
_X6 = np.vstack([_rng6.normal(c, 1.0, (400, 20)) for c in (0, 4, 8, 12)])
_ad6 = ad.AnnData(np.zeros((1600, 5), dtype=np.float32))
_ad6.obsm["X_pca"] = _X6
_ad6.obs_names = [f"c{i}" for i in range(1600)]
_M6 = MC.adaptive_rbf_kernel(_ad6, k=15)
_n6, _k6 = 1600, 40
_B6 = np.zeros((_n6, _k6))
_B6[_rng6.choice(_n6, _k6, replace=False), np.arange(_k6)] = 1.0
_A6 = _rng6.random((_k6, _n6)); _A6 /= _A6.sum(0)
_m6 = MC.MetacellModel(n_metacells=_k6, verbose=False); _m6.set_kernel(_M6)
_Binc = _m6._update_B(_A6, _B6.copy())

def _naive_update_B(M, A, B, fw):
    t1 = A @ A.T; t2 = MC._KX(M, A.T); cols = np.arange(B.shape[1])
    for t in range(fw):
        G = 2.0 * (MC._KX(M, B @ t1) - t2)
        am = np.argmin(G, 0)
        D = -B.copy(); D[am, cols] += 1.0
        B = B + 2.0 / (t + 2.0) * D
    return B

_Bnai = _naive_update_B(_M6, _A6, _B6.copy(), _m6.fw_iters)
_d6 = float(np.abs(np.asarray(_Binc) - _Bnai).max())
assert _d6 < 1e-10, _d6
print(f"  OK: incremental _update_B matches the naive version (max diff {_d6:.1e})")

# --- (2) the vectorized overdispersion estimator matches the per-bin loop ---
_nc6, _nb6 = 250, 30
_W6 = _rng6.poisson(40, (_nc6, _nb6)).astype(float)
_n_i6 = _W6.sum(1)
_p6 = _W6.sum(0) / _W6.sum()
_loop = np.array([CV._betabinom_phi(_W6[:, j], _n_i6, _p6[j]) for j in range(_nb6)])
_vec = CV._betabinom_phi_vec(_W6, _n_i6, _p6)
assert np.allclose(_loop, _vec, atol=1e-12), np.abs(_loop - _vec).max()
_mu6 = np.outer(_n_i6, _p6)
_la = np.array([CV._nb_alpha(_W6[:, j], _mu6[:, j]) for j in range(_nb6)])
_va = CV._nb_alpha_vec(_W6, _mu6)
assert np.allclose(_la, _va, atol=1e-12), np.abs(_la - _va).max()
print(f"  OK: vectorized phi / alpha match the loop version"
      f" (max diff {max(np.abs(_loop-_vec).max(), np.abs(_la-_va).max()):.1e})")

# --- (3) with the same seed, fit reproduces exactly ---
_a6a = _ad6.copy(); _a6b = _ad6.copy()
MC.MetacellModel(n_metacells=30, seed=7, verbose=False).fit(_a6a, max_iter=12, min_iter=4)
MC.MetacellModel(n_metacells=30, seed=7, verbose=False).fit(_a6b, max_iter=12, min_iter=4)
assert (_a6a.obs["SEACell"].values == _a6b.obs["SEACell"].values).all()
print("  OK: running twice with the same seed produces identical assignments")

print("\n########## v3.7: plotly figures (English labeling) ##########")
import metacellcnv_plotly as PL

# --- (1) group labels prefer anno_label, and never fall back to raw cell_type ---
_obs7 = pd.DataFrame({
    "cell_type": ["Myeloid_0", "Myeloid_1", "Epithelial", "Fibroblast"],
    "anno_label": ["Unclassified:low-depth", "Macrophage",
                   "Tumor:Epithelial", "Fibroblast"],
    "putative_malignant": ["normal", "normal", "malignant", "normal"],
    "n_cells": [30, 40, 25, 35], "cnv_score": [0.1, 0.2, 1.5, 0.3],
}, index=[f"SEACell-{i}" for i in range(4)])
_lab, _used = PL.resolve_group_labels(_obs7)
assert _used == "anno_label", _used
assert not any("Myeloid_" in v for v in _lab), list(_lab)
_pretty = PL.prettify_labels(_lab)
assert "Tumour Epithelial" in set(_pretty), list(_pretty)
assert "Unclassified (low depth)" in set(_pretty), list(_pretty)
print(f"  OK: group labels use {_used}, and Myeloid_* never appears in the figure")

# --- (2) without anno_label, it warns and falls back to cell_type ---
_obs7b = _obs7.drop(columns=["anno_label"])
_lab2, _used2 = PL.resolve_group_labels(_obs7b)
assert _used2 == "cell_type"
print("  OK: without anno_label, it warns and falls back to cell_type")

# --- (3) figures are in English, with no Japanese text mixed in ---
_fig_q = PL.fig_metacell_quality(_obs7)
_txt = _fig_q.to_plotly_json()
import json as _json
_s7 = _json.dumps(_txt, ensure_ascii=False)
_ja = [ch for ch in _s7 if "぀" <= ch <= "ヿ" or "一" <= ch <= "鿿"]
assert not _ja, f"figure contains Japanese text: {''.join(_ja[:20])}"
assert "Metacell quality" in _s7
print("  OK: no Japanese text mixed into figure strings")

# --- (4) the CNV heatmap drops non-numeric columns and still renders ---
_mat7 = pd.DataFrame({"chr1": [0.1, -0.2, 0.3, 0.0], "chr2": [-0.1, 0.2, -0.3, 0.1],
                      "cell_type": _obs7["cell_type"].values}, index=_obs7.index)
_fig_h = PL.fig_cnv_heatmap(_mat7, _obs7)
_z = _fig_h.data[0].z
assert _z.shape == (4, 2), _z.shape
print("  OK: renders even from a CNV table with an annotation column mixed in")

# --- (5) the lineage-tree figure bakes in "no structure" when applicable ---
_evt = pd.DataFrame({"loss:chrA": [1, 1, 0, 0], "loss:chrB": [1, 0, 0, 0]},
                    index=_obs7.index).astype(bool)
_root7, _ = LIN.dollo_tree(_evt, min_count=1, verbose=False)
_fig_t = PL.fig_lineage_tree(_root7, has_structure=False)
_st = _json.dumps(_fig_t.to_plotly_json(), ensure_ascii=False)
assert "not supported" in _st.lower(), "the no-structure warning is missing from the figure"
_fig_t2 = PL.fig_lineage_tree(_root7, has_structure=True)
assert "not supported" not in _json.dumps(_fig_t2.to_plotly_json()).lower()
print("  OK: the figure bakes in a warning whenever structure isn't supported")

# --- (6) everything combines into a single HTML file ---
_rep = PL.write_report({"Quality": _fig_q, "CNV": _fig_h}, TMP / "rep.html",
                       title="test")
_html = _rep.read_text(encoding="utf-8")
assert _html.count("plotly-graph-div") == 2, _html.count("plotly-graph-div")
assert "<h2>Quality</h2>" in _html and "<h2>CNV</h2>" in _html
print("  OK: the combined report comes out as a single HTML file")
