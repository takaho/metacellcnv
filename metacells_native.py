#!/usr/bin/env python3
"""metacells_native.py — Native metacell (archetype) construction, without a
SEACells dependency.

This module reimplements the SEACells archetypal-analysis algorithm
(adaptive-bandwidth RBF kernel + Frank-Wolfe factorization) natively. Doing
so avoids runtime monkey-patching of SEACells internals (fragile across
SEACells versions), memory blow-up from materializing a dense n x n kernel
matrix on large datasets, and non-determinism from SEACells' unseeded RNG
calls.

Key formulas (identical to the reference implementation):
  Adaptive-bandwidth RBF kernel M (symmetric, sparse):
      M_ij = exp(-||x_i - x_j||^2 / (sigma_i sigma_j))   for kNN edges (i,j)
      sigma_i = median neighbor distance for cell i (the k//2-th smallest)
  Archetypal decomposition:
      min_{A,B} ||M - M B A||_F^2
      B: n x k, columns on the simplex (metacells as convex combinations of cells)
      A: k x n, columns on the simplex (cells as convex combinations of metacells)
  Both are solved via Frank-Wolfe (no projection needed, since the gradient's
  argmin lands on a simplex vertex). Step size is either the fixed schedule
  2/(t+2) (default, matches the reference implementation) or exact line search.

Important: driving the objective further down degrades metacell quality.
The objective ||M - M B A||^2 favors solutions where one archetype exactly
reproduces a single cell and absorbs most other cells into one giant
metacell. The reference implementation stays usable only because it stops
optimization early (fixed step size, 50 inner iterations) — so `fw_iters` is
effectively a regularization parameter, not just a speed knob. This module
therefore defaults to line_search=False, fw_iters=50 (matching that
behavior) and always checks the size-distribution Gini via `size_balance()`,
warning when the result looks degenerate.

See README.md for a summary and validation results.
Greedy column subset selection reference: arXiv:1312.6838.
"""

from __future__ import annotations

import time
from typing import Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    from scrna_common import log, warn
except Exception:  # pragma: no cover - when used standalone
    def log(msg: str) -> None:
        print(f"[metacells] {msg}", flush=True)

    def warn(msg: str) -> None:
        print(f"[metacells][WARNING] {msg}", flush=True)


__version__ = "1.0"

DEFAULT_N_NEIGHBORS = 15
DEFAULT_FW_ITERS = 50
DEFAULT_CONVERGENCE_EPSILON = 1e-3
# Degeneracy check: how many multiples of the "even split" share 1/k the
# largest metacell is allowed to occupy. Empirically, healthy solutions were
# 1.8-3.3x while degenerate ones were 33-77x (reference SEACells: 3.3%,
# line search: 57.8%, with k=100 so an even split is 1%).
# Gini alone isn't sufficient: it can reach 0.47 even for healthy solutions.
MAX_SIZE_SHARE_RATIO = 10.0
MAX_SIZE_SHARE_ABS = 0.10
# Column block width for computing f_i = sum_j K_ij^2 without materializing K.
CSSP_BLOCK = 1024
# Lower bound on sigma_i, just to avoid division by zero; inert on real data.
MIN_BANDWIDTH = 1e-12


# ---------------------------------------------------------------------------
# Kernel construction
# ---------------------------------------------------------------------------

def _kth_smallest_nonzero(dist_csr: sp.csr_matrix, kth: int) -> np.ndarray:
    """Return the `kth` smallest (1-indexed) nonzero distance in each row.

    Computes the same quantity as SEACells' `kth_neighbor_distance`. Rows
    with fewer than `kth` nonzero entries would make the reference
    implementation return a 0 norm (division by zero in the RBF kernel);
    here we fall back to the row's max distance instead and warn.
    """
    n = dist_csr.shape[0]
    out = np.zeros(n, dtype=float)
    n_short = 0
    indptr, data = dist_csr.indptr, dist_csr.data
    for i in range(n):
        row = data[indptr[i]:indptr[i + 1]]
        row = row[row > 0]
        if row.size == 0:
            n_short += 1
            continue
        if row.size < kth:
            n_short += 1
            out[i] = row.max()
        else:
            out[i] = np.partition(row, kth - 1)[kth - 1]
    if n_short:
        warn(f"{n_short} cells had fewer than {kth} neighbors."
             " Using each row's max distance as bandwidth (SEACells would divide by zero here)")
    return np.maximum(out, MIN_BANDWIDTH)


def adaptive_rbf_kernel(adata, use_rep: str = "X_pca", k: int = DEFAULT_N_NEIGHBORS,
                        graph_construction: Literal["union", "intersection"] = "union",
                        block: int = 200_000,
                        reuse_distances: bool = True) -> sp.csr_matrix:
    """Build the adaptive-bandwidth RBF kernel M (same definition as SEACells' build_graph.rbf).

    The reference implementation computes distances from each row to all
    cells before masking by the kNN graph (O(n^2 d)). This computes
    distances only at the sparse graph's nonzero positions (O(nnz * d));
    the resulting values are identical.
    """
    import scanpy as sc

    if use_rep not in adata.obsm:
        raise KeyError(f"obsm['{use_rep}'] not found. Run PCA first")
    X = np.asarray(adata.obsm[use_rep], dtype=np.float64)
    n = X.shape[0]

    prev = (adata.uns.get("neighbors", {}) or {}).get("params", {}) or {}
    reusable = (reuse_distances and "distances" in adata.obsp
                and int(prev.get("n_neighbors", -1)) == int(k)
                and str(prev.get("use_rep", "")) == str(use_rep))
    if reusable:
        log(f"Reusing existing obsp['distances']"
            f" (n_neighbors={k}, use_rep={use_rep})")
    else:
        # scanpy's neighbor search falls back to pynndescent (approximate,
        # stochastic) for large n, so results vary slightly run to run even
        # on the same data. Reusing a cached graph when possible is more
        # reproducible.
        sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=k, knn=True)
    dist = adata.obsp["distances"].tocsr()

    sigma = _kth_smallest_nonzero(dist, max(1, k // 2))

    knn = dist.copy()
    knn.data[:] = 1.0
    knn = knn.tocsr()
    knn.setdiag(1.0)
    knn.eliminate_zeros()
    if graph_construction == "union":
        graph = ((knn + knn.T) > 0).astype(float)
    elif graph_construction in ("intersect", "intersection"):
        graph = knn.multiply(knn.T)
    else:
        raise ValueError("graph_construction must be 'union' or 'intersection'")
    graph = sp.csr_matrix(graph)

    coo = graph.tocoo()
    rows, cols = coo.row, coo.col
    vals = np.empty(rows.size, dtype=np.float64)
    for s in range(0, rows.size, block):
        e = min(s + block, rows.size)
        diff = X[rows[s:e]] - X[cols[s:e]]
        d2 = np.einsum("ij,ij->i", diff, diff)
        vals[s:e] = np.exp(-d2 / (sigma[rows[s:e]] * sigma[cols[s:e]]))
    M = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    M.eliminate_zeros()
    log(f"RBF kernel: {n:,} x {n:,}, nnz {M.nnz:,}"
        f" ({M.nnz / n:.1f} per row), k={k}, {graph_construction}")
    return M


# ---------------------------------------------------------------------------
# Factorized matrix products (never materialize K)
# ---------------------------------------------------------------------------

def _KX(M: sp.csr_matrix, Xd: np.ndarray) -> np.ndarray:
    """Evaluate K @ Xd without materializing K = M Mᵀ = M M (M is symmetric)."""
    return M @ (M @ Xd)


def _rss(M: sp.csr_matrix, A: np.ndarray, B: np.ndarray,
         m_fro_sq: float) -> float:
    """Compute ‖M − M B A‖_F exactly, without ever forming an n×n matrix.

    ‖M − MBA‖² = ‖M‖² − 2⟨M, MBA⟩ + ‖MBA‖²
      ⟨M, MBA⟩ = tr(Mᵀ M B A) = tr(K B A) = Σ (K B) ∘ Aᵀ        … only n×k
      ‖MBA‖²   = tr(Aᵀ Bᵀ K B A) = Σ (Bᵀ K B) ∘ (A Aᵀ)          … only k×k
    Mathematically equivalent to `kernel_matrix - kernel_matrix.dot(B).dot(A)`,
    which would require 80 GB dense at n=100,000.
    """
    C = _KX(M, B)                      # n×k
    cross = float(np.sum(C * A.T))
    P = B.T @ C                        # k×k
    quad = float(np.sum(P * (A @ A.T)))
    val = m_fro_sq - 2.0 * cross + quad
    return float(np.sqrt(max(val, 0.0)))


# ---------------------------------------------------------------------------
# Initialization
# ---------------------------------------------------------------------------

def init_greedy_cssp(M: sp.csr_matrix, k: int, block: int = CSSP_BLOCK,
                     verbose: bool = True) -> np.ndarray:
    """Select k cells by greedy adaptive column subset selection (arXiv:1312.6838).

    Maintains f_i = Σ_j K_ij² and g_i = K_ii for K = M Mᵀ (symmetric),
    greedily picking the column that maximizes score = f/g and deflating the
    residual. Same update rule as the reference implementation, but K is
    never materialized: f is computed in column blocks, and K's columns are
    obtained via matrix-vector products.
    """
    n = M.shape[0]
    if k >= n:
        raise ValueError(f"n_metacells {k} must be less than the number of cells {n}")

    # g_i = K_ii = Σ_j M_ij² (M symmetric)
    g = np.asarray(M.multiply(M).sum(axis=1)).ravel().astype(float)
    # f_i = Σ_j K_ij² computed in column blocks, staying sparse
    # (densifying would need an n x block array — 450 MB at n=27,411, block=2,048)
    f = np.empty(n, dtype=float)
    for s in range(0, n, block):
        e = min(s + block, n)
        Kb = (M @ M[:, s:e]).tocsc()          # n x (e-s), stays sparse
        f[s:e] = np.asarray(Kb.multiply(Kb).sum(axis=0)).ravel()
        del Kb

    omega = np.zeros((k, n), dtype=float)
    centers = np.zeros(k, dtype=int)
    g_safe = np.maximum(g, 1e-12)

    for j in range(k):
        p = int(np.argmax(f / g_safe))
        col = np.asarray((M @ M[:, p].toarray()).ravel(), dtype=float)  # K[:, p]
        delta = col - (omega[:, p][:, None] * omega).sum(axis=0)
        delta[p] = max(0.0, delta[p])
        o = delta / max(np.sqrt(delta[p]), 1e-6)
        oo = o * o
        term1 = float(np.dot(o, o)) * oo
        pl = np.zeros(n, dtype=float)
        if j:
            pl = omega[:j].T @ (omega[:j] @ o)
        term2 = o * (_KX(M, o.reshape(-1, 1)).ravel() - pl)
        f += -2.0 * term2 + term1
        g_safe = np.maximum(g + oo, 1e-12)
        g = g + oo
        omega[j, :] = o
        centers[j] = p

    if verbose:
        log(f"Greedy CSSP init: {len(np.unique(centers))}/{k} distinct")
    return centers


def init_maxmin(X: np.ndarray, k: int, seed: int = 0) -> np.ndarray:
    """Farthest-point (max-min) sampling; covers rare states well.

    Serves the same purpose as SEACells' waypoint initialization (max-min
    sampling on palantir diffusion components), but runs directly in PCA
    space without a palantir dependency.
    """
    n = X.shape[0]
    rng = np.random.default_rng(seed)
    first = int(rng.integers(n))
    picked = [first]
    d2 = np.einsum("ij,ij->i", X - X[first], X - X[first])
    for _ in range(1, k):
        nxt = int(np.argmax(d2))
        picked.append(nxt)
        dn = X - X[nxt]
        d2 = np.minimum(d2, np.einsum("ij,ij->i", dn, dn))
    return np.array(picked, dtype=int)


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

class MetacellModel:
    """Metacell construction via archetypal decomposition (no SEACells dependency).

    Usage mirrors SEACells:
        model = MetacellModel(n_metacells=300, use_rep="X_pca")
        model.fit(adata)                      # writes adata.obs['SEACell']
        model.hard_assignments()               # pd.DataFrame with column 'SEACell'
    """

    def __init__(self, n_metacells: int, use_rep: str = "X_pca",
                 n_neighbors: int = DEFAULT_N_NEIGHBORS,
                 fw_iters: int = DEFAULT_FW_ITERS,
                 line_search: bool = False,
                 convergence_epsilon: float = DEFAULT_CONVERGENCE_EPSILON,
                 l2_penalty: float = 0.0,
                 init: Literal["cssp", "maxmin", "mix"] = "cssp",
                 maxmin_fraction: float = 0.5,
                 float32: bool = False,
                 seed: int = 0,
                 graph_construction: str = "union",
                 verbose: bool = True):
        if int(n_metacells) < 2:
            raise ValueError("n_metacells must be at least 2")
        self.k = int(n_metacells)
        self.use_rep = use_rep
        self.n_neighbors = int(n_neighbors)
        self.fw_iters = int(fw_iters)
        self.line_search = bool(line_search)
        self.convergence_epsilon = float(convergence_epsilon)
        self.l2_penalty = float(l2_penalty)
        self.init = init
        self.maxmin_fraction = float(maxmin_fraction)
        self.float32 = bool(float32)
        self.seed = int(seed)
        self.graph_construction = graph_construction
        self.verbose = bool(verbose)

        self.M_: sp.csr_matrix | None = None
        self._M_csc = None
        self.A_: np.ndarray | None = None
        self.B_: np.ndarray | None = None
        self.archetypes_: np.ndarray | None = None
        self.RSS_iters: list[float] = []
        self.n_iter_ = 0
        self.converged_ = False
        self._m_fro_sq = None
        self._adata = None

    # --- Kernel -----------------------------------------------------------
    def build_kernel(self, adata) -> sp.csr_matrix:
        M = adaptive_rbf_kernel(adata, use_rep=self.use_rep, k=self.n_neighbors,
                                graph_construction=self.graph_construction)
        if self.float32:
            M = M.astype(np.float32)
        self.M_ = M
        self._M_csc = None
        self._m_fro_sq = float(M.multiply(M).sum())
        return M

    def set_kernel(self, M: sp.csr_matrix) -> None:
        """Use a precomputed kernel (e.g. to feed in SEACells' kernel for validation)."""
        M = sp.csr_matrix(M)
        self.M_ = M.astype(np.float32) if self.float32 else M
        self._M_csc = None
        self._m_fro_sq = float(self.M_.multiply(self.M_).sum())

    # --- Initialization -----------------------------------------------------
    def initialize(self, adata=None, initial_archetypes=None) -> None:
        M = self.M_
        if M is None:
            raise RuntimeError("Call build_kernel() or set_kernel() first")
        n = M.shape[0]
        dtype = np.float32 if self.float32 else np.float64

        if initial_archetypes is not None:
            arch = np.asarray(initial_archetypes, dtype=int)
        elif self.init == "cssp":
            arch = init_greedy_cssp(M, self.k, verbose=self.verbose)
        else:
            if adata is None:
                raise ValueError("adata is required for maxmin initialization")
            X = np.asarray(adata.obsm[self.use_rep], dtype=np.float64)
            if self.init == "maxmin":
                arch = init_maxmin(X, self.k, seed=self.seed)
            else:  # mix
                n_mm = max(1, int(round(self.k * self.maxmin_fraction)))
                mm = init_maxmin(X, n_mm, seed=self.seed)
                gr = init_greedy_cssp(M, self.k, verbose=False)
                cand = np.hstack([mm, gr])
                uniq, idx = np.unique(cand, return_index=True)
                arch = uniq[np.argsort(idx)][:self.k]
                if self.verbose:
                    log(f"Mixed init: {n_mm} max-min + greedy CSSP for the remainder")

        uniq, idx = np.unique(arch, return_index=True)
        arch = uniq[np.argsort(idx)]
        if len(arch) < self.k:
            # Fill any duplicates with greedy CSSP candidates (same policy as
            # the reference implementation)
            extra = init_greedy_cssp(M, min(n - 1, self.k * 2), verbose=False)
            for c in extra:
                if len(arch) >= self.k:
                    break
                if c not in arch:
                    arch = np.append(arch, c)
        self.archetypes_ = arch[:self.k]

        B = np.zeros((n, self.k), dtype=dtype)
        B[self.archetypes_, np.arange(self.k)] = 1.0

        rng = np.random.default_rng(self.seed)
        A = rng.random((self.k, n)).astype(dtype)
        A /= A.sum(0)
        A = self._update_A(B, A)

        self.A_, self.B_ = A, B
        self.RSS_iters = [_rss(M, A, B, self._m_fro_sq)]
        if self.verbose:
            log(f"Initial RSS {self.RSS_iters[0]:.4f}"
                f" / convergence threshold {self.convergence_epsilon * self.RSS_iters[0]:.5f}")

    # --- Frank-Wolfe --------------------------------------------------------
    def _update_A(self, B: np.ndarray, A: np.ndarray) -> np.ndarray:
        """Update A. Same gradient and vertex selection as the reference implementation.

        Factorized: t2 = (K B)ᵀ = (M (M B))ᵀ. The reference implementation
        multiplies the n×n sparse K by the dense n×k matrix on every
        iteration; applying M twice measured 12.2x faster (M has 25-60
        nonzeros per row vs. K's 994).
        """
        M = self.M_
        k, n = A.shape
        t2 = _KX(M, B).T                    # k×n
        t1 = t2 @ B                         # k×k
        cols = np.arange(n)
        for t in range(self.fw_iters):
            G = 2.0 * (t1 @ A - t2) - self.l2_penalty * A
            amins = np.argmin(G, axis=0)
            D = -A.copy()
            D[amins, cols] += 1.0           # D = e - A
            if self.line_search:
                num = -0.5 * float(np.sum(D * G))
                den = float(np.sum(D * (t1 @ D)))
                gamma = 0.0 if den <= 0 else min(1.0, max(0.0, num / den))
                if gamma <= 0.0:
                    break
            else:
                gamma = 2.0 / (t + 2.0)
            A = A + gamma * D
        return A

    def _update_B(self, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        """Update B. Same gradient and vertex selection as the reference
        implementation, but performs an **incremental update**.

        Recomputing G = 2 (K B t1 - K Aᵀ), K B t1 = M(M(B t1)) from scratch
        every inner iteration measured 380 ms/iteration at n=27,411, k=365 —
        19 s for 50 FW iterations, the dominant cost of metacell construction.

        Since B is updated as B <- B + gamma(E - B), where E has exactly one
        nonzero per column, D = K B t1 can instead be updated incrementally:
            D <- (1 - gamma) D + gamma (K E) t1
        K E = M (M E) is sparse (~994 nonzeros/column), so (K E) t1 is a
        sparse-times-dense product, bringing this down to ~20 ms/iteration.
        Mathematically identical; measured max absolute difference ~1e-10.
        """
        M = self.M_
        Mc = self._Mcsc()
        n, k = B.shape
        t1 = A @ A.T                        # k×k (constant across B's inner iterations)
        t2 = _KX(M, A.T)                    # n×k
        D = _KX(M, B @ t1)                  # n×k, maintained incrementally from here
        G = np.empty_like(D)
        cols = np.arange(k)
        B = np.array(B, copy=True, order="C")
        for t in range(self.fw_iters):
            np.subtract(D, t2, out=G)       # G/2; argmin is invariant to a constant factor
            amins = np.argmin(G, axis=0)
            if self.line_search:
                Dm = -B.copy()
                Dm[amins, cols] += 1.0
                num = -float(np.sum(Dm * G))
                den = float(np.sum((Dm.T @ _KX(M, Dm)) * t1))
                gamma = 0.0 if den <= 0 else min(1.0, max(0.0, num / den))
                if gamma <= 0.0:
                    break
            else:
                gamma = 2.0 / (t + 2.0)
            # K E is "K's amins columns"; selecting columns is faster than
            # building E and taking a sparse-sparse product
            KE = M @ Mc[:, amins]           # sparse (measured: 3.6% the nonzeros of dense n x k)
            KEt1 = KE @ t1                  # n×k dense
            B *= (1.0 - gamma)
            B[amins, cols] += gamma
            D *= (1.0 - gamma)
            D += gamma * KEt1
        return B

    def _Mcsc(self):
        """CSC copy for column selection; built once and reused."""
        if getattr(self, "_M_csc", None) is None or \
                self._M_csc.shape != self.M_.shape:
            self._M_csc = self.M_.tocsc()
        return self._M_csc

    def step(self) -> float:
        self.A_ = self._update_A(self.B_, self.A_)
        self.B_ = self._update_B(self.A_, self.B_)
        rss = _rss(self.M_, self.A_, self.B_, self._m_fro_sq)
        self.RSS_iters.append(rss)
        return rss

    def fit(self, adata=None, max_iter: int = 50, min_iter: int = 10,
            initial_archetypes=None, log_every: int = 10):
        """Iterate to solve for A, B and write assignments to adata.obs['SEACell']."""
        if self.M_ is None:
            if adata is None:
                raise ValueError("Either adata or a precomputed kernel is required")
            self.build_kernel(adata)
        if max_iter < min_iter:
            warn(f"max_iter {max_iter} < min_iter {min_iter}; lowering min_iter to match")
            min_iter = max_iter
        self._adata = adata
        self.initialize(adata=adata, initial_archetypes=initial_archetypes)

        threshold = self.convergence_epsilon * self.RSS_iters[0]
        t0 = time.time()
        n_iter = 0
        while (not self.converged_ and n_iter < max_iter) or n_iter < min_iter:
            n_iter += 1
            self.step()
            if self.verbose and (n_iter == 1 or n_iter % log_every == 0):
                log(f"iter {n_iter}: RSS {self.RSS_iters[-1]:.4f}"
                    f" (delta {abs(self.RSS_iters[-2] - self.RSS_iters[-1]):.5f},"
                    f" threshold {threshold:.5f}) {time.time() - t0:.1f}s")
            if abs(self.RSS_iters[-2] - self.RSS_iters[-1]) < threshold:
                self.converged_ = True
        self.n_iter_ = n_iter
        if self.verbose:
            log(f"{'converged' if self.converged_ else 'not converged'}: {n_iter} iters,"
                f" RSS {self.RSS_iters[0]:.4f} -> {self.RSS_iters[-1]:.4f},"
                f" {time.time() - t0:.1f}s")
        if not self.converged_:
            warn(f"Did not converge in {max_iter} iterations. Increase max_iter or"
                 " reconsider n_metacells")
        if adata is not None:
            adata.obs["SEACell"] = self.hard_assignments()["SEACell"].values
        self.balance_ = size_balance(self.hard_assignments()["SEACell"].values)
        if self.verbose:
            b = self.balance_
            log(f"Size distribution: median {b['size_median']:.0f}"
                f" [{b['size_min']}, {b['size_max']}]"
                f" / Gini {b['gini']:.3f} / max share {b['max_share']:.1%}"
                f" / singletons {b['n_singleton']}")
        share, even = self.balance_["max_share"], 1.0 / self.k
        if share > MAX_SIZE_SHARE_ABS and share > MAX_SIZE_SHARE_RATIO * even:
            warn(f"Metacell size distribution is degenerate: the largest metacell"
                 f" holds {share:.1%} of cells (even split would be {even:.1%}, {share/even:.0f}x)."
                 f" {self.balance_['n_singleton']} singleton metacells."
                 " The objective has been optimized too far and this result is"
                 " not usable as metacells. Set line_search=False and use a"
                 " smaller fw_iters (15-50)")
        return self

    # --- Output --------------------------------------------------------
    def hard_assignments(self, index=None) -> pd.DataFrame:
        """Assign by argmax of A (same rule and naming as SEACells)."""
        if self.A_ is None:
            raise RuntimeError("Call fit() first")
        lab = [f"SEACell-{i}" for i in self.A_.argmax(0)]
        if index is None and self._adata is not None:
            index = self._adata.obs_names
        df = pd.DataFrame({"SEACell": lab},
                          index=index if index is not None else np.arange(len(lab)))
        df.index.name = "index"
        return df

    def soft_assignments(self, n_top: int = 5):
        A = self.A_.T.copy()
        names = np.array([f"SEACell-{i}" for i in range(self.k)])
        labels, weights = [], []
        for _ in range(n_top):
            j = A.argmax(1)
            labels.append(names[j])
            weights.append(A[np.arange(A.shape[0]), j])
            A[np.arange(A.shape[0]), j] = -np.inf
        return np.vstack(labels).T, np.vstack(weights).T

    def archetype_matrix(self) -> np.ndarray:
        """Z = Bᵀ K (k×n), evaluated without materializing K."""
        return (_KX(self.M_, self.B_)).T


def size_balance(labels) -> dict:
    """Measure skew in the metacell size distribution.

    Returns the Gini coefficient and the fraction of cells held by the
    largest metacell. The archetypal objective ‖M − MBA‖² favors solutions
    where one archetype exactly reproduces a single cell and pushes the rest
    into one giant metacell; RSS and cell-type purity both misread this
    degeneracy as improvement, so size distribution must be checked
    separately.
    """
    sz = pd.Series(labels).value_counts().values.astype(float)
    n = sz.size
    srt = np.sort(sz)
    gini = float(2.0 * np.sum(np.arange(1, n + 1) * srt) / (n * srt.sum())
                 - (n + 1) / n) if n > 1 else 0.0
    return {"n_metacells": int(n), "gini": gini,
            "max_share": float(sz.max() / sz.sum()),
            "size_min": int(sz.min()), "size_median": float(np.median(sz)),
            "size_max": int(sz.max()),
            "n_singleton": int((sz == 1).sum())}


# ---------------------------------------------------------------------------
# Evaluation metrics (no palantir dependency)
# ---------------------------------------------------------------------------

def diffusion_components(X: np.ndarray, k: int = 30, n_comp: int = 10,
                         seed: int = 0, lam_tol: float = 1e-6) -> np.ndarray:
    """Diffusion components, computed natively with the same construction as
    palantir's run_diffusion_maps + determine_multiscale_space (adaptive-
    bandwidth affinity -> row normalization -> eigendecomposition -> scale by
    lambda/(1-lambda)). Absolute scale won't match palantir; compare by rank.

    Eigenvalues at lambda = 1 are always dropped. Besides the stationary
    component, a disconnected kNN graph also produces lambda = 1 eigenvalues
    for each connected component's indicator vector (observed on synthetic
    data). Keeping these pushes lambda/(1-lambda) to ~1e20 and lets that one
    component dominate the variance-averaged compactness metric.
    """
    from sklearn.neighbors import NearestNeighbors
    from scipy.sparse.linalg import eigs

    n = X.shape[0]
    k = int(min(k, max(2, n - 1)))
    nn = NearestNeighbors(n_neighbors=k).fit(X)
    d, idx = nn.kneighbors(X)
    adaptive = max(1, k // 3)
    sigma = np.maximum(d[:, adaptive], 1e-12)
    rows = np.repeat(np.arange(n), k)
    cols = idx.ravel()
    vals = np.exp(-(d.ravel() ** 2) / (sigma[rows] ** 2))
    W = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    W = W + W.T
    T = sp.diags(1.0 / np.maximum(np.asarray(W.sum(1)).ravel(), 1e-12)) @ W
    m = int(min(n_comp + 12, n - 2))
    vals_e, vecs = eigs(T, k=m, which="LM", v0=np.ones(n) / np.sqrt(n),
                        maxiter=10_000)
    order = np.argsort(-np.real(vals_e))
    lam = np.real(vals_e)[order]
    vec = np.real(vecs)[:, order]
    keep = lam < 1.0 - lam_tol
    n_drop = int((~keep).sum())
    if n_drop > 1:
        warn(f"{n_drop} eigenvalues at lambda=1 (kNN graph is disconnected)."
             " Dropping all of them when building diffusion components")
    lam, vec = lam[keep][:n_comp], vec[:, keep][:, :n_comp]
    if lam.size == 0:
        raise ValueError("No eigenvalues below lambda=1 were found. Increase k")
    scale = lam / np.maximum(1.0 - lam, 1e-12)
    return vec * scale


def compactness(adata, use_rep: str = "X_pca", key: str = "SEACell",
                space: Literal["diffusion", "pca"] = "diffusion") -> pd.DataFrame:
    """Metacell compactness; lower means cells within it are in more similar states.

    Defined as in SEACells, via mean variance of palantir diffusion
    components, computed here with native diffusion components
    (`space="pca"` uses raw PCA coordinates instead).
    """
    X = np.asarray(adata.obsm[use_rep], dtype=float)
    C = X if space == "pca" else diffusion_components(X)
    df = pd.DataFrame(C, index=adata.obs_names)
    df[key] = adata.obs[key].astype(str).values
    out = df.groupby(key).var(numeric_only=True).mean(axis=1)
    return pd.DataFrame({"compactness": out})


def separation(adata, use_rep: str = "X_pca", key: str = "SEACell",
               nth_nbr: int = 1, cluster: str | None = None,
               space: Literal["diffusion", "pca"] = "diffusion") -> pd.DataFrame:
    """Distance to the nearest other metacell; larger means better separated."""
    from sklearn.neighbors import NearestNeighbors

    X = np.asarray(adata.obsm[use_rep], dtype=float)
    C = X if space == "pca" else diffusion_components(X)
    df = pd.DataFrame(C, index=adata.obs_names)
    df[key] = adata.obs[key].astype(str).values
    cent = df.groupby(key).mean(numeric_only=True)
    nn = NearestNeighbors(n_neighbors=min(nth_nbr, len(cent) - 1)).fit(cent.values)
    dists, nbrs = nn.kneighbors()
    res = pd.DataFrame({"separation": dists[:, nth_nbr - 1]}, index=cent.index)
    if cluster is not None and cluster in adata.obs:
        maj = (adata.obs.groupby(adata.obs[key].astype(str))[cluster]
               .agg(lambda x: x.value_counts().index[0]))
        same = maj.values[nbrs[:, nth_nbr - 1]] == maj.reindex(cent.index).values
        res = res[same]
    return res


def celltype_purity(adata, col: str, key: str = "SEACell") -> pd.DataFrame:
    """The most common label within each metacell, and its fraction."""
    g = adata.obs.groupby(adata.obs[key].astype(str))[col]
    top = g.agg(lambda x: x.value_counts().index[0])
    frac = g.agg(lambda x: x.value_counts().iloc[0] / len(x))
    return pd.DataFrame({col: top, f"{col}_purity": frac})
