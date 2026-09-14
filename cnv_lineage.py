#!/usr/bin/env python3
"""cnv_lineage.py — Build a tumor cell lineage from CNV calls (Dollo parsimony).

Chromosomal loss is treated as irreversible: a region lost in one cell
cannot reappear in its descendants without independently reacquiring the
lost sequence. Dollo parsimony formalizes this directly — a derived trait
may be gained at most once on the tree, but lost any number of times — so
CNV losses are modeled as Dollo traits. Gains (amplifications) don't share
this asymmetry (the same region can be amplified independently more than
once, and can also regress), so this module treats loss as the primary
trait and gain as a separate, secondary one (`event_kind` selects which).

Pipeline:
1. `cnv_events()`               log2 ratios (metacell x bin) -> a binary
                                 loss/gain event matrix, relative to normal
                                 metacells.
2. `four_gamete_violations()`   tests whether a perfect phylogeny (no
                                 back-mutation) is even possible before a
                                 tree is built: if patterns 11/10/01 all
                                 co-occur for a pair of traits, no single
                                 tree can explain them.
3. `dollo_tree()`                builds a tree from containment relationships
                                 between event carrier sets (carrier set of
                                 trait j subset of trait i's), dropping
                                 traits that violate the four-gamete test
                                 (lowest support first) and reporting how
                                 many were dropped.
4. `nj_tree()`                   a neighbor-joining tree for comparison,
                                 without the Dollo constraint.
5. `lineage_null_test()`         permutation test (row/column sums of the
                                 event matrix preserved) for whether the
                                 observed parsimony improvement could be due
                                 to chance. **No tree is presented as a
                                 lineage if there's no structure.**
6. `to_newick()`                 Newick string output.

Produces a Newick string plus a per-branch event table (which region was
lost on which branch). See README.md for a summary and validation results.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import pandas as pd

try:
    from scrna_common import log, warn
except Exception:  # pragma: no cover
    def log(m: str) -> None:
        print(f"[lineage] {m}", flush=True)

    def warn(m: str) -> None:
        print(f"[lineage][WARNING] {m}", flush=True)


__version__ = "1.0"

# Log2-ratio threshold for calling loss/gain, set from the normal-metacell
# distribution's quantiles (1st/99th percentile) by default rather than a
# fixed value, since spread varies a lot with depth and expression program.
DEFAULT_EVENT_ALPHA = 0.05
# Upper bound on the fraction of normal metacells an event may appear in.
# A tumor-derived loss should not appear in normal cells.
MAX_NORMAL_RATE = 0.05
# Minimum carrier count to accept an event (fewer carriers = noise)
MIN_EVENT_SUPPORT = 3
# Minimum count for the four-gamete test to call a pattern "present"
MIN_GAMETE_COUNT = 2
# Relative tolerance for the four-gamete test: what fraction of the smaller
# carrier set is allowed as slack
GAMETE_TOL = 0.10
# Minimum effect size (how much the laminarity deficit drops vs. the null)
# to call structure present. Measured: 0.66 for a planted lineage, -0.001
# for random data, 0.008 on real Case1 data.
MIN_STRUCTURE_EFFECT = 0.20


@dataclass
class Node:
    name: str
    children: list["Node"] = field(default_factory=list)
    members: list[str] = field(default_factory=list)   # metacells attached directly to this node
    events: list[str] = field(default_factory=list)    # events that occurred on this branch
    support: int = 0


# ---------------------------------------------------------------------------
# 1. Event matrix
# ---------------------------------------------------------------------------

def cnv_events(profiles: pd.DataFrame, normal_mask: np.ndarray, *,
               kind: Literal["loss", "gain", "both"] = "loss",
               alpha: float = DEFAULT_EVENT_ALPHA,
               min_support: int = MIN_EVENT_SUPPORT,
               max_normal_rate: float = MAX_NORMAL_RATE
               ) -> tuple[pd.DataFrame, dict]:
    """Convert a log2-ratio table (metacell x region) into a binary event matrix.

    **Thresholds are set by order statistics.** For each region, normal
    metacell values are sorted and the k-th smallest (k = floor(max_normal_rate
    x n_normal)) is used as the loss threshold (k-th largest for gain). This
    means:

      * The threshold tolerates up to k contaminating tumor metacells in the
        normal reference without being pulled by them. Tumor contamination
        of the normal reference is the most common CNV-analysis failure mode
        — infercnv hides it entirely by rounding reference values to 0 — and
        a plain quantile without this trimming lets a few contaminants pull
        the threshold down enough that **no real loss gets called at all**
        (measured: with 2 contaminants out of 60, an untrimmed quantile
        called zero chr1 events).
      * By construction, the probability of a normal metacell crossing the
        threshold is exactly max_normal_rate, so the support-count floor
        below can be derived from the binomial distribution.

    The support floor is set from the noise level: under the null, a
    malignant metacell's carrier count follows Binom(n_malignant,
    max_normal_rate), and the observed count must exceed the
    Bonferroni-corrected quantile over regions. A fixed `min_support` alone
    lets noise through (measured: with 40 regions and 115 malignant
    metacells, a fixed threshold of 3 left 14 spurious events across 33
    unrelated regions).
    """
    from scipy import stats as _st

    if profiles.isna().any().any():
        profiles = profiles.fillna(0.0)
    nm = np.asarray(normal_mask, bool)
    n_norm, n_mal = int(nm.sum()), int((~nm).sum())
    if n_norm < 10:
        raise ValueError(f"Only {n_norm} normal metacells available")
    if n_mal < 3:
        raise ValueError(f"Only {n_mal} malignant metacells available")
    V = profiles.values
    srt = np.sort(V[nm], axis=0)
    n_reg = profiles.shape[1]

    k = int(np.floor(max_normal_rate * n_norm))
    k = int(min(max(k, 0), n_norm - 2))
    lo = srt[k]
    hi = srt[n_norm - 1 - k]
    eff_rate = (k + 1) / (n_norm + 1)          # actual contamination rate tolerated

    noise_thr = int(_st.binom.isf(alpha / max(n_reg, 1), n_mal, eff_rate)) + 1
    support_thr = max(min_support, noise_thr)

    cols, names, rej = [], [], {"low_support": 0, "seen_in_normal": 0}

    def _add(v: np.ndarray, label: str) -> None:
        if int((~nm & v).sum()) < support_thr:
            rej["low_support"] += 1
            return
        if int((nm & v).sum()) > k:
            rej["seen_in_normal"] += 1
            return
        cols.append(v); names.append(label)

    if kind in ("loss", "both"):
        for j2, c in enumerate(profiles.columns):
            _add(V[:, j2] < lo[j2], f"loss:{c}")
    if kind in ("gain", "both"):
        for j2, c in enumerate(profiles.columns):
            _add(V[:, j2] > hi[j2], f"gain:{c}")
    if not cols:
        raise ValueError(
            "No events satisfied the criteria"
            f" (low support {rej['low_support']}, also seen in normal {rej['seen_in_normal']};"
            f" support floor {support_thr}). Loosen alpha or raise max_normal_rate")
    B = pd.DataFrame(np.column_stack(cols), index=profiles.index, columns=names)
    info = {"n_events": int(B.shape[1]), "n_metacells": int(B.shape[0]),
            "alpha": alpha, "min_support": min_support,
            "support_threshold": int(support_thr), "noise_threshold": int(noise_thr),
            "max_normal_carriers": int(k), "effective_normal_rate": float(eff_rate),
            "kind": kind, "rejected": rej,
            "support": {c: int(B[c].sum()) for c in B.columns}}
    log(f"Event matrix: {B.shape[0]} metacells x {B.shape[1]} events"
        f" (loss {sum(c.startswith('loss') for c in B.columns)} /"
        f" gain {sum(c.startswith('gain') for c in B.columns)})"
        f" / normal contamination tolerated {k} ({eff_rate:.1%})"
        f" / support floor {support_thr} (noise-derived floor {noise_thr})"
        f" / dropped: low support {rej['low_support']}, also seen in normal {rej['seen_in_normal']}")
    return B, info


# ---------------------------------------------------------------------------
# 2. Whether a perfect phylogeny is possible (four-gamete test)
# ---------------------------------------------------------------------------

def laminarity_deficit(B: pd.DataFrame) -> int:
    """Deficit from laminarity: directly measures how well Dollo's required structure holds.

    Under Dollo (gain at most once), event i's carrier set S_i corresponds
    to a single subtree, so for any two events S_i and S_j must be either
    **nested or disjoint** (a laminar family). The degree to which they are
    not is measured as

        D = Σ_{i<j} min(|S_i ∩ S_j|, |S_i \\ S_j|, |S_j \\ S_i|)

    A pair contributes 0 if any of the three terms is 0 (nested or
    disjoint); D = 0 for a perfectly nested family.

    **Computed without building a tree.** Using a tree's own parsimony score
    as the test statistic ties the result to whether tree-building happens
    to drop events, which degenerates a support-preserving permutation test
    (observed and null collapse to the same value — this was observed in
    practice). D depends only on carrier sets, so it compares correctly
    against a support-preserving null.
    """
    X = B.values.astype(bool)
    A = X.T.astype(np.int64)
    n11 = A @ A.T
    col = X.sum(0).astype(np.int64)
    n10 = col[:, None] - n11
    n01 = col[None, :] - n11
    M = np.minimum(np.minimum(n11, n10), n01)
    np.fill_diagonal(M, 0)
    return int(M.sum() // 2)


def four_gamete_violations(B: pd.DataFrame,
                           min_count: int = MIN_GAMETE_COUNT,
                           tol: float = GAMETE_TOL) -> dict:
    """For each pair of traits, count whether patterns 11/10/01 all co-occur.

    If all three occur, neither trait is an ancestor of the other nor are
    they mutually exclusive, so no single tree explains them without
    allowing back-mutation or independent gains. Satisfying this condition
    for every pair is necessary and sufficient for a perfect phylogeny to
    exist (Gusfield 1991). **Check this before drawing a tree.**
    """
    X = B.values.astype(bool)
    m = X.shape[1]
    n11 = X.T.astype(np.int32) @ X.astype(np.int32)
    col = X.sum(0).astype(np.int32)
    n10 = col[:, None] - n11
    n01 = col[None, :] - n11
    # Tolerance is the larger of the absolute and relative bounds. Calling
    # 2-3 noise-driven overlaps between 50-carrier events a "violation"
    # would break even a cleanly nested lineage (measured: only 1/8 events
    # accepted on a planted lineage without this). Uses tol times the
    # smaller carrier set as the floor.
    smaller = np.minimum(col[:, None], col[None, :])
    thr = np.maximum(min_count, np.ceil(tol * smaller).astype(np.int64))
    ok = (n11 >= thr) & (n10 >= thr) & (n01 >= thr)
    np.fill_diagonal(ok, False)
    n_pairs = m * (m - 1) // 2
    n_bad = int(ok.sum() // 2)
    # Number of violating partners per trait
    per = pd.Series(ok.sum(1), index=B.columns)
    out = {"n_events": m, "n_pairs": n_pairs, "n_violating_pairs": n_bad,
           "violation_rate": n_bad / max(n_pairs, 1),
           "per_event_violations": per.to_dict()}
    log(f"Four-gamete test: {n_bad:,}/{n_pairs:,} violating pairs"
        f" ({out['violation_rate']:.1%})")
    if out["violation_rate"] > 0.05:
        warn(f"{out['violation_rate']:.1%} of pairs violate the four-gamete"
             " condition. A large fraction of the data can't be explained by"
             " a single tree, so branch interpretation should be cautious."
             " Likely causes: (a) the same region was lost independently"
             " more than once, (b) a metacell mixes cells from multiple"
             " lineages, (c) event calling is picking up noise")
    return out


# ---------------------------------------------------------------------------
# 3. Tree via Dollo parsimony
# ---------------------------------------------------------------------------

def dollo_tree(B: pd.DataFrame, *, min_count: int = MIN_GAMETE_COUNT,
               tol: float = GAMETE_TOL, root_name: str = "diploid_root",
               verbose: bool = True) -> tuple[Node, dict]:
    """Build a tree from containment relationships.

    Under Dollo (gain at most once), the set of metacells carrying trait i
    corresponds to a single subtree, so any two carrier sets must be either
    nested or disjoint (a laminar family). Pairs that aren't are four-gamete
    violations and are dropped, lowest support first.

    Returns the number of traits dropped and their total support; if this
    is large, forcing a tree is itself distorting the conclusion.
    """
    X = B.values.astype(bool)
    order = np.argsort(-X.sum(0))          # most-supported first = closest to root
    names = list(B.columns)
    kept: list[int] = []
    dropped: list[str] = []
    sets: list[set[int]] = []
    for j in order:
        sj = set(np.where(X[:, j])[0])
        conflict = False
        for si in sets:
            t = max(min_count, int(np.ceil(tol * min(len(sj), len(si)))))
            if (len(sj & si) >= t and len(sj - si) >= t and len(si - sj) >= t):
                conflict = True
                break
        if conflict:
            dropped.append(names[j])
        else:
            kept.append(j); sets.append(sj)

    # Build a containment forest: each trait's parent = the smallest trait
    # that properly contains it
    idx = {j: set(np.where(X[:, j])[0]) for j in kept}
    parent: dict[int, int | None] = {}
    for j in kept:
        best, best_size = None, None
        for i in kept:
            if i == j:
                continue
            slack = max(min_count, int(np.ceil(tol * len(idx[j]))))
            if (len(idx[j] - idx[i]) < slack and len(idx[i]) > len(idx[j])):
                if best_size is None or len(idx[i]) < best_size:
                    best, best_size = i, len(idx[i])
        parent[j] = best

    nodes = {j: Node(name=names[j], events=[names[j]], support=len(idx[j]))
             for j in kept}
    root = Node(name=root_name)
    for j in kept:
        (nodes[parent[j]] if parent[j] is not None else root).children.append(nodes[j])

    # Attach each metacell to the most derived trait it carries
    for r, mc in enumerate(B.index):
        carried = [j for j in kept if X[r, j]]
        if not carried:
            root.members.append(str(mc)); continue
        deepest = min(carried, key=lambda j: len(idx[j]))
        nodes[deepest].members.append(str(mc))

    stats = {
        "n_events_kept": len(kept), "n_events_dropped": len(dropped),
        "dropped": dropped,
        "dropped_support": int(sum(int(X[:, names.index(d)].sum()) for d in dropped)),
        "n_internal_nodes": len(kept),
        "n_metacells_at_root": len(root.members),
        "parsimony_score": _parsimony_score(B, root),
    }
    if verbose:
        log(f"Dollo tree: kept {len(kept)} events / dropped {len(dropped)}"
            f" / internal nodes {len(kept)} / metacells at root {len(root.members)}")
    if dropped and verbose:
        warn(f"Dropped {len(dropped)} events that violate the four-gamete"
             f" condition (total support {stats['dropped_support']}). The"
             " tree is simplified by that much")
    return root, stats


def _subtree_members(n: Node) -> list[str]:
    out = list(n.members)
    for c in n.children:
        out += _subtree_members(c)
    return out


def _parsimony_score(B: pd.DataFrame, root: Node) -> int:
    """Number of changes (Dollo cost) needed for this tree to explain **all events**.

    For each event j's carrier set S_j, finds the smallest subtree fully
    containing S_j:
        cost = 1 (gain) + (that subtree's metacell count - |S_j|) (losses)
    **Events dropped while building the tree are still counted here**, to
    close the loophole where dropping more events would otherwise lower the
    score — without this, "trees that drop more events look more
    parsimonious" (this was observed to actually happen).
    """
    subs = {id(n): set(_subtree_members(n)) for n in _iter_nodes(root)}
    order = sorted(_iter_nodes(root), key=lambda n: len(subs[id(n)]))
    X = B.values.astype(bool)
    total = 0
    for j in range(X.shape[1]):
        S = {str(m) for m, v in zip(B.index, X[:, j]) if v}
        if not S:
            continue
        host = None
        for n in order:                      # smallest subtree first, take the first that contains it
            if S <= subs[id(n)]:
                host = n
                break
        if host is None:
            total += 1 + (len(subs[id(root)]) - len(S))
        else:
            total += 1 + (len(subs[id(host)]) - len(S))
    return int(total)


def _iter_nodes(n: Node):
    yield n
    for c in n.children:
        yield from _iter_nodes(c)


# ---------------------------------------------------------------------------
# 4. Comparison baseline: neighbor joining
# ---------------------------------------------------------------------------

def nj_tree(B: pd.DataFrame) -> tuple[Node, dict]:
    """Comparison baseline that ignores the Dollo constraint: Hamming distance + neighbor joining.

    Used to measure the benefit of imposing irreversibility. Returns only tree topology.
    """
    X = B.values.astype(float)
    n = X.shape[0]
    D = np.abs(X[:, None, :] - X[None, :, :]).sum(2)
    labels = [Node(name=str(m), members=[str(m)]) for m in B.index]
    active = list(range(n))
    Dm = D.copy().astype(float)
    nxt = 0
    while len(active) > 2:
        k = len(active)
        sub = Dm[np.ix_(active, active)]
        u = sub.sum(1) / (k - 2)
        Q = sub - u[:, None] - u[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = np.unravel_index(np.argmin(Q), Q.shape)
        ai, aj = active[i], active[j]
        nxt += 1
        new = Node(name=f"nj{nxt}", children=[labels[ai], labels[aj]])
        labels.append(new)
        newd = 0.5 * (Dm[ai, :] + Dm[aj, :] - Dm[ai, aj])
        Dm = np.pad(Dm, ((0, 1), (0, 1)), constant_values=0.0)
        Dm[-1, :-1] = newd; Dm[:-1, -1] = newd; Dm[-1, -1] = 0.0
        active = [a for a in active if a not in (ai, aj)] + [len(labels) - 1]
    root = Node(name="nj_root", children=[labels[a] for a in active])
    return root, {"method": "neighbor-joining", "n_leaves": n}


# ---------------------------------------------------------------------------
# 5. Testing whether structure is present
# ---------------------------------------------------------------------------

def lineage_null_test(B: pd.DataFrame, *, n_perm: int = 200,
                      method: Literal["independent", "curveball", "both"] = "both",
                      normal_mask: np.ndarray | None = None,
                      seed: int = 0) -> dict:
    """Permutation test for whether event co-occurrence shows lineage structure.

    Two null models are available; **both have tradeoffs, so both run by default.**

    ``independent``
        Shuffles each event's carriers independently (support preserved,
        co-occurrence destroyed). Directly asks "is the nesting just chance?"
        Doesn't preserve per-metacell event counts, so it's somewhat
        liberal if events cluster on a few metacells.
    ``curveball``
        Preserves both row sums (events per metacell) and column sums
        (support) simultaneously. More conservative, but **degenerates to
        zero null variance if there are too few distinct row patterns to
        swap** (measured: SD 0.0, p=1.0 on synthetic data with 3 planted
        clones). Falls back to `independent` with a warning when degenerate.

    If there's no significant structure, the tree should not be read as a lineage.
    """
    rng = np.random.default_rng(seed)
    obs = laminarity_deficit(B)
    X = B.values.astype(bool)
    nm = (np.zeros(X.shape[0], bool) if normal_mask is None
          else np.asarray(normal_mask, bool))
    pool = np.where(~nm)[0] if nm.any() else np.arange(X.shape[0])

    def _score(Y):
        return laminarity_deficit(pd.DataFrame(Y, index=B.index, columns=B.columns))

    out: dict = {"observed": int(obs), "n_perm": int(n_perm)}
    runs = ("independent", "curveball") if method == "both" else (method,)
    for meth in runs:
        vals = []
        for _ in range(n_perm):
            if meth == "curveball":
                Y = _curveball(X, rng, n_steps=5 * X.size)
            else:
                Y = np.zeros_like(X)
                for j2 in range(X.shape[1]):
                    k = int(X[pool, j2].sum())
                    if k:
                        Y[rng.choice(pool, k, replace=False), j2] = True
                    # Keep the normal-side carriers unchanged (a handful, by
                    # construction of the threshold)
                    Y[np.where(nm & X[:, j2])[0], j2] = True
            v = _score(Y)
            if v is not None:
                vals.append(v)
        a = np.asarray(vals, float)
        if a.size < 10:
            continue
        pv = float((np.sum(a <= obs) + 1) / (a.size + 1))
        out[meth] = {"null_mean": float(a.mean()), "null_sd": float(a.std()),
                     "p_value": pv, "degenerate": bool(a.std() < 1e-9)}
        log(f"  null[{meth}]: laminarity deficit observed {obs:,} / null {a.mean():,.0f}"
            f" +/- {a.std():.0f} -> p={pv:.4f}"
            + ("  <- zero variance, unusable" if a.std() < 1e-9 else ""))

    use = "independent"
    if "curveball" in out and not out["curveball"]["degenerate"] and \
            "independent" in out:
        use = "curveball"          # prefer the more conservative one when usable
    elif "curveball" in out and out["curveball"]["degenerate"]:
        warn("The curveball null degenerated (too few distinct row patterns)."
             " Falling back to the independent null")
    if use not in out:
        use = next(iter(k for k in ("independent", "curveball") if k in out))
    chosen = out[use]
    out["null_used"] = use
    out["p_value"] = chosen["p_value"]
    out["null_mean"] = chosen["null_mean"]
    out["null_sd"] = chosen["null_sd"]
    effect = 1.0 - obs / max(chosen["null_mean"], 1e-9)
    out["effect_size"] = float(effect)
    out["min_effect"] = MIN_STRUCTURE_EFFECT
    # Don't decide on the p-value alone: with hundreds of metacells even a
    # <1% deviation can be p<0.05 without being a lineage-readable structure
    # (measured on real Case1 data: effect size 0.008 at p=0.005).
    out["has_structure"] = bool(chosen["p_value"] < 0.05
                                and effect >= MIN_STRUCTURE_EFFECT)
    log(f"Lineage structure test[{use}]: p={out['p_value']:.4f} /"
        f" effect size {effect:+.3f} (floor {MIN_STRUCTURE_EFFECT})"
        f" -> structure {'present' if out['has_structure'] else 'not supported'}")
    if not out["has_structure"]:
        warn("Event co-occurrence is consistent with chance. Any tree drawn"
             " from this should not be read as a lineage")
    return out


def _curveball(X: np.ndarray, rng, n_steps: int) -> np.ndarray:
    """Shuffle a binary matrix preserving row and column sums (curveball algorithm)."""
    Y = X.copy()
    n, m = Y.shape
    rows = [set(np.where(Y[i])[0]) for i in range(n)]
    for _ in range(n_steps):
        i, j = rng.integers(n), rng.integers(n)
        if i == j:
            continue
        a, b = rows[i] - rows[j], rows[j] - rows[i]
        if not a or not b:
            continue
        k = min(len(a), len(b))
        sw = rng.permutation(list(a))[:k]
        sb = rng.permutation(list(b))[:k]
        for x, y in zip(sw, sb):
            rows[i].discard(x); rows[i].add(y)
            rows[j].discard(y); rows[j].add(x)
    Z = np.zeros_like(Y)
    for i in range(n):
        Z[i, list(rows[i])] = True
    return Z


# ---------------------------------------------------------------------------
# 6. Newick
# ---------------------------------------------------------------------------

def to_newick(root: Node, *, with_members: bool = True) -> str:
    """Newick format. Branch length is the number of events on that branch."""
    def esc(s: str) -> str:
        return str(s).replace(" ", "_").replace(",", "_").replace(":", "_") \
                     .replace("(", "[").replace(")", "]").replace(";", "_")

    def rec(n: Node) -> str:
        kids = [rec(c) for c in n.children]
        if with_members:
            kids += [f"{esc(m)}:0" for m in n.members]
        blen = max(len(n.events), 0)
        label = esc(n.name)
        if kids:
            return f"({','.join(kids)}){label}:{blen}"
        return f"{label}:{blen}"

    return rec(root) + ";"


def branch_table(root: Node) -> pd.DataFrame:
    """Table of events per branch and the metacells attached to each."""
    rows = []

    def rec(n: Node, depth: int, parent: str | None):
        sub = set(n.members)
        for c in n.children:
            sub |= set(_all_members(c))
        rows.append({"node": n.name, "parent": parent, "depth": depth,
                     "events": ";".join(n.events),
                     "n_events": len(n.events),
                     "n_metacells_here": len(n.members),
                     "n_metacells_subtree": len(sub),
                     "metacells": ";".join(n.members)})
        for c in n.children:
            rec(c, depth + 1, n.name)

    rec(root, 0, None)
    return pd.DataFrame(rows)


_all_members = _subtree_members


# ---------------------------------------------------------------------------
# Run the full pipeline
# ---------------------------------------------------------------------------

def build_lineage(profiles: pd.DataFrame, normal_mask: np.ndarray, *,
                  kind: str = "loss", alpha: float = DEFAULT_EVENT_ALPHA,
                  min_support: int = MIN_EVENT_SUPPORT,
                  n_perm: int = 200, seed: int = 0,
                  out_dir: str | Path | None = None) -> dict:
    """Run event calling -> condition checks -> significance test -> tree -> Newick, end to end.

    **If the significance test doesn't support structure, a tree is still
    built but `has_structure=False` is set.** Callers should check this flag
    before deciding whether to display the figure.
    """
    B, einfo = cnv_events(profiles, normal_mask, kind=kind, alpha=alpha,
                          min_support=min_support)
    fg = four_gamete_violations(B)
    null = lineage_null_test(B, n_perm=n_perm, seed=seed,
                             normal_mask=np.asarray(normal_mask, bool))
    root, stats = dollo_tree(B)
    nwk = to_newick(root)
    tab = branch_table(root)
    res = {"events": einfo, "four_gamete": fg, "null_test": null,
           "tree": stats, "newick": nwk,
           "has_structure": null["has_structure"]}
    if out_dir:
        d = Path(out_dir); d.mkdir(parents=True, exist_ok=True)
        (d / "cnv_lineage.nwk").write_text(nwk + "\n", encoding="utf-8")
        tab.to_csv(d / "cnv_lineage_branches.csv", index=False)
        B.astype(int).to_csv(d / "cnv_lineage_events.csv")
        (d / "cnv_lineage_report.json").write_text(
            json.dumps({k: v for k, v in res.items() if k != "newick"},
                       indent=2, ensure_ascii=False, default=float),
            encoding="utf-8")
        log(f"Wrote: {d}/cnv_lineage.nwk and 3 other files")
    return res


# ---------------------------------------------------------------------------
# Entry point from AnnData
# ---------------------------------------------------------------------------

def profiles_from_adata(cnv_adata, key: str = "cnv",
                        aggregate: Literal["window", "chromosome"] = "chromosome"
                        ) -> pd.DataFrame:
    """Turn obsm['X_cnv'] into a metacell x region table.

    With `aggregate="chromosome"`, windows are averaged per chromosome. This
    is the default because Case1's continuity test didn't support
    sub-chromosomal structure (p=0.09); using raw windows would give
    thousands of regions, letting noisy events dominate the tree shape.
    """
    import scipy.sparse as _sp

    X = cnv_adata.obsm[f"X_{key}"]
    X = np.asarray(X.todense()) if _sp.issparse(X) else np.asarray(X)
    chr_pos = (cnv_adata.uns.get(key, {}) or {}).get("chr_pos", {})
    if aggregate == "window" or not chr_pos:
        cols = [f"w{i}" for i in range(X.shape[1])]
        return pd.DataFrame(X, index=cnv_adata.obs_names.astype(str), columns=cols)
    names = list(chr_pos.keys())
    starts = list(chr_pos.values()) + [X.shape[1]]
    data = {c: X[:, s:e].mean(1) for c, s, e in zip(names, starts[:-1], starts[1:])
            if e > s}
    return pd.DataFrame(data, index=cnv_adata.obs_names.astype(str))
