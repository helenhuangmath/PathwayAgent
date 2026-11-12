"""
Vectorised enrichment statistics on the sparse library matrix (design §7.2).

All functions take the library matrix restricted to the selected sets (rows) and
work on gene *column indices*, so ORA over ~20k sets is one sparse mat-vec.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse, stats


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def _mask(n_genes: int, idx: np.ndarray) -> np.ndarray:
    m = np.zeros(n_genes, dtype=np.float32)
    m[idx] = 1.0
    return m


def _hits(M: sparse.csr_matrix, rows: np.ndarray, q_mask: np.ndarray, genes: np.ndarray) -> list[list[str]]:
    sub = sparse.csr_matrix(M[rows].multiply(q_mask[None, :]))
    sub.eliminate_zeros()   # multiply() keeps explicit zeros -> would leak non-query genes
    return [genes[sub[i].indices].tolist() for i in range(sub.shape[0])]


def ora(M: sparse.csr_matrix, rows: np.ndarray, universe: np.ndarray, query: np.ndarray,
        genes: np.ndarray, min_size: int, max_size: int) -> pd.DataFrame:
    """
    One-sided hypergeometric test (== Fisher exact, greater) for every set in `rows`.
    Background = universe genes annotated in at least one of these sets (the
    clusterProfiler / g:Profiler 'annotated' convention), so results match
    clusterProfiler::enricher on the same sets. Only sets with >=1 hit are returned.
    """
    n_genes = M.shape[1]
    sub = M[rows]
    annotated = np.asarray(sub.sum(axis=0)).ravel() > 0
    u_mask = _mask(n_genes, universe) * annotated
    universe = np.where(u_mask > 0)[0]
    q = np.intersect1d(query, universe)
    q_mask = _mask(n_genes, q)
    K = np.asarray(sub @ u_mask).ravel()
    k = np.asarray(sub @ q_mask).ravel()
    N, n = len(universe), len(q)
    keep = (K >= min_size) & (K <= max_size) & (k >= 1)
    rows_k, K, k = rows[keep], K[keep], k[keep]
    if len(rows_k) == 0 or n == 0:
        return pd.DataFrame()
    p = stats.hypergeom.sf(k - 1, N, K, n)
    a, b, c, d = k, n - k, K - k, N - K - n + k
    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return pd.DataFrame({
        "row": rows_k, "overlap": k.astype(int), "term_size": K.astype(int),
        "query_size": n, "universe_size": N, "pval": p, "effect": odds,
        "hit_genes": _hits(M, rows_k, q_mask, genes),
    })


def mhg(M: sparse.csr_matrix, rows: np.ndarray, universe: np.ndarray, ranked: np.ndarray,
        genes: np.ndarray, min_size: int, max_size: int, max_rank: int = 5000) -> pd.DataFrame:
    """
    Minimum-hypergeometric test over prefixes of a ranked list (as in g:Profiler
    'ordered' / GOrilla). The minimum over the k hit positions is Bonferroni-corrected
    by k (conservative approximation of the exact mHG null).
    `ranked` = gene columns in rank order (best first); genes after it are unranked.
    """
    n_genes = M.shape[1]
    sub = M[rows].tocsr()
    annotated = np.asarray(sub.sum(axis=0)).ravel() > 0
    u_mask = _mask(n_genes, universe) * annotated
    universe = np.where(u_mask > 0)[0]
    in_u = set(universe.tolist())
    ranked = np.array([g for g in ranked if g in in_u], dtype=int)[:max_rank]
    if len(ranked) == 0:
        return pd.DataFrame()
    pos = np.full(n_genes, -1, dtype=np.int64)
    pos[ranked] = np.arange(1, len(ranked) + 1)
    K = np.asarray(sub @ u_mask).ravel()
    N = len(universe)
    set_ids, js, ps = [], [], []
    for i in np.where((K >= min_size) & (K <= max_size))[0]:
        p_i = pos[sub[i].indices]
        p_i = np.sort(p_i[p_i > 0])
        if len(p_i):
            set_ids.append(np.full(len(p_i), i))
            js.append(np.arange(1, len(p_i) + 1))
            ps.append(p_i)
    if not set_ids:
        return pd.DataFrame()
    si, j, pp = np.concatenate(set_ids), np.concatenate(js), np.concatenate(ps)
    pv = stats.hypergeom.sf(j - 1, N, K[si], pp)
    df = pd.DataFrame({"i": si, "j": j, "pos": pp, "p": pv})
    best = df.loc[df.groupby("i")["p"].idxmin()]
    n_hits = df.groupby("i").size()
    best = best.set_index("i")
    p_corr = np.minimum(1.0, best["p"].values * n_hits.loc[best.index].values)
    rank_genes = genes[ranked]
    out = pd.DataFrame({
        "row": rows[best.index.values], "overlap": best["j"].values.astype(int),
        "term_size": K[best.index.values].astype(int), "query_size": best["pos"].values.astype(int),
        "universe_size": N, "pval": p_corr,
        "effect": (best["j"].values / best["pos"].values) / (K[best.index.values] / N),
    })
    hits = []
    for i, cut in zip(best.index.values, best["pos"].values):
        members = set(sub[i].indices.tolist())
        hits.append([rank_genes[t] for t in range(int(cut)) if ranked[t] in members])
    out["hit_genes"] = hits
    return out


def rank_test(M: sparse.csr_matrix, rows: np.ndarray, scored: np.ndarray, scores: np.ndarray,
              genes: np.ndarray, min_size: int, max_size: int) -> pd.DataFrame:
    """
    Two-sided Wilcoxon rank-sum (Mann-Whitney, normal approximation) of set members vs
    non-members among the scored genes. Effect = AUC (0.5 = no shift); direction = sign.
    """
    m = len(scored)
    ranks = stats.rankdata(scores)
    r_vec = np.zeros(M.shape[1], dtype=np.float64)
    r_vec[scored] = ranks
    s_mask = _mask(M.shape[1], scored).astype(np.float64)
    sub = M[rows].astype(np.float64)
    n1 = np.asarray(sub @ s_mask).ravel()
    R1 = np.asarray(sub @ r_vec).ravel()
    keep = (n1 >= min_size) & (n1 <= max_size) & (n1 < m)
    n1, R1, rows_k = n1[keep], R1[keep], rows[keep]
    if len(rows_k) == 0:
        return pd.DataFrame()
    n2 = m - n1
    U = R1 - n1 * (n1 + 1) / 2
    mu, sd = n1 * n2 / 2, np.sqrt(n1 * n2 * (m + 1) / 12)
    z = (U - mu) / sd
    p = 2 * stats.norm.sf(np.abs(z))
    auc = U / (n1 * n2)
    # members in the top/bottom of the ranking (leading genes) for the report
    order = np.argsort(-scores)
    top = set(scored[order[: max(1, m // 10)]].tolist())
    bottom = set(scored[order[-max(1, m // 10):]].tolist())
    hits = []
    for r, zz in zip(rows_k, z):
        idx = M[r].indices
        pick = top if zz > 0 else bottom
        hits.append(genes[[g for g in idx if g in pick]].tolist())
    return pd.DataFrame({"row": rows_k, "overlap": n1.astype(int), "term_size": n1.astype(int),
                         "query_size": m, "universe_size": m, "pval": p, "effect": auc,
                         "direction": np.sign(z).astype(int), "hit_genes": hits})


def gsea_blitz(signature: pd.Series, library: dict[str, list[str]], min_size: int, max_size: int,
               permutations: int = 2000, seed: int = 1) -> pd.DataFrame:
    """Pre-ranked GSEA via blitzGSEA (analytical null; seconds per library)."""
    import blitzgsea as blitz

    sig = pd.DataFrame({0: signature.index.values, 1: signature.values.astype(float)})
    res = blitz.gsea(sig, library, permutations=permutations, min_size=min_size, max_size=max_size,
                     processes=1, seed=seed, verbose=False, progress=False)
    if res is None or len(res) == 0:
        return pd.DataFrame()
    res = res.reset_index().rename(columns={"Term": "set_id", "index": "set_id"})
    return res
