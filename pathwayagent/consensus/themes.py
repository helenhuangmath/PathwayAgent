"""
Themes: cluster significant terms across databases by shared hit genes (design §9.3).

Metascape/DAVID-style: Cohen's kappa between terms' hit-gene vectors over the query
genes, average-linkage hierarchical clustering cut at kappa 0.3. Theme score rewards
cross-database agreement: ThemeScore = max member score * (1 + 0.1 * (n_databases - 1)).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

KAPPA_CUT = 0.3


def kappa_matrix(X: np.ndarray) -> np.ndarray:
    """Pairwise Cohen's kappa between rows of a binary matrix (terms x genes)."""
    X = X.astype(float)
    n = X.shape[1]
    a = X @ X.T                                  # both 1
    s = X.sum(axis=1)
    b = s[:, None] - a                           # row i only
    c = s[None, :] - a                           # row j only
    d = n - a - b - c
    po = (a + d) / n
    pe = ((a + b) * (a + c) + (c + d) * (b + d)) / n ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        k = (po - pe) / (1 - pe)
    k[~np.isfinite(k)] = 0
    np.fill_diagonal(k, 1.0)
    return k


def build_themes(cons: pd.DataFrame, max_terms: int = 150, min_tier3: int = 40) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (terms with theme_id, theme table)."""
    if cons.empty:
        return cons.assign(theme_id=np.nan), pd.DataFrame()
    cand = cons[cons["tier"].isin([1, 2])]
    t3 = cons[cons["tier"] == 3].head(max(min_tier3, max_terms - len(cand)))
    cand = pd.concat([cand, t3]).head(max_terms)
    cand = cand[cand["n_hit_genes"] > 0]
    cons = cons.copy()
    cons["theme_id"] = np.nan
    if cand.empty:
        return cons, pd.DataFrame()
    genes = sorted({g for h in cand["hit_genes"] for g in h.split(";") if g})
    gpos = {g: i for i, g in enumerate(genes)}
    X = np.zeros((len(cand), len(genes)), dtype=np.int8)
    for i, h in enumerate(cand["hit_genes"]):
        for g in h.split(";"):
            if g:
                X[i, gpos[g]] = 1
    if len(cand) == 1:
        labels = np.array([1])
    else:
        K = kappa_matrix(X)
        D = np.clip(1 - K, 0, 2)
        np.fill_diagonal(D, 0)
        Z = linkage(squareform(D, checks=False), method="average")
        labels = fcluster(Z, t=1 - KAPPA_CUT, criterion="distance")
    cand = cand.assign(_cl=labels)
    themes = []
    order = (cand.groupby("_cl")["consensus_score"].max().sort_values(ascending=False).index)
    remap = {cl: i + 1 for i, cl in enumerate(order)}
    for cl in order:
        m = cand[cand["_cl"] == cl].sort_values("consensus_score", ascending=False)
        lead = m.iloc[0]
        dbs = sorted(set(m["database"]))
        gene_counts = pd.Series([g for h in m["hit_genes"] for g in h.split(";") if g]).value_counts()
        themes.append({
            "theme_id": remap[cl], "theme_name": lead["term_name"], "lead_term": lead["term_key"],
            "theme_score": round(float(lead["consensus_score"]) * (1 + 0.1 * (len(dbs) - 1)), 2),
            "best_tier": int(m["tier"].replace(0, 9).min()), "n_terms": len(m),
            "databases": "|".join(dbs), "n_databases": len(dbs),
            "best_padj": float(m["best_padj"].min()),
            "member_terms": " | ".join(m["term_name"].head(12)),
            "member_keys": "|".join(m["term_key"]),
            "core_genes": ";".join(gene_counts.index[: 25]),
            "gene_frequency": ";".join(f"{g}:{c}" for g, c in gene_counts.head(40).items()),
            "n_genes": int(len(gene_counts)),
        })
    th = pd.DataFrame(themes).sort_values(["best_tier", "theme_score"], ascending=[True, False])
    th["theme_rank"] = np.arange(1, len(th) + 1)
    th["tier_label"] = th["best_tier"].map({1: "ROBUST", 2: "SUPPORTED", 3: "SINGLE-TOOL"})
    key2theme = {k: remap[c] for k, c in zip(cand["term_key"], cand["_cl"])}
    cons["theme_id"] = cons["term_key"].map(key2theme)
    return cons, th.reset_index(drop=True)


def term_similarity_edges(cons: pd.DataFrame, keys: list[str], cut: float = KAPPA_CUT) -> list[tuple]:
    """Kappa edges among the given terms (for the enrichment-map network)."""
    sub = cons.set_index("term_key").loc[keys]
    genes = sorted({g for h in sub["hit_genes"] for g in h.split(";") if g})
    if not genes or len(keys) < 2:
        return []
    gpos = {g: i for i, g in enumerate(genes)}
    X = np.zeros((len(keys), len(genes)), dtype=np.int8)
    for i, h in enumerate(sub["hit_genes"]):
        for g in h.split(";"):
            if g:
                X[i, gpos[g]] = 1
    K = kappa_matrix(X)
    return [(keys[i], keys[j], float(K[i, j])) for i in range(len(keys)) for j in range(i + 1, len(keys))
            if K[i, j] >= cut]
