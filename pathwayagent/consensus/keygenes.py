"""
Key genes for follow-up (design §9.4).

  m_gT   = fraction of significant (tool x member-term) results in theme T whose hits contain g
  Gene_g = sum_T ThemeScore_T * m_gT  * (1 + |z(score_g)|  if scores were given)

DRIVER: removing g from the query makes the theme's best locally tested term lose
significance (leave-one-gene-out hypergeometric test, exact and cheap).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def key_genes(records: pd.DataFrame, cons: pd.DataFrame, themes: pd.DataFrame, query,
              alpha: float = 0.05, interest: set | None = None, top: int = 100) -> pd.DataFrame:
    if themes is None or themes.empty:
        return pd.DataFrame()
    use = themes[themes["best_tier"].isin([1, 2])]
    fallback = use.empty
    if fallback:
        use = themes.head(10)
    interest = {g.upper() for g in (interest or set())}
    sig = records[records["padj"] < alpha]
    scores = pd.Series(query.scores) if query.scores else None
    z = None
    if scores is not None and len(scores) > 2 and scores.std() > 0:
        z = (scores - scores.mean()) / scores.std()

    gene_score: dict[str, float] = {}
    gene_themes: dict[str, list] = {}
    gene_tools: dict[str, set] = {}
    driver_of: dict[str, list] = {}
    for th in use.itertuples():
        keys = th.member_keys.split("|")
        rs = sig[sig["term_key"].isin(keys)]
        if rs.empty:
            continue
        n_res = len(rs)
        counts: dict[str, int] = {}
        for tool, h in zip(rs["tool"], rs["hit_genes"]):
            for g in set(x for x in h.split(";") if x):
                counts[g] = counts.get(g, 0) + 1
                gene_tools.setdefault(g, set()).add(tool)
        for g, c in counts.items():
            m = c / n_res
            w = 1 + abs(float(z.get(g, 0))) if z is not None else 1.0
            gene_score[g] = gene_score.get(g, 0) + th.theme_score * m * w
            gene_themes.setdefault(g, []).append(th.theme_name)
        # driver analysis on the theme's best local-ORA term
        lo = records[(records["tool"] == "local_ora") & (records["term_key"].isin(keys))].sort_values("pval")
        if len(lo):
            r = lo.iloc[0]
            k, K, n, N = int(r.overlap), int(r.term_size), int(r.query_size), int(r.universe_size)
            ratio = (r.padj / r.pval) if r.pval > 0 else 1.0
            if r.padj < alpha:
                p_new = stats.hypergeom.sf(k - 2, N, K, n - 1)
                if min(1.0, p_new * ratio) >= alpha:
                    for g in (x for x in r.hit_genes.split(";") if x):
                        driver_of.setdefault(g, []).append(th.theme_name)

    if not gene_score:
        return pd.DataFrame()
    rows = []
    for g, s in gene_score.items():
        rows.append({
            "gene": g, "gene_score": round(s, 2), "n_themes": len(gene_themes[g]),
            "themes": " | ".join(gene_themes[g][:6]), "n_tools": len(gene_tools.get(g, [])),
            "tools": "|".join(sorted(gene_tools.get(g, []))),
            "input_score": float(scores[g]) if scores is not None and g in scores else np.nan,
            "driver_of": " | ".join(driver_of.get(g, [])),
            "gene_of_interest": g.upper() in interest,
        })
    kg = pd.DataFrame(rows).sort_values("gene_score", ascending=False).head(top).reset_index(drop=True)
    kg.insert(0, "rank", np.arange(1, len(kg) + 1))
    kg["why"] = [
        (f"{'Driver of ' + r.driver_of + '; ' if r.driver_of else ''}"
         f"in {r.n_themes} theme(s) ({r.themes[:120]}); hit in {r.n_tools} tool(s)"
         + (f"; score={r.input_score:.2f}" if not pd.isna(r.input_score) else ""))
        for r in kg.itertuples()]
    kg.attrs["fallback"] = fallback
    return kg
