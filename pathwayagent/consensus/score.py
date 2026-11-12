"""
Term-level consensus scoring (design §9.1-9.2).

Per canonical term t and per method family f, over the tools *eligible* for t
(the tool ran and covers t's database):

  sig_f = median_k  min(10, -log10 padj_k) / 10        (padj_k = 1 if k did not report t)
  sup_f = mean_k    [padj_k < alpha]
  S_sig = mean_f sig_f     S_sup = mean_f sup_f        (family-balanced: 8 ORA tools != 8 votes)
  S_eff = median over reporting tools of the within-tool percentile of the effect size
  D     = direction consistency of significant directional results (1 if none)
  spec  = 1 / (1 + log10(size/200)) for size > 200, x0.7 for very general GO terms (depth <= 3)

  ConsensusScore = 100 * spec * (0.40 S_sig + 0.30 S_sup + 0.15 S_eff + 0.15 D)

Rank consistency is reported separately as a Robust Rank Aggregation p-value
(Kolde et al. 2012) over the tools' normalised ranks.
Tiers: 1 ROBUST, 2 SUPPORTED, 3 SINGLE-TOOL, 0 not significant.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps

from ..core.stats import bh

WEIGHTS = {"sig": 0.40, "sup": 0.30, "eff": 0.15, "dir": 0.15}
TIER_LABEL = {1: "ROBUST", 2: "SUPPORTED", 3: "SINGLE-TOOL", 0: "n.s."}


def rra_pvalue(r: np.ndarray) -> float:
    """RRA rho-score p-value for a vector of normalised ranks in (0, 1]."""
    r = np.sort(np.clip(r, 1e-12, 1.0))
    n = len(r)
    if n == 0:
        return 1.0
    j = np.arange(1, n + 1)
    rho = sps.beta.cdf(r, j, n - j + 1).min()
    return float(min(1.0, rho * n))


def _specificity(size: float, depth: float | None) -> float:
    s = 1.0
    if size and size > 200:
        s = 1.0 / (1.0 + np.log10(size / 200.0))
    if depth is not None and not pd.isna(depth) and depth <= 3:
        s *= 0.7
    return s


def consensus(records: pd.DataFrame, tool_info: dict, library, alpha: float = 0.05,
              weights: dict = WEIGHTS, drop_tool: str | None = None) -> pd.DataFrame:
    """
    records   : harmonised ResultRecords for ONE query (all tools)
    tool_info : tool -> {"family": str, "databases": [..]} for tools that ran (ok/empty)
    """
    tools = {t: v for t, v in tool_info.items() if t != drop_tool}
    rec = records[records["tool"].isin(tools)]
    if rec.empty:
        return pd.DataFrame()

    # best record per (term, tool)
    rec = rec.sort_values("padj").drop_duplicates(["term_key", "tool"])
    padj = rec.pivot(index="term_key", columns="tool", values="padj")
    terms = padj.index
    tool_list = list(padj.columns.union(pd.Index(list(tools))))
    padj = padj.reindex(columns=tool_list)

    # database of each term (majority over reports), library metadata when available
    meta = (rec.groupby("term_key")
            .agg(database=("database", lambda x: x.mode().iat[0]),
                 reported_name=("term_name", lambda x: x.mode().iat[0]),
                 reported_size=("term_size", "median"))
            .reindex(terms))
    lib_rows = [library.key_pos.get(k) for k in terms]
    meta["term_name"] = [library.sets["name"].iat[i] if i is not None else n
                         for i, n in zip(lib_rows, meta["reported_name"])]
    meta["term_size"] = [library.sets["size"].iat[i] if i is not None else s
                         for i, s in zip(lib_rows, meta["reported_size"])]
    meta["go_depth"] = [library.sets["go_depth"].iat[i] if i is not None else np.nan for i in lib_rows]
    meta["url"] = [library.sets["url"].iat[i] if i is not None else "" for i in lib_rows]
    meta["in_library"] = [i is not None for i in lib_rows]

    # eligibility: tool ran and covers the term's database
    elig = pd.DataFrame(False, index=terms, columns=tool_list)
    for t in tool_list:
        elig[t] = meta["database"].isin(tools[t]["databases"]).values
    elig = elig | padj.notna()              # a reported term is always eligible
    p_full = padj.fillna(1.0).where(elig)   # NaN = not eligible
    mlog = np.minimum(10.0, -np.log10(p_full.clip(lower=1e-300))) / 10.0
    sig = (p_full < alpha).astype(float).where(elig)

    families = {}
    for t in tool_list:
        families.setdefault(tools[t]["family"], []).append(t)
    fam_sig, fam_sup, fam_any = [], [], []
    for fam, ts in families.items():
        fam_sig.append(mlog[ts].median(axis=1, skipna=True))
        fam_sup.append(sig[ts].mean(axis=1, skipna=True))
        fam_any.append((sig[ts] > 0).any(axis=1) & elig[ts].any(axis=1))
    S_sig = pd.concat(fam_sig, axis=1).mean(axis=1, skipna=True).fillna(0)
    S_sup = pd.concat(fam_sup, axis=1).mean(axis=1, skipna=True).fillna(0)
    n_fam_elig = pd.concat([elig[ts].any(axis=1) for ts in families.values()], axis=1).sum(axis=1)
    n_fam_sig = pd.concat(fam_any, axis=1).sum(axis=1)

    # effect: within-tool percentile of effect size (|NES| for GSEA), median over reporting tools
    rr = rec.copy()
    rr["eff_abs"] = rr["effect"].abs()
    rr["eff_pct"] = rr.groupby("tool")["eff_abs"].rank(pct=True)
    S_eff = rr.groupby("term_key")["eff_pct"].median().reindex(terms).fillna(0)

    # direction consistency among significant directional results
    dr = rec[(rec["direction"] != 0) & (rec["padj"] < alpha)]
    dsum = dr.groupby("term_key")["direction"].agg(["sum", "count"]).reindex(terms)
    D = (dsum["sum"].abs() / dsum["count"]).fillna(1.0)
    direction = np.sign(dsum["sum"].fillna(0)).astype(int)

    spec = pd.Series([_specificity(s, d) for s, d in zip(meta["term_size"], meta["go_depth"])], index=terms)
    score = 100 * spec * (weights["sig"] * S_sig + weights["sup"] * S_sup +
                          weights["eff"] * S_eff + weights["dir"] * D)

    # rank consistency (RRA) over eligible tools: normalised rank, unreported -> 1
    n_per_db = library.sets.groupby("database").size()
    nr = rec.assign(nrank=rec["rank"] / rec["database"].map(n_per_db).fillna(rec["n_tested"]).clip(lower=1))
    nrank = nr.pivot(index="term_key", columns="tool", values="nrank").reindex(index=terms, columns=tool_list)
    nrank = nrank.clip(upper=1.0).fillna(1.0).where(elig)
    rra_p = np.array([rra_pvalue(row[~np.isnan(row)]) for row in nrank.values])

    n_elig = elig.sum(axis=1)
    n_sig = (sig > 0).sum(axis=1)
    best_p = padj.min(axis=1)
    tier = np.zeros(len(terms), dtype=int)
    multi_fam_ok = (n_fam_elig <= 1) | (n_fam_sig >= 2)
    t1 = (n_sig >= 2) & (S_sup >= 0.6) & (S_sig >= 0.2) & (D >= 0.8) & multi_fam_ok
    t2 = (n_sig >= 2) & (S_sup >= 0.3)
    t3 = n_sig >= 1
    tier[t3.values] = 3
    tier[t2.values] = 2
    tier[t1.values] = 1

    sig_tools = sig.apply(lambda r: "|".join(sorted(c for c in r.index if r[c] == 1.0)), axis=1)
    hits = (rec[rec["padj"] < alpha].assign(g=lambda d: d["hit_genes"].str.split(";"))
            .groupby("term_key")["g"].apply(lambda s: sorted({x for l in s for x in l if x})))
    hits_any = (rec.assign(g=lambda d: d["hit_genes"].str.split(";"))
                .groupby("term_key")["g"].apply(lambda s: sorted({x for l in s for x in l if x})))
    hit_genes = [";".join(hits.get(k, hits_any.get(k, []))) for k in terms]

    out = pd.DataFrame({
        "term_key": terms, "term_name": meta["term_name"].values, "database": meta["database"].values,
        "term_size": meta["term_size"].values, "consensus_score": score.values.round(2),
        "tier": tier, "tier_label": [TIER_LABEL[t] for t in tier],
        "n_tools_sig": n_sig.values, "n_tools_eligible": n_elig.values,
        "tool_support": [f"{a}/{b}" for a, b in zip(n_sig.values, n_elig.values)],
        "n_families_sig": n_fam_sig.values, "n_families_eligible": n_fam_elig.values,
        "S_sig": S_sig.values.round(3), "S_sup": S_sup.values.round(3), "S_eff": S_eff.values.round(3),
        "direction_consistency": D.values.round(3), "direction": direction.values,
        "specificity": spec.values.round(3), "best_padj": best_p.values,
        "median_padj": padj.fillna(1.0).where(elig).median(axis=1).values,
        "rra_p": rra_p, "tools_significant": sig_tools.values, "hit_genes": hit_genes,
        "n_hit_genes": [len(h.split(";")) if h else 0 for h in hit_genes],
        "url": meta["url"].values, "in_library": meta["in_library"].values,
    })
    out["rra_q"] = bh(out["rra_p"].values)
    for t in tool_list:  # per-tool padj columns for the heatmap / tables
        out[f"padj__{t}"] = padj[t].values
    out["_tier_order"] = out["tier"].replace(0, 9)
    out = (out.sort_values(["_tier_order", "consensus_score", "best_padj"], ascending=[True, False, True])
           .drop(columns="_tier_order").reset_index(drop=True))
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    return out


def stability(records: pd.DataFrame, tool_info: dict, library, base: pd.DataFrame,
              alpha: float = 0.05, top: int = 50) -> pd.DataFrame:
    """Leave-one-tool-out: rank range of the top terms and which tool they depend on."""
    top_keys = base.loc[base["tier"] > 0, "term_key"].head(top).tolist()
    if len(tool_info) < 2 or not top_keys:
        base["loo_rank_range"] = ""
        base["depends_on"] = ""
        return base
    sub = records[records["term_key"].isin(top_keys)]
    ranks = {k: [] for k in top_keys}
    lost_by = {k: [] for k in top_keys}
    for t in tool_info:
        c = consensus(sub, tool_info, library, alpha, drop_tool=t)
        if c.empty:
            continue
        pos = {k: i + 1 for i, k in enumerate(c["term_key"])}
        tiers = dict(zip(c["term_key"], c["tier"]))
        for k in top_keys:
            ranks[k].append(pos.get(k, len(top_keys) + 1))
            if tiers.get(k, 0) == 0:
                lost_by[k].append(t)
    base["loo_rank_range"] = [f"{min(ranks[k])}-{max(ranks[k])}" if k in ranks and ranks[k] else ""
                              for k in base["term_key"]]
    base["depends_on"] = ["|".join(lost_by.get(k, [])) for k in base["term_key"]]
    return base
