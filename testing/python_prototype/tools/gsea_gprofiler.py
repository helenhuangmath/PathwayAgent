"""Ordered GSEA via g:Profiler (gprofiler ordered query mode)."""
from __future__ import annotations
import requests
import pandas as pd
import numpy as np


def run_gprofiler_gsea(gene_list, scores, organism, fdr_cutoff,
                        min_genes, max_genes) -> pd.DataFrame:
    """
    Run g:Profiler in ordered-query mode (GSEA-like).
    Genes must be ordered by descending significance/score.
    """
    scores_arr = np.array(scores, dtype=float)
    order = np.argsort(scores_arr)[::-1]
    ordered_genes = [gene_list[i] for i in order]

    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    payload = {
        "organism":   organism,
        "query":      ordered_genes,
        "ordered":    True,
        "sources":    ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"],
        "user_threshold":   fdr_cutoff,
        "significance_threshold_method": "fdr",
    }
    try:
        resp = requests.post(url, json=payload, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:
        print(f"  [gProfiler-GSEA] API call failed: {exc}")
        return pd.DataFrame()

    results = data.get("result", [])
    if not results:
        return pd.DataFrame()

    rows = []
    for r in results:
        overlap = r.get("intersection_size", 0)
        term_size = r.get("term_size", 0)
        if term_size > max_genes or overlap < min_genes:
            continue
        inters = r.get("intersections") or []
        flat_g: list = []
        for item in inters:
            if isinstance(item, list):
                flat_g.extend(str(x) for x in item)
            else:
                flat_g.append(str(item))
        rows.append({
            "pathway_name":   r.get("name", ""),
            "pathway_id":     r.get("native", ""),
            "source":         f"gProfiler_GSEA:{r.get('source','')}",
            "database":       r.get("source", ""),
            "tool":           "gProfiler_GSEA",
            "fdr":            r.get("p_value", 1.0),
            "pval":           r.get("p_value", 1.0),
            "overlap":        overlap,
            "term_size":      term_size,
            "genes":          ";".join(flat_g),
            "odds_ratio":     None,
            "combined_score": None,
            "nes":            None,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df[df["fdr"] < fdr_cutoff].reset_index(drop=True)
