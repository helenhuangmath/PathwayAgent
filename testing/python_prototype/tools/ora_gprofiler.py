"""ORA via g:Profiler REST API (gprofiler2 / REST fallback)."""
from __future__ import annotations
import requests
import pandas as pd


def run_gprofiler_ora(gene_list, organism, fdr_cutoff, min_genes, max_genes) -> pd.DataFrame:
    """
    Query g:Profiler gost endpoint directly via REST (no R dependency).
    Sources queried: GO:BP, GO:MF, GO:CC, KEGG, REAC, WP, HPA, MIRNA.
    """
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    payload = {
        "organism":  organism,
        "query":     gene_list,
        "sources":   ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"],
        "user_threshold":   fdr_cutoff,
        "significance_threshold_method": "fdr",
        "no_evidences": False,
    }
    try:
        resp = requests.post(url, json=payload, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:
        print(f"  [gProfiler] API call failed: {exc}")
        return pd.DataFrame()

    results = data.get("result", [])
    if not results:
        return pd.DataFrame()

    rows = []
    for r in results:
        overlap = len(r.get("intersections", []) or [])
        if not overlap:
            overlap = r.get("intersection_size", 0)
        term_size = r.get("term_size", 0)
        if term_size > max_genes or overlap < min_genes:
            continue
        inters = r.get("intersections") or []
        # intersections may be a list-of-lists (grouped) or flat list
        flat: list = []
        for item in inters:
            if isinstance(item, list):
                flat.extend(str(x) for x in item)
            else:
                flat.append(str(item))
        genes = ";".join(flat)
        rows.append({
            "pathway_name":   r.get("name", ""),
            "pathway_id":     r.get("native", ""),
            "source":         f"gProfiler:{r.get('source','')}",
            "database":       r.get("source", ""),
            "tool":           "gProfiler",
            "fdr":            r.get("p_value", 1.0),
            "pval":           r.get("p_value", 1.0),
            "overlap":        overlap,
            "term_size":      term_size,
            "genes":          genes,
            "odds_ratio":     None,
            "combined_score": None,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df[df["fdr"] < fdr_cutoff].reset_index(drop=True)
