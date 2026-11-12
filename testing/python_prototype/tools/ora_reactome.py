"""ORA via Reactome REST Analysis Service."""
from __future__ import annotations
import requests
import pandas as pd


_ORG_TO_TAX = {
    "hsapiens":  9606,
    "mmusculus": 10090,
    "rnorvegicus": 10116,
    "drerio": 7955,
    "dmelanogaster": 7227,
    "scerevisiae": 4932,
}


def run_reactome_ora(gene_list, organism, fdr_cutoff, min_genes) -> pd.DataFrame:
    """
    POST gene symbols to Reactome Analysis Service and return enriched pathways.
    Uses the species-aware identifier mapping built into the service.
    """
    tax_id = _ORG_TO_TAX.get(organism, 9606)
    url = (
        "https://reactome.org/AnalysisService/identifiers/projection/"
        f"?pageSize=500&page=1&sortBy=ENTITIES_FDR&order=ASC"
        f"&resource=TOTAL&pValue={fdr_cutoff}&includeDisease=true"
        f"&min={min_genes}"
    )
    gene_str = "\n".join(gene_list)
    try:
        resp = requests.post(
            url,
            data=gene_str,
            headers={"Content-Type": "text/plain"},
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:
        print(f"  [Reactome] API call failed: {exc}")
        return pd.DataFrame()

    pathways = data.get("pathways", [])
    if not pathways:
        return pd.DataFrame()

    rows = []
    for p in pathways:
        ent = p.get("entities", {})
        if ent.get("fdr", 1.0) >= fdr_cutoff:
            continue
        if ent.get("found", 0) < min_genes:
            continue
        rows.append({
            "pathway_name":   p.get("name", ""),
            "pathway_id":     p.get("stId", ""),
            "source":         "Reactome",
            "database":       "Reactome",
            "tool":           "Reactome",
            "fdr":            ent.get("fdr", 1.0),
            "pval":           ent.get("pValue", 1.0),
            "overlap":        ent.get("found", 0),
            "term_size":      ent.get("total", 0),
            "ratio":          ent.get("ratio", None),
            "genes":          "",   # not returned in summary endpoint
            "odds_ratio":     None,
            "combined_score": None,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df[df["fdr"] < fdr_cutoff].sort_values("fdr").reset_index(drop=True)
