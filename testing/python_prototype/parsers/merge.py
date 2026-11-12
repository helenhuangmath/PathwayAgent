"""
Merge multi-tool ORA/GSEA results and compute a consensus score.

Consensus score = Σ_tools [ -log10(FDR_tool) × tool_weight ]
where tool_weight = 1.0 for ORA, 1.3 for GSEA (directional evidence).
The score is then normalised to [0, 1] across all pathways.
Pathway name matching is fuzzy — we strip suffixes/punctuation and
lowercase before grouping so "KEGG_2021: Apoptosis" and
"Reactome: Apoptosis" cluster under "apoptosis".
"""
from __future__ import annotations
import re
import numpy as np
import pandas as pd


_TOOL_WEIGHTS = {
    "Enrichr":          1.0,
    "gProfiler":        1.0,
    "gProfiler_ORA":    1.0,
    "Reactome":         1.0,
    "Reactome_ORA":     1.0,
    "gseapy_Preranked": 1.3,
    "gProfiler_GSEA":   1.2,
}

_STRIP = re.compile(
    r"\(GO:\d+\)|GO:\d+|WP\d+|R-HSA-\d+|hsa\d+|"
    r"_(20\d\d)(_\w+)?[:\s]*|_20\d\d\b|"
    r"^KEGG[_\s]+(20\d\d[_\s]*\w*[:\s]*)?|"
    r"^Reactome[:\s]+|^WikiPathway[s]?[:\s]+"
    r"|^GO_\w+_Process[:\s]+|^MSigDB_Hallmark_\d+[:\s]+",
    re.IGNORECASE,
)


def _normalise_name(name: str) -> str:
    n = re.sub(_STRIP, " ", str(name)).strip()
    n = re.sub(r"\s+", " ", n).lower()
    # remove trailing year patterns
    n = re.sub(r"\s+20\d\d$", "", n)
    return n


def merge_results(raw: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Concatenate all tool DataFrames, adding a canonical pathway_key column."""
    frames = []
    for tool_name, df in raw.items():
        if df is None or df.empty:
            continue
        df = df.copy()
        df["tool"] = df.get("tool", tool_name)
        if "pathway_name" not in df.columns:
            continue
        df["pathway_key"] = df["pathway_name"].apply(_normalise_name)
        frames.append(df)

    if not frames:
        return pd.DataFrame()

    merged = pd.concat(frames, ignore_index=True)
    merged["fdr"]  = pd.to_numeric(merged.get("fdr",  1.0), errors="coerce").fillna(1.0)
    merged["pval"] = pd.to_numeric(merged.get("pval", 1.0), errors="coerce").fillna(1.0)
    merged["overlap"] = pd.to_numeric(merged.get("overlap", 0), errors="coerce").fillna(0)
    return merged


def consensus_score(df: pd.DataFrame, fdr_cutoff: float = 0.05) -> pd.DataFrame:
    """
    Add consensus_score, n_tools, and best_fdr columns.
    Returns a de-duplicated DataFrame (one row per pathway_key × source).
    """
    if df.empty or "pathway_key" not in df.columns:
        return df

    def _score_row(row):
        w = _TOOL_WEIGHTS.get(row.get("tool", ""), 1.0)
        fdr = max(float(row.get("fdr", 1.0)), 1e-300)
        return w * (-np.log10(fdr))

    df = df.copy()
    df["_raw_score"] = df.apply(_score_row, axis=1)

    # Aggregate per canonical pathway: sum of scores, min FDR, total tools
    agg = (
        df.groupby("pathway_key")
        .agg(
            consensus_score=("_raw_score", "sum"),
            fdr=("fdr", "min"),
            pval=("pval", "min"),
            n_tools=("tool", "nunique"),
            best_overlap=("overlap", "max"),
        )
        .reset_index()
    )

    # Normalise scores to [0, 1]
    s_max = agg["consensus_score"].max()
    if s_max > 0:
        agg["consensus_score"] = agg["consensus_score"] / s_max

    # Bring back representative pathway_name and merged gene list
    # (pick the row with lowest FDR per key)
    best = (
        df.sort_values("fdr")
        .drop_duplicates(subset="pathway_key")
        [["pathway_key", "pathway_name", "genes", "source", "tool", "database"]]
    )
    agg = agg.merge(best, on="pathway_key", how="left")

    # Aggregate all genes across tools
    def _merge_genes(x):
        all_g: set = set()
        for v in x:
            if isinstance(v, str) and v:
                all_g.update(re.split(r"[;,/]", v))
        return ";".join(sorted(all_g - {""}))

    gene_agg = (
        df.groupby("pathway_key")["genes"]
        .apply(_merge_genes)
        .reset_index()
        .rename(columns={"genes": "all_genes"})
    )
    agg = agg.merge(gene_agg, on="pathway_key", how="left")

    # Also collect which tools found this pathway
    tools_agg = (
        df.groupby("pathway_key")["tool"]
        .apply(lambda x: "|".join(sorted(set(x.dropna().astype(str)))))
        .reset_index()
        .rename(columns={"tool": "tools_found"})
    )
    agg = agg.merge(tools_agg, on="pathway_key", how="left")

    return agg.sort_values("consensus_score", ascending=False).reset_index(drop=True)
