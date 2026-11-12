"""ORA via Enrichr API (gseapy wrapper). Batches all libraries in parallel."""
from __future__ import annotations
import concurrent.futures
import pandas as pd
import gseapy as gp


def _run_one_library(gene_list: list, lib: str, fdr: float, min_overlap: int) -> pd.DataFrame:
    """Run Enrichr for a single library and return a normalised DataFrame."""
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=lib,
            outdir=None,
            no_plot=True,
            cutoff=1.0,          # collect all; filter later
        )
        df = enr.results.copy()
        if df is None or df.empty:
            return pd.DataFrame()

        # Normalise columns
        df = df.rename(columns={
            "Term":              "pathway_name",
            "Adjusted P-value":  "fdr",
            "P-value":           "pval",
            "Odds Ratio":        "odds_ratio",
            "Combined Score":    "combined_score",
            "Genes":             "genes",
            "Gene_set":          "gene_set",
        })
        df["source"]     = f"Enrichr:{lib}"
        df["tool"]       = "Enrichr"
        df["database"]   = lib
        df["overlap"]    = df["Overlap"].apply(
            lambda x: int(str(x).split("/")[0]) if "/" in str(x) else 0
        )
        df = df[df["overlap"] >= min_overlap]
        df = df[df["fdr"] < fdr]

        keep = ["pathway_name", "fdr", "pval", "odds_ratio",
                "combined_score", "overlap", "genes", "source", "tool", "database"]
        existing = [c for c in keep if c in df.columns]
        return df[existing].reset_index(drop=True)
    except Exception:
        return pd.DataFrame()


def run_enrichr_batch(gene_list, libraries, fdr_cutoff, min_genes) -> pd.DataFrame:
    """Run all Enrichr libraries in parallel and concatenate."""
    frames = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(8, len(libraries))) as pool:
        futures = {pool.submit(_run_one_library, gene_list, lib, fdr_cutoff, min_genes): lib
                   for lib in libraries}
        for fut in concurrent.futures.as_completed(futures):
            df = fut.result()
            if df is not None and not df.empty:
                frames.append(df)

    if not frames:
        return pd.DataFrame()

    combined = pd.concat(frames, ignore_index=True)
    return combined
