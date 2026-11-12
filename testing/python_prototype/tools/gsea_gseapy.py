"""Pre-ranked GSEA via gseapy across MSigDB, KEGG, GO, Reactome."""
from __future__ import annotations
import os
import concurrent.futures
import pandas as pd
import numpy as np
import gseapy as gp


_GSEA_LIBS = [
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human",
    "Reactome_2022",
    "GO_Biological_Process_2023",
    "WikiPathway_2023_Human",
]


def _gsea_one_library(rnk: pd.Series, lib: str, min_size: int,
                       max_size: int, outdir: str) -> pd.DataFrame:
    try:
        res = gp.prerank(
            rnk=rnk,
            gene_sets=lib,
            min_size=min_size,
            max_size=max_size,
            permutation_num=500,
            outdir=outdir,
            seed=42,
            no_plot=True,
            verbose=False,
        )
        df = res.res2d.copy()
        if df is None or df.empty:
            return pd.DataFrame()

        df = df.rename(columns={
            "Term":  "pathway_name",
            "NOM p-val": "pval",
            "FDR q-val": "fdr",
            "NES":  "nes",
            "ES":   "es",
            "Lead_genes": "genes",
        })
        df["source"]   = f"GSEA:{lib}"
        df["database"] = lib
        df["tool"]     = "gseapy_Preranked"
        df["overlap"]  = df.get("Tag %", "").apply(
            lambda x: int(x.split("%")[0]) if isinstance(x, str) and "%" in x else 0
        )
        df["odds_ratio"]     = None
        df["combined_score"] = df.get("nes", 0.0)

        keep = ["pathway_name", "fdr", "pval", "nes", "es", "overlap",
                "genes", "source", "tool", "database",
                "odds_ratio", "combined_score"]
        existing = [c for c in keep if c in df.columns]
        df = df[existing]
        df["fdr"]  = pd.to_numeric(df["fdr"],  errors="coerce").fillna(1.0)
        df["pval"] = pd.to_numeric(df["pval"], errors="coerce").fillna(1.0)
        return df.reset_index(drop=True)
    except Exception:
        return pd.DataFrame()


def run_gsea_preranked(gene_list, scores, fdr_cutoff,
                        min_genes, max_genes, outdir) -> pd.DataFrame:
    """Run pre-ranked GSEA on all libraries in parallel."""
    scores_arr = np.array(scores, dtype=float)
    rnk = pd.Series(scores_arr, index=gene_list).sort_values(ascending=False)
    # Drop duplicates — keep highest score
    rnk = rnk[~rnk.index.duplicated(keep="first")]

    gsea_out = os.path.join(outdir, "gsea_tmp")
    os.makedirs(gsea_out, exist_ok=True)

    frames = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures = {
            pool.submit(_gsea_one_library, rnk, lib, min_genes, max_genes, gsea_out): lib
            for lib in _GSEA_LIBS
        }
        for fut in concurrent.futures.as_completed(futures):
            df = fut.result()
            if df is not None and not df.empty:
                frames.append(df)

    if not frames:
        return pd.DataFrame()

    combined = pd.concat(frames, ignore_index=True)
    combined = combined[combined["fdr"] < fdr_cutoff]
    return combined.reset_index(drop=True)
