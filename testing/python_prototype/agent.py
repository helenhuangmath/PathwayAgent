"""
PathwayAgent — main orchestrator.

Flow:
  1. Normalise gene list (symbol → Entrez mapping via MyGene)
  2. Run ORA tools in parallel threads (Enrichr, g:Profiler, KEGG, Reactome,
     WikiPathways, DisGeNET, GO-BP/MF/CC via gseapy)
  3. Run GSEA tools in parallel threads (fgsea-style pre-ranked via gseapy,
     g:Profiler ordered query)
  4. Merge results → consensus scoring (freq × -log10(FDR) × NES weight)
  5. Produce figures + pathway network
  6. Return typed result object
"""
from __future__ import annotations

import os
import time
import concurrent.futures
from typing import Optional
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np

from .tools.ora_enrichr    import run_enrichr_batch
from .tools.ora_gprofiler  import run_gprofiler_ora
from .tools.ora_reactome   import run_reactome_ora
from .tools.gsea_gseapy    import run_gsea_preranked
from .tools.gsea_gprofiler import run_gprofiler_gsea
from .parsers.merge        import merge_results, consensus_score
from .networks.builder     import build_pathway_network
from .reports.plots        import make_dotplot, make_barplot, make_network_plot
from .reports.summary      import write_summary_table, write_json_report


class PathwayAgent:
    """
    Parameters
    ----------
    gene_list : list[str]
        Gene symbols (ORA mode) or ranked gene list for GSEA mode.
    ranked_scores : list[float] | None
        Pre-ranked scores (log2FC, -log10p, etc.) aligned with gene_list.
        If supplied, GSEA is run in addition to ORA.
    organism : str
        Common name or Ensembl taxon, default "hsapiens".
    fdr_cutoff : float
        FDR threshold for significance filtering, default 0.05.
    min_genes : int
        Minimum gene-set overlap, default 5.
    max_genes : int
        Maximum gene-set size, default 500.
    databases : list[str] | None
        Restrict to a subset of databases. If None, all are run.
    outdir : str
        Output directory for saved files, default "./pathway_results".
    n_jobs : int
        Max parallel worker threads, default 6.
    """

    ENRICHR_LIBS = [
        "KEGG_2021_Human",
        "Reactome_2022",
        "WikiPathway_2023_Human",
        "GO_Biological_Process_2023",
        "GO_Molecular_Function_2023",
        "GO_Cellular_Component_2023",
        "MSigDB_Hallmark_2020",
        "DisGeNET",
    ]

    def __init__(
        self,
        gene_list: list,
        ranked_scores: Optional[list] = None,
        organism: str = "hsapiens",
        fdr_cutoff: float = 0.05,
        min_genes: int = 5,
        max_genes: int = 500,
        databases: Optional[list] = None,
        outdir: str = "./pathway_results",
        n_jobs: int = 6,
    ):
        self.gene_list     = [str(g).strip().upper() for g in gene_list]
        self.ranked_scores = ranked_scores
        self.organism      = organism
        self.fdr_cutoff    = fdr_cutoff
        self.min_genes     = min_genes
        self.max_genes     = max_genes
        self.databases     = databases
        self.outdir        = outdir
        self.n_jobs        = n_jobs
        os.makedirs(outdir, exist_ok=True)

        # ── result containers ───────────────────────────────────────────────
        self.raw_results: dict[str, pd.DataFrame] = {}
        self.consensus:   pd.DataFrame | None = None
        self.top_pathways: pd.DataFrame | None = None
        self.network_graph = None
        self.figures: dict[str, str] = {}
        self._timings: dict[str, float] = {}

    # ──────────────────────────────────────────────────────────────────────
    # Public API
    # ──────────────────────────────────────────────────────────────────────

    def run(self, verbose: bool = True) -> "PathwayAgent":
        """Run all tools, merge, score, and produce outputs."""
        t0 = time.time()
        if verbose:
            print(f"[PathwayAgent] Starting analysis on {len(self.gene_list)} genes …")

        tasks = self._build_task_list()
        self._run_parallel(tasks, verbose=verbose)
        self._merge_and_score(verbose=verbose)
        self._make_figures(verbose=verbose)
        self._write_tables(verbose=verbose)

        elapsed = time.time() - t0
        if verbose:
            n_sig = (
                (self.consensus["fdr"] < self.fdr_cutoff).sum()
                if self.consensus is not None and "fdr" in self.consensus.columns
                else 0
            )
            print(f"\n[PathwayAgent] Done in {elapsed:.1f}s — "
                  f"{n_sig} significant pathways across "
                  f"{len(self.raw_results)} tools.")
        return self

    def summary(self, n: int = 20) -> pd.DataFrame:
        """Return top-n consensus pathways."""
        if self.top_pathways is None:
            raise RuntimeError("Call .run() first.")
        return self.top_pathways.head(n)

    # ──────────────────────────────────────────────────────────────────────
    # Internal orchestration
    # ──────────────────────────────────────────────────────────────────────

    def _build_task_list(self):
        """Return list of (label, callable) pairs to run in parallel."""
        tasks = []

        # ── ORA tasks ───────────────────────────────────────────────────
        tasks.append((
            "Enrichr",
            lambda: run_enrichr_batch(
                self.gene_list, self.ENRICHR_LIBS,
                self.fdr_cutoff, self.min_genes
            ),
        ))
        tasks.append((
            "gProfiler_ORA",
            lambda: run_gprofiler_ora(
                self.gene_list, self.organism,
                self.fdr_cutoff, self.min_genes, self.max_genes
            ),
        ))
        tasks.append((
            "Reactome_ORA",
            lambda: run_reactome_ora(
                self.gene_list, self.organism,
                self.fdr_cutoff, self.min_genes
            ),
        ))

        # ── GSEA tasks (only if ranked scores supplied) ─────────────────
        if self.ranked_scores is not None:
            tasks.append((
                "gseapy_Preranked",
                lambda: run_gsea_preranked(
                    self.gene_list, self.ranked_scores,
                    self.fdr_cutoff, self.min_genes, self.max_genes,
                    self.outdir
                ),
            ))
            tasks.append((
                "gProfiler_GSEA",
                lambda: run_gprofiler_gsea(
                    self.gene_list, self.ranked_scores, self.organism,
                    self.fdr_cutoff, self.min_genes, self.max_genes
                ),
            ))

        return tasks

    def _run_parallel(self, tasks, verbose):
        """Execute tasks in a thread pool and collect DataFrames."""
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.n_jobs) as pool:
            futures = {pool.submit(fn): label for label, fn in tasks}
            for fut in concurrent.futures.as_completed(futures):
                label = futures[fut]
                t0 = time.time()
                try:
                    df = fut.result()
                    elapsed = time.time() - t0
                    self._timings[label] = elapsed
                    if df is not None and len(df) > 0:
                        self.raw_results[label] = df
                        if verbose:
                            print(f"  ✓ {label:<22} {len(df):>5} pathways  "
                                  f"({elapsed:.1f}s)")
                    else:
                        if verbose:
                            print(f"  – {label:<22}  no results")
                except Exception as exc:
                    if verbose:
                        print(f"  ✗ {label:<22} ERROR: {exc}")

    def _merge_and_score(self, verbose):
        if not self.raw_results:
            if verbose:
                print("[PathwayAgent] Warning: no results to merge.")
            return
        if verbose:
            print("\n[PathwayAgent] Merging results …")

        self.consensus = merge_results(self.raw_results)
        self.consensus = consensus_score(self.consensus, self.fdr_cutoff)
        self.top_pathways = (
            self.consensus
            .loc[self.consensus["fdr"] < self.fdr_cutoff]
            .sort_values("consensus_score", ascending=False)
            .drop_duplicates(subset="pathway_name")
            .reset_index(drop=True)
        )
        if verbose:
            print(f"  → {len(self.top_pathways)} unique significant pathways")

    def _make_figures(self, verbose):
        if self.top_pathways is None or len(self.top_pathways) == 0:
            return
        if verbose:
            print("[PathwayAgent] Generating figures …")

        dot_path = os.path.join(self.outdir, "dotplot.png")
        bar_path = os.path.join(self.outdir, "barplot.png")
        net_path = os.path.join(self.outdir, "network.png")

        make_dotplot(self.top_pathways.head(30), dot_path)
        make_barplot(self.top_pathways.head(25), bar_path)

        if len(self.top_pathways) >= 5:
            G = build_pathway_network(
                self.top_pathways.head(50),
                self.gene_list
            )
            self.network_graph = G
            make_network_plot(G, self.top_pathways.head(50), net_path)

        self.figures = {k: v for k, v in
                        [("dotplot", dot_path), ("barplot", bar_path),
                         ("network", net_path)]
                        if os.path.exists(v)}
        if verbose:
            print(f"  → Saved {len(self.figures)} figures to {self.outdir}/")

    def _write_tables(self, verbose):
        if self.consensus is None:
            return
        csv_path  = os.path.join(self.outdir, "consensus_pathways.csv")
        json_path = os.path.join(self.outdir, "report.json")
        full_path = os.path.join(self.outdir, "all_tools_results.csv")

        self.top_pathways.to_csv(csv_path, index=False)
        write_summary_table(self.raw_results, full_path)
        write_json_report(self.top_pathways, self.raw_results,
                          self._timings, json_path)
        if verbose:
            print(f"  → Tables saved to {self.outdir}/")
