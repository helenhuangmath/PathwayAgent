"""
PathwayAgent — A fast, multi-tool pathway analysis orchestrator.

Runs ORA and GSEA across: gseapy (Enrichr API), g:Profiler, Reactome,
KEGG, WikiPathways, MSigDB, GO (BP/MF/CC). Produces per-tool tables,
ranked cross-tool consensus, dot/bar plots, and a pathway network figure.
Designed for < 60 s total runtime on a standard gene list.
"""
from .agent import PathwayAgent

__all__ = ["PathwayAgent"]
