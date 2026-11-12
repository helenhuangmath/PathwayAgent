"""Paths, organisms, database registry and defaults shared by every module."""
from __future__ import annotations

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
# Large generated state (gene-set library, run results) lives outside the Git tree when possible.
# Resolution order: environment variable, then a ``large_data/`` folder beside the repository
# (the layout used in development), then a folder inside the repository (ignored by Git).
LARGE_DATA_DIR = Path(os.environ.get("PATHWAYAGENT_DATA", PROJECT_ROOT.parent / "large_data"))


def _default_dir(env: str, name: str) -> Path:
    if os.environ.get(env):
        return Path(os.environ[env])
    if LARGE_DATA_DIR.is_dir():
        return LARGE_DATA_DIR / name
    return PROJECT_ROOT / name


LIBRARY_DIR = _default_dir("PATHWAYAGENT_LIBRARY", "library")
RESULTS_DIR = _default_dir("PATHWAYAGENT_RESULTS", "results")
CACHE_PATH = Path(os.environ.get("PATHWAYAGENT_CACHE", LIBRARY_DIR / "cache.sqlite"))
RSCRIPT = os.environ.get("PATHWAYAGENT_RSCRIPT", "/usr/local/bin/Rscript")

ORGANISMS = {
    "mmusculus": {"taxon": 10090, "label": "Mouse", "kegg": "mmu", "reactome": "Mus musculus",
                  "orgdb": "org.Mm.eg.db", "reactomepa": "mouse", "msigdb_species": "Mm"},
    "hsapiens": {"taxon": 9606, "label": "Human", "kegg": "hsa", "reactome": "Homo sapiens",
                 "orgdb": "org.Hs.eg.db", "reactomepa": "human", "msigdb_species": "Hs"},
}

# Canonical database types. Every adapter tags its records with one of these so the
# consensus knows which tools were *eligible* to report a given term (design §9.1).
DATABASES = {
    "GO_BP": "GO Biological Process",
    "GO_MF": "GO Molecular Function",
    "GO_CC": "GO Cellular Component",
    "KEGG": "KEGG pathways",
    "REACTOME": "Reactome pathways",
    "WIKIPATHWAYS": "WikiPathways",
    "HALLMARK": "MSigDB Hallmark",
    "DISGENET": "DisGeNET (disease genes)",
}
DEFAULT_DATABASES = ["GO_BP", "GO_MF", "GO_CC", "KEGG", "REACTOME", "WIKIPATHWAYS", "HALLMARK", "DISGENET"]

# Method families (design §4). Consensus counts families, not tools.
FAMILIES = {
    "F1_ORA": "Over-representation (hypergeometric / Fisher)",
    "F2_GSEA": "Pre-ranked GSEA",
    "F3_RANK": "Rank-based tests (mHG, Wilcoxon)",
    "F4_NETWORK": "Network / PPI-based",
}

DEFAULTS = {
    "min_size": 10,
    "max_size": 500,
    "alpha": 0.05,
    "deadline_s": 60.0,
    "top_n_plot": 20,
}

PRESETS = {
    "quick": ["local_ora", "local_mhg", "local_gsea", "local_rank"],
    "standard": ["local_ora", "local_mhg", "local_gsea", "local_rank",
                 "enrichr", "gprofiler", "reactome", "string", "clusterprofiler", "fgsea"],
    "everything": ["local_ora", "local_mhg", "local_gsea", "local_rank",
                   "enrichr", "gprofiler", "gprofiler_ordered", "reactome", "string",
                   "clusterprofiler", "clusterprofiler_go", "reactomepa", "fgsea", "metascape"],
}


def library_path(organism: str) -> Path:
    return LIBRARY_DIR / organism
