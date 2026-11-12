"""Adapter registry: tool name -> adapter class."""
from .local import LocalGSEA, LocalMHG, LocalORA, LocalRank
from .metascape import Metascape
from .rtools import FGSEA, ClusterProfiler, ClusterProfilerGO, ReactomePA
from .web import Enrichr, GProfiler, GProfilerOrdered, Reactome, StringDB

REGISTRY = {cls.name: cls for cls in (
    LocalORA, LocalMHG, LocalRank, LocalGSEA,
    Enrichr, GProfiler, GProfilerOrdered, Reactome, StringDB,
    ClusterProfiler, ClusterProfilerGO, ReactomePA, FGSEA, Metascape,
)}

# Tools that are not implemented yet (design §6) — listed so the UI/report can say so honestly.
PLANNED = {
    "webgestalt": "WebGestaltR not installed in R (install.packages('WebGestaltR')); ORA/GSEA/NTA planned",
    "david": "DAVID web service requires registration; planned (optional, rate-limited)",
    "great": "needs genomic regions (BED) input + rGREAT; planned",
}
