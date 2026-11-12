"""
Build a pathway-gene bipartite network and a pathway-pathway similarity
network (Jaccard overlap of gene sets). Returns a networkx Graph.
"""
from __future__ import annotations
import re
import networkx as nx
import pandas as pd
import numpy as np


def _parse_genes(gene_str: str) -> set:
    if not isinstance(gene_str, str) or not gene_str:
        return set()
    return set(re.split(r"[;,/\s]+", gene_str.strip())) - {""}


def build_pathway_network(
    pathways: pd.DataFrame,
    query_genes: list[str],
    jaccard_threshold: float = 0.1,
    max_pathways: int = 50,
) -> nx.Graph:
    """
    Build a pathway-similarity network (pathway nodes connected by shared genes).

    Nodes = pathways, coloured by database.
    Edges = Jaccard index ≥ jaccard_threshold.
    Node size = consensus_score.
    """
    df = pathways.head(max_pathways).copy()
    query_set = set(g.upper() for g in query_genes)

    G = nx.Graph()

    # ── add pathway nodes ──────────────────────────────────────────────
    gene_sets: dict[str, set] = {}
    for _, row in df.iterrows():
        name  = str(row.get("pathway_name", row.get("pathway_key", "unknown")))
        genes = _parse_genes(str(row.get("all_genes", row.get("genes", ""))))
        # intersect with query set so edges only reflect direct evidence
        genes = genes & query_set if query_set else genes
        gene_sets[name] = genes

        G.add_node(
            name,
            consensus_score=float(row.get("consensus_score", 0)),
            fdr=float(row.get("fdr", 1.0)),
            n_tools=int(row.get("n_tools", 1)),
            database=str(row.get("database", "")),
            genes=";".join(genes),
        )

    # ── add edges by Jaccard ───────────────────────────────────────────
    names = list(gene_sets.keys())
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = gene_sets[names[i]], gene_sets[names[j]]
            union = a | b
            if not union:
                continue
            j_idx = len(a & b) / len(union)
            if j_idx >= jaccard_threshold:
                G.add_edge(names[i], names[j],
                           weight=j_idx,
                           shared_genes=";".join(a & b))

    return G
