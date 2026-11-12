"""
Publication-grade figures for PathwayAgent:
  - dotplot  (FDR × overlap, coloured by -log10 FDR)
  - barplot  (consensus_score bar chart)
  - network  (pathway similarity network)
"""
from __future__ import annotations
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import networkx as nx


# ── shared style ────────────────────────────────────────────────────────────
FONT = {"family": "DejaVu Sans", "size": 9}
plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})

# DB → colour map
_DB_PALETTE = {
    "GO:BP":      "#4C72B0",
    "GO:MF":      "#55A868",
    "GO:CC":      "#C44E52",
    "KEGG":       "#DD8452",
    "REAC":       "#8172B2",
    "Reactome":   "#8172B2",
    "Reactome_2022": "#8172B2",
    "WP":         "#937860",
    "MSigDB_Hallmark_2020": "#DA8BC3",
    "DisGeNET":   "#8C8C8C",
    "WikiPathway_2023_Human": "#937860",
    "KEGG_2021_Human": "#DD8452",
    "GO_Biological_Process_2023": "#4C72B0",
    "GO_Molecular_Function_2023": "#55A868",
    "GO_Cellular_Component_2023": "#C44E52",
    "_default":   "#666666",
}


def _db_color(db: str) -> str:
    for key in _DB_PALETTE:
        if key in str(db):
            return _DB_PALETTE[key]
    return _DB_PALETTE["_default"]


def _truncate(s: str, n: int = 52) -> str:
    s = str(s)
    return (s[:n] + "…") if len(s) > n else s


# ── dotplot ─────────────────────────────────────────────────────────────────

def make_dotplot(df: pd.DataFrame, path: str, top_n: int = 25):
    df = df.head(top_n).copy()
    if df.empty:
        return

    df["_mlog10fdr"] = -np.log10(df["fdr"].clip(1e-300, 1))
    df["_size"]      = np.sqrt(df["best_overlap"].clip(1)).clip(2, 20) * 8
    df["_color"]     = df["database"].apply(_db_color)
    df["_label"]     = df["pathway_name"].apply(lambda x: _truncate(x, 52))
    df = df.sort_values("consensus_score", ascending=True)

    fig, ax = plt.subplots(figsize=(9, max(5, len(df) * 0.32 + 1)))
    sc = ax.scatter(
        df["_mlog10fdr"], range(len(df)),
        s=df["_size"],
        c=df["_mlog10fdr"],
        cmap="RdYlGn",
        alpha=0.85,
        edgecolors="grey",
        linewidths=0.4,
    )
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["_label"], fontsize=8)
    ax.set_xlabel("−log₁₀(FDR)", fontsize=10)
    ax.set_title("Top Enriched Pathways — Dot Plot", fontsize=12, fontweight="bold")
    ax.axvline(-np.log10(0.05), color="grey", linestyle="--", lw=0.8, alpha=0.5)
    cb = plt.colorbar(sc, ax=ax, shrink=0.5, pad=0.01)
    cb.set_label("−log₁₀(FDR)", fontsize=8)

    # legend for dot size
    for s_val, label in [(5, "5 genes"), (25, "25 genes"), (100, "100 genes")]:
        ax.scatter([], [], s=np.sqrt(s_val) * 8, c="grey",
                   alpha=0.7, label=label, edgecolors="grey", linewidths=0.4)
    ax.legend(title="Overlap", fontsize=7, title_fontsize=8,
              loc="lower right", frameon=False)

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── barplot ─────────────────────────────────────────────────────────────────

def make_barplot(df: pd.DataFrame, path: str, top_n: int = 25):
    df = df.head(top_n).copy()
    if df.empty:
        return

    df["_label"]  = df["pathway_name"].apply(lambda x: _truncate(x, 55))
    df["_color"]  = df["database"].apply(_db_color)
    df = df.sort_values("consensus_score", ascending=True)

    fig, ax = plt.subplots(figsize=(9, max(5, len(df) * 0.35 + 1)))
    bars = ax.barh(range(len(df)), df["consensus_score"],
                   color=df["_color"], edgecolor="white", height=0.7)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["_label"], fontsize=8)
    ax.set_xlabel("Consensus Score (normalised)", fontsize=10)
    ax.set_title("Top Pathways — Consensus Score", fontsize=12, fontweight="bold")

    # Annotate n_tools
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(row["consensus_score"] + 0.005, i,
                f"n={int(row['n_tools'])}",
                va="center", fontsize=7, color="#444")

    # Database legend
    seen_db = df["database"].unique()
    patches = [mpatches.Patch(color=_db_color(d), label=d[:30])
               for d in seen_db if d]
    ax.legend(handles=patches, fontsize=7, loc="lower right",
              frameon=False, title="Database", title_fontsize=8)

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── network plot ─────────────────────────────────────────────────────────────

def make_network_plot(G: nx.Graph, df: pd.DataFrame, path: str):
    if G is None or G.number_of_nodes() == 0:
        return

    # Remove isolated nodes with no edges if graph is dense
    connected = [n for n in G.nodes() if G.degree(n) > 0]
    sub = G.subgraph(connected) if len(connected) > 3 else G

    if sub.number_of_nodes() == 0:
        sub = G

    fig, ax = plt.subplots(figsize=(12, 10))

    # Layout — try spring with weight, fallback to kamada-kawai
    try:
        pos = nx.spring_layout(sub, weight="weight", seed=42, k=1.5)
    except Exception:
        pos = nx.kamada_kawai_layout(sub)

    # Node properties
    node_scores = [sub.nodes[n].get("consensus_score", 0.1) for n in sub.nodes()]
    node_sizes  = [max(100, s * 3000) for s in node_scores]
    node_colors = [_db_color(sub.nodes[n].get("database", "")) for n in sub.nodes()]
    node_n_tools = [sub.nodes[n].get("n_tools", 1) for n in sub.nodes()]

    # Edge width by Jaccard
    edge_weights = [sub[u][v].get("weight", 0.1) * 5 for u, v in sub.edges()]

    nx.draw_networkx_edges(sub, pos, ax=ax,
                            width=edge_weights, alpha=0.35, edge_color="#aaaaaa")
    nx.draw_networkx_nodes(sub, pos, ax=ax,
                            node_size=node_sizes, node_color=node_colors,
                            alpha=0.85, edgecolors="white", linewidths=0.8)

    # Label only top-20 by score
    top_names = (
        df.sort_values("consensus_score", ascending=False)
        .head(20)["pathway_name"]
        .tolist()
    )
    label_nodes = {n: _truncate(n, 30) for n in sub.nodes() if n in top_names}
    nx.draw_networkx_labels(sub, pos, labels=label_nodes, ax=ax,
                             font_size=6.5, font_color="#111",
                             bbox=dict(boxstyle="round,pad=0.15",
                                       fc="white", alpha=0.6, lw=0))

    ax.set_title("Pathway Similarity Network\n(node size = consensus score, "
                 "edge weight = Jaccard overlap)",
                 fontsize=12, fontweight="bold")
    ax.axis("off")

    # Legend patches
    seen_db = list({sub.nodes[n].get("database", "") for n in sub.nodes()})
    patches = [mpatches.Patch(color=_db_color(d), label=d[:28])
               for d in sorted(seen_db) if d]
    if patches:
        ax.legend(handles=patches, fontsize=8, loc="lower left",
                  frameon=True, title="Database")

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
