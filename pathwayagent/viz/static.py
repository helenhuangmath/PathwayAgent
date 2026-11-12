"""
Publication figures (design §10).

Every figure is drawn at its final printed size (<= 7.2 in = Nature double column) with
9-10 pt Arial text, and written as a vector PDF (editable TrueType text) plus a PNG preview for
the HTML report. Legends are stacked in a single right-hand column.

Colour roles (journal palettes, fixed order; categorical colours always sit next to a text label):
  categorical (themes)      -> CAT (NPG / NEJM / AAAS hues, order validated for colour-vision separation)
  significance (-log10 P)   -> single-hue red ramp
  signed effects (NES)      -> diverging blue-white-red (red = up)
  membership / overlap      -> single-hue navy ramp
"""
from __future__ import annotations

import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402
from matplotlib.text import Text  # noqa: E402

from ..consensus.themes import kappa_matrix  # noqa: E402

# ───────────────────────────── style ─────────────────────────────
FULL_W = 7.2          # Nature double column (183 mm)
CAT = ["#E64B35", "#3C5488", "#00A087", "#E18727", "#4DBBD5", "#631879", "#7876B1", "#7E6148",
       "#EE4C97", "#7F9A3A"]
INK, INK2, GRID, EDGE_BLUE, NA_GREY = "#000000", "#4D4D4D", "#E6E6E6", "#9ECAE1", "#BDBDBD"
REDS = LinearSegmentedColormap.from_list("npg_reds", ["#FDE0D7", "#F7A48F", "#E64B35", "#B2182B", "#67000D"])
NAVY = LinearSegmentedColormap.from_list("npg_navy", ["#FFFFFF", "#C6D0E6", "#8491B4", "#3C5488", "#1B2A4A"])
RDBU = LinearSegmentedColormap.from_list("rdbu", ["#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7",
                                                  "#FDDBC7", "#EF8A62", "#B2182B"])
TIER_NAMES = {1: "ROBUST", 2: "SUPPORTED", 3: "SINGLE-TOOL"}
TIER_EDGE = {1: (INK, 1.5, "solid"), 2: (INK, 0.7, "solid"), 3: ("#7F7F7F", 0.7, (0, (2, 1.5)))}
TOOL_LABEL = {"local_ora": "Local ORA", "local_mhg": "Local mHG", "local_gsea": "Local GSEA",
              "local_rank": "Local Wilcoxon", "enrichr": "Enrichr", "gprofiler": "g:Profiler",
              "gprofiler_ordered": "g:Profiler ordered", "reactome": "Reactome", "string": "STRING",
              "clusterprofiler": "clusterProfiler", "clusterprofiler_go": "clusterProfiler GO",
              "reactomepa": "ReactomePA", "fgsea": "fgsea", "metascape": "Metascape"}
DB_LABEL = {"GO_BP": "GO Biological Process", "GO_MF": "GO Molecular Function", "GO_CC": "GO Cellular Component",
            "KEGG": "KEGG", "REACTOME": "Reactome", "WIKIPATHWAYS": "WikiPathways", "HALLMARK": "MSigDB Hallmark",
            "DISGENET": "DisGeNET"}
MLOG_CAP = 8.0        # shared range for every -log10 P encoding

plt.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9.5,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8.5, "legend.title_fontsize": 9,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK, "ytick.color": INK,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "mathtext.sf": "Arial", "mathtext.cal": "Arial", "mathtext.default": "regular",
    "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.titlepad": 6,
    "savefig.pad_inches": 0.0,
})
MLOG_LABEL = "−log$_{10}$ adj. P"


# ───────────────────────────── text helpers ─────────────────────────────

def cluster_style(c: int) -> tuple[str, str]:
    """(colour, marker) of map cluster c: fixed by the cluster id, so identical across runs and lists."""
    return CAT[c % len(CAT)], "os^D"[(c // len(CAT)) % 4]


def tool_label(t: str) -> str:
    return TOOL_LABEL.get(t, t)


_ACRONYMS = {"MAPK", "TNF", "NF-KB", "JAK", "STAT", "TGF", "WNT", "PI3K", "AKT", "MTOR", "RNA", "DNA", "ECM", "GPCR",
             "IL", "NK", "VEGF", "EGFR", "ERBB", "PPAR", "AMPK", "CAMP", "CGMP", "FOXO", "P53", "TCR", "BCR", "MHC",
             "HPV", "HIV", "ABC", "SNARE", "NOD", "TOLL", "HIF-1", "JAK-STAT", "PD-L1", "PD-1", "GABA", "AGE-RAGE",
             "ATP", "NADPH", "CGAS-STING", "RIG-I-LIKE", "FC", "EPSILON", "GAMMA", "HTLV-I", "EBV", "HCMV", "KSHV",
             "SARS-COV-2", "COVID-19", "II", "III", "IV", "I", "T", "B", "G2/M", "E2F", "MYC", "KRAS", "UV", "UPR"}
_STOP = {"of", "via", "by", "to", "in", "and", "the", "for", "on", "or", "with", "from", "at", "a", "an"}


def disp(name: str) -> str:
    """Sentence case for display. Title Case / ALL CAPS names (Enrichr, KEGG) are lowered word by word,
    keeping acronyms and mixed-case tokens (MAPK, IL-2/STAT5, ATPase, V-type, p53)."""
    s = re.sub(r"\s+", " ", str(name).strip())
    words = s.split(" ")
    content = [w for w in words if w.lower() not in _STOP] or words
    all_caps = sum(w.isupper() for w in content) / len(content) > 0.7
    capitalised = sum(w[:1].isupper() for w in content) / len(content) >= 0.6

    def low(w: str) -> str:
        if all_caps:
            return w if w in _ACRONYMS or re.search(r"\d", w) else w.lower()
        parts = re.split(r"([-/])", w)
        simple = all(p in "-/" or p == "" or (p[:1].isupper() and p[1:].islower()) or p.islower() for p in parts)
        return w.lower() if simple and not (len(w) == 1 and w.isupper()) else w

    if all_caps or capitalised:
        s = " ".join(low(w) for w in words)
        s = re.sub(r"\bdn\b", "down", s)
    return s[:1].upper() + s[1:]


def _wrap(s: str, width: int = 26, lines: int = 2) -> str:
    parts = textwrap.wrap(str(s), width, break_long_words=False)
    if len(parts) > lines:
        parts = parts[:lines]
        parts[-1] = parts[-1][: width - 1].rstrip() + "…"
    return "\n".join(parts)


def _short(s: str, n: int = 55) -> str:
    s = str(s)
    return s if len(s) <= n else s[: n - 1].rstrip() + "…"


def _term_labels(names, dbs, fmt) -> list[str]:
    """Display names for a list of terms; a name shown more than once gets its database appended."""
    shown = [disp(n) for n in names]
    dup = {n for n in shown if shown.count(n) > 1}
    return [fmt(f"{n} ({DB_LABEL.get(d, d)})" if n in dup else n) for n, d in zip(shown, dbs)]


def _mlog(p) -> np.ndarray:
    return np.minimum(MLOG_CAP, -np.log10(np.clip(np.asarray(p, dtype=float), 1e-300, 1.0)))


def _save(fig, path: Path) -> str:
    """Write <path>.pdf (vector, publication) and <path>.png (report preview); return the PNG path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"))            # page = figsize exactly (<= 7.2 in wide)
    fig.savefig(path.with_suffix(".png"), dpi=220)
    plt.close(fig)
    return str(path.with_suffix(".png"))


def _box(ax):
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.6)


# ───────────────────────────── layout: plot + legend column ─────────────────

def _figure(w: float, h: float, legend_w: float = 1.35):
    """Figure with the plot on the left and a legend column (axis-off) on the right."""
    fig = plt.figure(figsize=(min(w, FULL_W), max(h, 2.6)), layout="constrained")
    fig.get_layout_engine().set(w_pad=0.03, h_pad=0.03, wspace=0.02)
    gs = fig.add_gridspec(1, 2, width_ratios=[min(w, FULL_W) - legend_w, legend_w])
    ax = fig.add_subplot(gs[0])
    lax = fig.add_subplot(gs[1])
    lax.axis("off")
    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    return fig, ax, lax


class LegendColumn:
    """Stack size legends, colour bars and notes top-to-bottom in the legend axis."""

    def __init__(self, fig, lax, x0: float = 0.12):
        self.fig, self.lax, self.y, self.x0 = fig, lax, 0.99, x0
        fig.canvas.draw()
        self.h_in = lax.get_window_extent().height / fig.dpi   # legend axis height in inches

    def _dy(self, inches: float) -> float:
        return inches / max(self.h_in, 1e-6)

    def title(self, text: str):
        self.lax.text(self.x0, self.y, text, fontsize=9, fontweight="bold", va="top", ha="left")
        self.y -= self._dy(0.2 + 0.13 * text.count("\n"))

    def sizes(self, title: str, values, to_size, marker="o", color="white", fmt="{:g}"):
        self.title(title)
        for v in values:
            s = to_size(v)
            r_in = np.sqrt(s) / 72 / 2
            step = max(0.2, 2 * r_in + 0.07)
            yc = self.y - self._dy(step / 2)
            self.lax.scatter([self.x0 + 0.12], [yc], s=s, marker=marker, facecolor=color, edgecolor=INK,
                             linewidths=0.6, clip_on=False)
            self.lax.text(self.x0 + 0.34, yc, fmt.format(v), va="center", fontsize=8.5)
            self.y -= self._dy(step)
        self.y -= self._dy(0.12)

    def colorbar(self, title: str, norm, cmap, height_in: float = 0.9):
        if title:
            self.title(title)
        h = self._dy(height_in)
        cax = self.lax.inset_axes([self.x0, self.y - h, 0.12, h])
        cb = self.fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
        cb.outline.set_linewidth(0.5)
        cb.ax.tick_params(labelsize=8.5, width=0.5, length=2)
        self.y -= h + self._dy(0.22)

    def markers(self, title: str, items):
        """items: list of (label, dict(marker kwargs))."""
        self.title(title)
        for label, kw in items:
            yc = self.y - self._dy(0.1)
            self.lax.scatter([self.x0 + 0.12], [yc], clip_on=False, **kw)
            self.lax.text(self.x0 + 0.34, yc, label, va="center", fontsize=8.5)
            self.y -= self._dy(0.2)
        self.y -= self._dy(0.1)

    def note(self, text: str):
        self.lax.text(self.x0, self.y, text, fontsize=8.5, va="top", ha="left", color=INK2, linespacing=1.4)
        self.y -= self._dy(0.16 * (text.count("\n") + 1) + 0.1)


# ───────────────────────────── per-tool figures ─────────────────────────────

def tool_dotplot(rec: pd.DataFrame, title: str, path: Path, top: int = 15) -> str | None:
    """clusterProfiler-style dot plot: x = gene ratio, size = hit genes, colour = -log10 adj. P."""
    d = rec.sort_values("padj").drop_duplicates("term_key").head(top).copy()
    if d.empty:
        return None
    qs = pd.to_numeric(d["query_size"], errors="coerce")
    d["ratio"] = np.where(qs > 0, d["overlap"] / qs, np.nan)
    use_ratio = bool(d["ratio"].notna().all())
    d["x"] = d["ratio"] if use_ratio else _mlog(d["padj"])
    d["c"] = _mlog(d["padj"])
    d = d.sort_values("x")
    fig, ax, lax = _figure(FULL_W, 0.24 * len(d) + 0.9)
    to_size = lambda k: 10 + 7 * min(float(k), 40)  # noqa: E731
    norm = Normalize(0, MLOG_CAP)
    ax.scatter(d["x"], np.arange(len(d)), s=[to_size(k) for k in d["overlap"].fillna(1)], c=d["c"], cmap=REDS,
               norm=norm, edgecolors=INK, linewidths=0.5, zorder=3)
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels([_wrap(disp(n), 48, 2) for n in d["term_name"]], linespacing=0.95)
    ax.set_ylim(-0.7, len(d) - 0.3)
    xr = d["x"].max() - d["x"].min()
    ax.set_xlim(d["x"].min() - 0.12 * xr - 1e-3, d["x"].max() + 0.12 * xr + 1e-3)
    ax.set_xlabel("Gene ratio (hits / query genes)" if use_ratio else MLOG_LABEL)
    ax.grid(axis="both", color=GRID, lw=0.4, zorder=0)
    _box(ax)
    ax.set_title(title)
    L = LegendColumn(fig, lax)
    kmax = int(max(3, d["overlap"].max()))
    L.sizes("Hit genes", sorted({max(1, round(kmax / 4)), max(2, round(kmax / 2)), kmax}), to_size)
    L.colorbar(MLOG_LABEL, norm, REDS)
    if d["padj"].min() >= 0.05:
        L.note("No term reaches\nadj. P < 0.05")
    return _save(fig, path)


def cluster_terms(d: pd.DataFrame, cut: float = 0.3) -> np.ndarray:
    """Kappa / average-linkage clustering of terms by hit genes (for per-tool networks)."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    genes = sorted({g for h in d["hit_genes"] for g in str(h).split(";") if g})
    if len(d) < 2 or not genes:
        return np.ones(len(d), dtype=int)
    gpos = {g: i for i, g in enumerate(genes)}
    X = np.zeros((len(d), len(genes)), dtype=np.int8)
    for i, h in enumerate(d["hit_genes"]):
        for g in str(h).split(";"):
            if g:
                X[i, gpos[g]] = 1
    D = np.clip(1 - kappa_matrix(X), 0, 2)
    np.fill_diagonal(D, 0)
    return fcluster(linkage(squareform(D, checks=False), "average"), t=1 - cut, criterion="distance")


def tool_network(rec: pd.DataFrame, title: str, path: Path, top: int = 40, max_term_size: int = 1000) -> str | None:
    """Enrichment map of one tool's top terms (very generic terms removed), clustered by shared hit genes."""
    r = rec[pd.to_numeric(rec["term_size"], errors="coerce").fillna(0) <= max_term_size]
    d = r[r["padj"] < 0.05].sort_values("padj").drop_duplicates("term_key").head(top)
    if len(d) < 3:
        d = r.sort_values("padj").drop_duplicates("term_key").head(15)
    d = d[d["hit_genes"].astype(str).str.len() > 0].copy()
    if len(d) < 2:
        return None
    d["cl"] = cluster_terms(d)
    terms = pd.DataFrame({"term_key": d["term_key"].values, "term_name": d["term_name"].values,
                          "padj": d["padj"].values,
                          "size": pd.to_numeric(d["term_size"], errors="coerce").fillna(d["overlap"]).values,
                          "hit_genes": d["hit_genes"].values, "tier": 2,
                          "direction": d["direction"].values if "direction" in d else 0})
    groups = []
    for _, m in terms.assign(cl=d["cl"].values).groupby("cl"):
        m = m.sort_values("padj")
        groups.append({"label": disp(m.iloc[0]["term_name"]), "keys": m["term_key"].tolist(),
                       "best": float(m["padj"].min())})
    groups = sorted(groups, key=lambda g: (-len(g["keys"]), g["best"]))
    png, _ = enrichment_map(terms, groups[:12], title, path, show_tier=False, n_hidden=max(0, len(groups) - 12))
    return png


# ───────────────────────────── enrichment map (network) ─────────────────────

def _pack(boxes: list[tuple[float, float]], row_w: float, gap: float = 0.3):
    """Shelf-pack boxes (width, height incl. label) into rows of width row_w; returns box centres."""
    centres, x, y, row_h = [], 0.0, 0.0, 0.0
    for w, h in boxes:
        if x + w > row_w + 1e-9 and x > 0:
            y -= row_h + gap
            x, row_h = 0.0, 0.0
        centres.append((x + w / 2, y - h / 2))
        x += w + gap
        row_h = max(row_h, h)
    return centres


def _label_size(text: str, fs: float = 7.5) -> tuple[float, float]:
    """Approximate (width, height) in inches of multi-line Arial text."""
    lines = text.split("\n")
    return max(len(l) for l in lines) * fs * 0.5 / 72, len(lines) * fs * 1.2 / 72


def enrichment_map(terms: pd.DataFrame, groups: list[dict], title: str, path: Path,
                   show_tier: bool = True, n_hidden: int = 0) -> tuple[str | None, dict]:
    """
    EnrichmentMap / Metascape-style network, 7.2 in wide (1 data unit = 1 inch on the page).
      terms : term_key, term_name, padj, size, hit_genes, tier, direction (+ optional score, database, support)
      groups: [{"label", "keys"}], one per theme: outlined cluster with its label above
    Node size = gene-set size; fill = -log10 adj. P (red), or signed by direction (blue/red);
    outline = consensus tier; edges = shared hit genes (kappa >= 0.3) within a theme, and
    kappa >= 0.5 between themes (faint). Legend strip below the network.
    """
    terms = terms.drop_duplicates("term_key").set_index("term_key")
    groups = [dict(g, keys=[k for k in g["keys"] if k in terms.index]) for g in groups]
    groups = [g for g in groups if g["keys"]]
    if not groups:
        return None, {}
    # multi-term themes first (big to small), single-term themes packed together at the end
    groups = sorted(groups, key=lambda g: (len(g["keys"]) == 1, 0))
    keys = [k for g in groups for k in g["keys"]]
    genes = sorted({x for h in terms.loc[keys, "hit_genes"] for x in str(h).split(";") if x})
    gpos = {g: i for i, g in enumerate(genes)}
    X = np.zeros((len(keys), max(1, len(genes))), dtype=np.int8)
    for i, h in enumerate(terms.loc[keys, "hit_genes"]):
        for x in str(h).split(";"):
            if x:
                X[i, gpos[x]] = 1
    K = kappa_matrix(X)
    kidx = {k: i for i, k in enumerate(keys)}

    boxes, locals_, rads = [], [], []
    for g in groups:
        n = len(g["keys"])
        g["_text"] = _wrap(g["label"], 24, 3)
        lw_, lh_ = _label_size(g["_text"])
        if n == 1:
            lh_ = _label_size("x\nx\nx")[1]          # all single-term boxes share one height
        if n == 1:
            loc, r = {g["keys"][0]: np.zeros(2)}, 0.12
        else:
            G = nx.Graph()
            G.add_nodes_from(g["keys"])
            for a_ in range(n):
                for b_ in range(a_ + 1, n):
                    w = K[kidx[g["keys"][a_]], kidx[g["keys"][b_]]]
                    if w >= 0.3:
                        G.add_edge(g["keys"][a_], g["keys"][b_], weight=w)
            loc = nx.spring_layout(G, weight="weight", seed=11, k=1.4 / np.sqrt(n), iterations=200)
            arr = np.array(list(loc.values()))
            arr = arr - arr.mean(axis=0)
            arr = arr / (np.abs(arr).max() + 1e-9)
            loc = dict(zip(loc.keys(), arr))
            r = 0.16 + 0.135 * np.sqrt(n)
        rads.append(r)
        locals_.append(loc)
        boxes.append((max(2 * r + 0.25, lw_ + 0.1), 2 * r + 0.2 + lh_ + 0.06))
    plot_w = FULL_W - 0.15
    centres = _pack(boxes, plot_w)
    pos = {}
    for g, (cx, cy), (bw, bh), r, loc in zip(groups, centres, boxes, rads, locals_):
        node_c = np.array([cx, cy - bh / 2 + r + 0.1])            # cluster sits at the bottom of its box
        for k, v in loc.items():
            pos[k] = node_c + v * max(0.0, r - 0.12)
        g["_label_x"], g["_label_y"] = cx, node_c[1] + r + 0.08
        g["_node_c"], g["_r"] = node_c, r

    direction = terms.loc[keys, "direction"].fillna(0).astype(float).values if "direction" in terms \
        else np.zeros(len(keys))
    mlog = _mlog(terms.loc[keys, "padj"].values)
    if np.abs(direction).sum() > 0:
        val = mlog * np.sign(direction)
        norm, cmap, clabel = TwoSlopeNorm(0, -MLOG_CAP, MLOG_CAP), RDBU, "Signed " + MLOG_LABEL
    else:
        val, norm, cmap, clabel = mlog, Normalize(0, MLOG_CAP), REDS, MLOG_LABEL

    x_lo = min(c[0] - b[0] / 2 for c, b in zip(centres, boxes)) - 0.05
    y_lo = min(c[1] - b[1] / 2 for c, b in zip(centres, boxes)) - 0.05
    y_hi = max(c[1] + b[1] / 2 for c, b in zip(centres, boxes)) + 0.05
    plot_h = y_hi - y_lo
    legend_h, title_h = 0.95, 0.3
    fig = plt.figure(figsize=(FULL_W, plot_h + legend_h + title_h))
    ax = fig.add_axes([0.075 / FULL_W, (legend_h) / (plot_h + legend_h + title_h),
                       (FULL_W - 0.15) / FULL_W, plot_h / (plot_h + legend_h + title_h)])
    group_of = {k: gi for gi, g in enumerate(groups) for k in g["keys"]}

    segs_in, w_in = [], []
    between: dict[tuple, list] = {}
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            gi, gj = group_of[keys[i]], group_of[keys[j]]
            if gi == gj and K[i, j] >= 0.3:
                segs_in.append([pos[keys[i]], pos[keys[j]]])
                w_in.append(0.4 + 1.8 * K[i, j])
            elif gi != gj and K[i, j] >= 0.5:
                between.setdefault(tuple(sorted((gi, gj))), []).append(K[i, j])
    # one bundled edge per pair of themes (width = number of links), drawn beneath the theme discs
    segs_out, w_out = [], []
    def _crosses(p0, p1, skip):
        d = p1 - p0
        for gk, g in enumerate(groups):
            if gk in skip:
                continue
            t = np.clip(np.dot(g["_node_c"] - p0, d) / max(np.dot(d, d), 1e-9), 0, 1)
            if np.linalg.norm(p0 + t * d - g["_node_c"]) < g["_r"] + 0.12:
                return True
        return False

    for (gi, gj), ks in between.items():
        if _crosses(groups[gi]["_node_c"], groups[gj]["_node_c"], {gi, gj}):
            continue
        segs_out.append([groups[gi]["_node_c"], groups[gj]["_node_c"]])
        w_out.append(0.5 + 0.5 * np.sqrt(len(ks)))
    if segs_out:
        ax.add_collection(LineCollection(segs_out, colors=EDGE_BLUE, linewidths=w_out, alpha=0.45, zorder=0.1))
    if segs_in:
        ax.add_collection(LineCollection(segs_in, colors=EDGE_BLUE, linewidths=w_in, alpha=0.95, zorder=1))

    for g in groups:
        if len(g["keys"]) > 1:
            ax.add_patch(Ellipse(g["_node_c"], 2 * g["_r"] + 0.12, 2 * g["_r"] + 0.12, facecolor="white", lw=0.6,
                                 ec=INK, zorder=0.5))
        ax.text(g["_label_x"], g["_label_y"], g["_text"], ha="center", va="bottom", fontsize=8.5, color=INK,
                zorder=4, linespacing=1.0, bbox=dict(fc="white", ec="none", pad=0.4))

    xs = np.array([pos[k][0] for k in keys])
    ys = np.array([pos[k][1] for k in keys])
    sizes = pd.to_numeric(terms.loc[keys, "size"], errors="coerce").fillna(20).clip(5, 500).values
    to_s = lambda n: 8 + 4.5 * np.sqrt(n)  # noqa: E731
    tiers = terms.loc[keys, "tier"].fillna(3).astype(int).values if show_tier else np.full(len(keys), 2)
    for t in (3, 2, 1):
        idx = np.where(tiers == t)[0]
        if not len(idx):
            continue
        ec, lw, ls = TIER_EDGE[t]
        ax.scatter(xs[idx], ys[idx], s=to_s(sizes[idx]), c=val[idx], cmap=cmap, norm=norm, edgecolors=ec,
                   linewidths=lw, linestyles=ls, zorder=3)
    ax.set_xlim(x_lo, x_lo + (FULL_W - 0.15))
    ax.set_ylim(y_lo, y_hi)
    ax.axis("off")
    fig.text(0.075 / FULL_W, 1 - 0.06 / (plot_h + legend_h + title_h), title, fontsize=10, fontweight="bold",
             va="top", ha="left")
    if n_hidden:
        fig.text(1 - 0.075 / FULL_W, 1 - 0.06 / (plot_h + legend_h + title_h),
                 f"{len(groups)} of {len(groups) + n_hidden} themes shown; full list in themes table",
                 fontsize=8.5, color=INK2, va="top", ha="right")

    # legend strip (four blocks) along the bottom
    H = plot_h + legend_h + title_h
    lg = fig.add_axes([0.075 / FULL_W, 0.02 / H, (FULL_W - 0.15) / FULL_W, (legend_h - 0.1) / H])
    lg.axis("off")
    lg.set_xlim(0, FULL_W - 0.15)
    lg.set_ylim(0, legend_h - 0.1)
    top = legend_h - 0.18
    lg.text(0, top, "Gene-set size", fontsize=9, fontweight="bold", va="center")
    for i, n in enumerate((10, 100, 500)):
        lg.scatter(0.15 + i * 0.45, 0.38, s=to_s(n), facecolor="white", edgecolor=INK, lw=0.6)
        lg.text(0.15 + i * 0.45, 0.08, str(n), ha="center", fontsize=8.5)
    lg.text(1.75, top, clabel, fontsize=9, fontweight="bold", va="center")
    cax = lg.inset_axes([1.75, 0.3, 1.3, 0.12], transform=lg.transData)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=8.5, width=0.5, length=2)
    x_t = 3.55
    if show_tier:
        lg.text(x_t, top, "Consensus tier", fontsize=9, fontweight="bold", va="center")
        for i, t in enumerate((1, 2, 3)):
            ec, lw, ls = TIER_EDGE[t]
            lg.scatter(x_t + 0.08, 0.5 - i * 0.2, s=30, facecolor="white", edgecolor=ec, lw=lw, linestyle=ls)
            lg.text(x_t + 0.22, 0.5 - i * 0.2, TIER_NAMES[t], fontsize=8.5, va="center")
        x_e = 5.2
    else:
        x_e = 3.55
    lg.text(x_e, top, "Shared hit genes", fontsize=9, fontweight="bold", va="center")
    for i, (kappa, lw) in enumerate(((0.3, 0.4 + 1.8 * 0.3), (1.0, 2.2))):
        lg.plot([x_e, x_e + 0.35], [0.45 - i * 0.25] * 2, color=EDGE_BLUE, lw=lw)
        lg.text(x_e + 0.45, 0.45 - i * 0.25, f"κ = {kappa:g}", va="center", fontsize=8.5)
    png = _save(fig, path)

    def _get(k, col, default=""):
        return terms.loc[k, col] if col in terms.columns else default

    cy = {"nodes": [{"data": {"id": k, "label": disp(terms.loc[k, "term_name"]), "tier": int(tiers[i]),
                              "score": float(_get(k, "score", 0) or 0), "db": str(_get(k, "database")),
                              "theme": groups[group_of[k]]["label"], "padj": float(terms.loc[k, "padj"]),
                              "genes": str(terms.loc[k, "hit_genes"]), "support": str(_get(k, "support"))},
                     "position": {"x": float(pos[k][0] * 110), "y": float(-pos[k][1] * 110)}}
                    for i, k in enumerate(keys)],
          "edges": [{"data": {"source": keys[i], "target": keys[j], "weight": round(float(K[i, j]), 3)}}
                    for i in range(len(keys)) for j in range(i + 1, len(keys)) if K[i, j] >= 0.3]}
    return png, cy


def theme_network(cons: pd.DataFrame, themes: pd.DataFrame, title: str, path: Path,
                  n_themes: int = 12, per_theme: int = 8) -> tuple[str | None, dict]:
    if themes is None or themes.empty:
        return None, {}
    c = cons.set_index("term_key")
    groups = []
    for th in themes.head(n_themes).itertuples():
        ks = [k for k in th.member_keys.split("|") if k in c.index][:per_theme]
        groups.append({"label": f"T{th.theme_rank} {disp(th.theme_name)}", "keys": ks})
    keys = [k for g in groups for k in g["keys"]]
    sub = c.loc[keys].reset_index()
    terms = pd.DataFrame({"term_key": sub["term_key"], "term_name": sub["term_name"], "padj": sub["best_padj"],
                          "size": sub["term_size"], "hit_genes": sub["hit_genes"], "tier": sub["tier"],
                          "direction": sub["direction"], "score": sub["consensus_score"],
                          "database": sub["database"], "support": sub["tool_support"]})
    return enrichment_map(terms, groups, title, path, show_tier=True, n_hidden=max(0, len(themes) - n_themes))


# ───────────────────────────── consensus figures ────────────────────────────

def consensus_bar(cons: pd.DataFrame, themes: pd.DataFrame, title: str, path: Path,
                  n_themes: int = 10, per_theme: int = 4) -> str | None:
    """Top terms grouped by theme: one colour per theme with a side strip + T# label; bar label = hit genes."""
    if cons.empty or themes is None or themes.empty:
        return None
    c = cons[cons["tier"] > 0].set_index("term_key")
    rows = []
    for gi, th in enumerate(themes.head(n_themes).itertuples()):
        ks = [k for k in th.member_keys.split("|") if k in c.index]
        sub = c.loc[ks].sort_values("consensus_score", ascending=False).head(per_theme)
        for _, r in sub.iterrows():
            rows.append({**r.to_dict(), "g": gi, "theme_rank_": th.theme_rank})
    if not rows:
        return None
    d = pd.DataFrame(rows)
    y = np.arange(len(d))[::-1]
    fig, ax, lax = _figure(FULL_W, 0.2 * len(d) + 0.75, legend_w=1.25)
    alpha = d["tier"].map({1: 1.0, 2: 0.65, 3: 0.45}).values
    for yi, (_, r), a in zip(y, d.iterrows(), alpha):
        col = CAT[r["g"] % len(CAT)]
        ax.barh(yi, r["consensus_score"], height=0.74, color=col, alpha=a, edgecolor=col, linewidth=0.5)
        ax.text(r["consensus_score"] + 0.8, yi, f"{int(r['n_hit_genes'])}", va="center", fontsize=8.5)
    ax.set_yticks(y)
    ax.set_yticklabels(_term_labels(d["term_name"], d["database"], lambda n: _short(n, 60)))
    ax.set_ylim(-0.6, len(d) - 0.4)
    xmax = float(np.ceil((d["consensus_score"].max() + 9) / 10) * 10)
    ax.set_xlim(0, xmax)
    for gi in d["g"].unique():
        ys = y[d["g"].values == gi]
        col = CAT[gi % len(CAT)]
        ax.plot([xmax * 1.015, xmax * 1.015], [ys.min() - 0.38, ys.max() + 0.38], color=col, lw=3,
                solid_capstyle="butt", clip_on=False)
        ax.text(xmax * 1.035, ys.mean(), f"T{d.loc[d['g'] == gi, 'theme_rank_'].iat[0]}", va="center",
                ha="left", fontsize=9, fontweight="bold")
        if ys.min() > 0:
            ax.axhline(ys.min() - 0.5, color=GRID, lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("Consensus score")
    ax.set_title(title)
    L = LegendColumn(fig, lax, x0=0.28)
    L.markers("Tier (opacity)", [(TIER_NAMES[t], dict(s=45, marker="s", color="#7F7F7F", alpha=a, edgecolor="none"))
                                 for t, a in ((1, 1.0), (2, 0.65), (3, 0.45))])
    L.note("Bar label:\nnumber of hit genes\nT#: theme")
    return _save(fig, path)


def evidence_bubble(cons: pd.DataFrame, themes: pd.DataFrame, tools: list[str], tool_dbs: dict, title: str,
                    path: Path, n_themes: int = 8, per_theme: int = 4) -> str | None:
    """
    Term x tool evidence matrix: rows grouped by theme; bubble size and colour = -log10 adj. P (one
    shared 0-8 scale); outlined bubble = significant; grey dot = tested, not significant;
    '–' = database not covered. Tools that cover none of the shown terms are dropped.
    """
    cols = [t for t in tools if f"padj__{t}" in cons.columns]
    if cons.empty or not cols or themes is None or themes.empty:
        return None
    c = cons.set_index("term_key")
    rows = []
    for gi, th in enumerate(themes.head(n_themes).itertuples()):
        ks = [k for k in th.member_keys.split("|") if k in c.index and c.loc[k, "tier"] > 0]
        for k in sorted(ks, key=lambda k: -c.loc[k, "consensus_score"])[:per_theme]:
            rows.append((gi, th.theme_rank, k))
    if not rows:
        return None
    cols = [t for t in cols if any(c.loc[k, "database"] in tool_dbs.get(t, []) for _, _, k in rows)]
    n, m = len(rows), len(cols)
    fig, ax, lax = _figure(FULL_W, 0.27 * n + 1.2, legend_w=1.3)
    to_s = lambda v: 6 + 14 * v  # noqa: E731
    norm = Normalize(0, MLOG_CAP)
    any_na = False
    for i, (gi, _, k) in enumerate(rows):
        yy = n - 1 - i
        for j, t in enumerate(cols):
            p = c.loc[k, f"padj__{t}"]
            covered = c.loc[k, "database"] in tool_dbs.get(t, [])
            if pd.isna(p) and not covered:
                any_na = True
                ax.text(j, yy, "–", ha="center", va="center", fontsize=9, color=NA_GREY)
            elif pd.isna(p) or p >= 0.05:
                ax.scatter(j, yy, s=6, c=NA_GREY, edgecolors="none", zorder=2)
            else:
                v = float(_mlog(p))
                ax.scatter(j, yy, s=to_s(v), c=[v], cmap=REDS, norm=norm, edgecolors=INK, linewidths=0.5, zorder=3)
    ax.set_xticks(range(m))
    ax.set_xticklabels([tool_label(t) for t in cols], rotation=40, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(_term_labels([c.loc[k, "term_name"] for _, _, k in rows],
                                    [c.loc[k, "database"] for _, _, k in rows], lambda s: _wrap(s, 46, 2))[::-1],
                       linespacing=0.95)
    ax.set_xlim(-0.6, m - 0.4)
    ax.set_ylim(-0.6, n - 0.4)
    prev = None
    for i, (gi, rank, _) in enumerate(rows):
        if gi != prev:
            yy = n - 1 - i
            if i:
                ax.axhline(yy + 0.5, color=GRID, lw=0.6)
            ax.annotate(f"T{rank}", xy=(1, yy), xycoords=("axes fraction", "data"), xytext=(4, 0),
                        textcoords="offset points", fontsize=8.5, fontweight="bold", ha="left", va="center")
            prev = gi
    _box(ax)
    ax.set_title(title)
    L = LegendColumn(fig, lax, x0=0.2)
    L.sizes(MLOG_LABEL + "\n(size and colour)", [2, 4, 8], to_s, color="#F7A48F")
    L.colorbar("", norm, REDS, height_in=0.7)
    key = [("P$_{adj}$ < 0.05", dict(s=30, facecolor="#F7A48F", edgecolor=INK, linewidths=0.5)),
           ("tested, n.s.", dict(s=6, color=NA_GREY))]
    if any_na:
        key.append(("not covered", dict(s=30, marker="_", color=NA_GREY)))
    L.markers("Key", key)
    return _save(fig, path)


def theme_tool_squares(records: pd.DataFrame, themes: pd.DataFrame, tools: list[str], title: str,
                       path: Path, n_themes: int = 16) -> str | None:
    """Themes (rows) x tools (columns); square size = number of the theme's terms significant in that tool."""
    if themes is None or themes.empty:
        return None
    th = themes.head(n_themes)
    sig = records[records["padj"] < 0.05]
    tools = [t for t in tools if t in set(sig["tool"])]
    if not tools:
        return None
    M = np.zeros((len(th), len(tools)))
    for i, r in enumerate(th.itertuples()):
        ks = set(r.member_keys.split("|"))
        for j, t in enumerate(tools):
            M[i, j] = sig[(sig["tool"] == t) & sig["term_key"].isin(ks)]["term_key"].nunique()
    fig, ax, lax = _figure(min(FULL_W, 3.3 + 0.45 * len(tools) + 1.3), 0.3 * len(th) + 1.2, legend_w=1.3)
    to_s = lambda v: 0 if v == 0 else 14 + 26 * np.sqrt(v)  # noqa: E731
    for i in range(len(th)):
        for j in range(len(tools)):
            if M[i, j]:
                ax.scatter(j, len(th) - 1 - i, marker="s", s=to_s(M[i, j]), c="#3C5488", edgecolors="none")
    ax.set_yticks(range(len(th)))
    ax.set_yticklabels([_wrap(f"T{r.theme_rank}  {disp(r.theme_name)}", 38, 2) for r in th.itertuples()][::-1],
                       linespacing=0.95)
    ax.set_xticks(range(len(tools)))
    ax.set_xticklabels([tool_label(t) for t in tools], rotation=40, ha="right")
    ax.set_xlim(-0.6, len(tools) - 0.4)
    ax.set_ylim(-0.6, len(th) - 0.4)
    ax.grid(color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    _box(ax)
    ax.set_title(title)
    vmax = int(max(1, M.max()))
    L = LegendColumn(fig, lax, x0=0.2)
    L.sizes("Significant\nterms per theme", sorted({1, max(1, vmax // 2), vmax}), to_s, marker="s", color="#3C5488")
    return _save(fig, path)


def tool_overlap(records: pd.DataFrame, tools: list[str], title: str, path: Path) -> str | None:
    sig = {t: set(records.loc[(records["tool"] == t) & (records["padj"] < 0.05), "term_key"]) for t in tools}
    tools = [t for t in tools if sig.get(t)]
    if len(tools) < 2:
        return None
    J = np.array([[len(sig[a] & sig[b]) / max(len(sig[a] | sig[b]), 1) for b in tools] for a in tools])
    n = len(tools)
    fig, ax = plt.subplots(figsize=(min(FULL_W, 2.3 + 0.4 * n + 0.7), 0.4 * n + 1.1), layout="constrained")
    im = ax.imshow(J, cmap=NAVY, vmin=0, vmax=1)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{J[i, j]:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if J[i, j] > 0.55 else INK)
    ax.set_xticks(range(n))
    ax.set_xticklabels([tool_label(t) for t in tools], rotation=40, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels([f"{tool_label(t)} (n = {len(sig[t])})" for t in tools])
    ax.set_xticks(np.arange(-0.5, n), minor=True)
    ax.set_yticks(np.arange(-0.5, n), minor=True)
    ax.grid(which="minor", color="white", lw=1.0)
    ax.tick_params(which="minor", length=0)
    ax.set_title(title)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Jaccard index of significant terms")
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=8.5, width=0.5, length=2)
    return _save(fig, path)


def _gene_theme_fraction(genes: list[str], th: pd.DataFrame) -> np.ndarray:
    F = np.zeros((len(genes), len(th)))
    for j, r in enumerate(th.itertuples()):
        freq = dict(x.rsplit(":", 1) for x in str(r.gene_frequency).split(";") if ":" in x)
        for i, g in enumerate(genes):
            F[i, j] = int(freq.get(g, 0)) / max(r.n_terms, 1)
    return F


def _theme_ticklabels(th: pd.DataFrame) -> list[str]:
    return [_wrap(f"T{r.theme_rank} {disp(r.theme_name)}", 26, 3) for r in th.itertuples()]


def key_gene_bubble(key_genes: pd.DataFrame, themes: pd.DataFrame, title: str, path: Path,
                    n_genes: int = 25, n_themes: int = 8) -> str | None:
    """
    Key genes x themes: size = fraction of the theme's terms containing the gene; colour = the gene's
    input score (diverging if signed) or gene score (navy). Marks: '*' driver (only when it
    discriminates), '^' gene of interest.
    """
    if key_genes is None or key_genes.empty or themes is None or themes.empty:
        return None
    kg = key_genes.head(n_genes)
    th = themes.head(n_themes)
    genes = kg["gene"].tolist()
    F = _gene_theme_fraction(genes, th)
    if kg["input_score"].notna().any():
        vals = kg["input_score"].fillna(0).values.astype(float)
        vmax = max(1e-9, float(np.abs(vals).max()))
        norm, cmap = (TwoSlopeNorm(0, -vmax, vmax), RDBU) if vals.min() < 0 else (Normalize(0, vmax), REDS)
        clab = "Input score"
    else:
        vals = kg["gene_score"].values.astype(float)
        norm, cmap, clab = Normalize(0, float(vals.max())), NAVY, "Gene score"
    fig, ax, lax = _figure(min(FULL_W, 1.4 + 0.42 * len(th) + 1.3), 0.19 * len(genes) + 2.0, legend_w=1.3)
    to_s = lambda f: 8 + 90 * f  # noqa: E731
    for i in range(len(genes)):
        for j in range(len(th)):
            if F[i, j] > 0:
                ax.scatter(j, len(genes) - 1 - i, s=to_s(F[i, j]), c=[vals[i]], cmap=cmap, norm=norm,
                           edgecolors=INK, linewidths=0.5, zorder=3)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes[::-1], fontstyle="italic")
    ax.set_xticks(range(len(th)))
    ax.set_xticklabels(_theme_ticklabels(th), rotation=90, ha="center", va="top", fontsize=8.5, linespacing=0.95)
    ax.set_xlim(-0.6, len(th) - 0.4)
    ax.set_ylim(-0.6, len(genes) - 0.4)
    drivers = kg["driver_of"].astype(str).str.len() > 0
    mark_drivers = 0 < drivers.sum() < len(kg)
    marks = []
    for i, r in enumerate(kg.itertuples()):
        mk = ("*" if mark_drivers and drivers.iat[i] else "") + ("^" if r.gene_of_interest else "")
        if mk:
            ax.text(len(th) - 0.35, len(genes) - 1 - i, mk, va="center", fontsize=10, clip_on=False)
            marks.append(mk)
    _box(ax)
    ax.set_title(title)
    L = LegendColumn(fig, lax, x0=0.22)
    L.sizes("Fraction of\ntheme terms", [0.25, 0.5, 1.0], to_s)
    L.colorbar(clab, norm, cmap, height_in=0.7)
    notes = (["* driver gene"] if mark_drivers else []) + (["^ gene of interest"] if any("^" in m for m in marks) else [])
    if notes:
        L.note("\n".join(notes))
    return _save(fig, path)


def gene_theme_heatmap(key_genes: pd.DataFrame, themes: pd.DataFrame, title: str, path: Path,
                       n_genes: int = 30, n_themes: int = 10) -> str | None:
    """Membership heatmap: themes in T# order, genes ordered by hierarchical clustering."""
    from scipy.cluster.hierarchy import leaves_list, linkage
    if key_genes is None or key_genes.empty or themes is None or themes.empty:
        return None
    th = themes.head(n_themes)
    genes = key_genes.head(n_genes)["gene"].tolist()
    M = _gene_theme_fraction(genes, th)
    ro = leaves_list(linkage(M, "average")) if len(genes) > 2 else np.arange(len(genes))
    M = M[ro]
    fig, ax, lax = _figure(min(FULL_W, 1.3 + 0.36 * len(th) + 1.3), 0.17 * len(genes) + 2.0, legend_w=1.3)
    ax.imshow(M, aspect="auto", cmap=NAVY, vmin=0, vmax=1)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([genes[i] for i in ro], fontstyle="italic")
    ax.set_xticks(range(len(th)))
    ax.set_xticklabels(_theme_ticklabels(th), rotation=90, ha="center", va="top", fontsize=8.5, linespacing=0.95)
    ax.set_xticks(np.arange(-0.5, len(th)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes)), minor=True)
    ax.grid(which="minor", color="white", lw=0.8)
    ax.tick_params(which="minor", length=0)
    _box(ax)
    ax.set_title(title)
    L = LegendColumn(fig, lax, x0=0.2)
    L.colorbar("Fraction of theme\nterms with gene", Normalize(0, 1), NAVY, height_in=0.9)
    return _save(fig, path)


def metascape_clusters(cons: pd.DataFrame, themes: pd.DataFrame, clusters: list[dict], title: str, path: Path,
                       n_clusters: int = 20, n_themes: int = 12) -> str | None:
    """
    Metascape's own term clusters (rows, by best P) against the consensus themes (columns): Jaccard index
    of their hit genes. The right-hand column gives the consensus tier of each cluster's summary term.
    """
    if not clusters or cons.empty or themes is None or themes.empty:
        return None
    c = cons.drop_duplicates("term_key").set_index("term_key")
    ms = pd.DataFrame(clusters)
    hits = dict(zip(ms["term_key"], ms["hit_genes"]))

    def genes(keys, source):
        return {g for k in keys for g in str(source.get(k, "")).split(";") if g}

    rows = []
    for gid, m in ms.groupby("GROUP_ID"):
        m = m.sort_values("LogP")
        rows.append((float(m["LogP"].iat[0]), int(gid), str(m["Description"].iat[0]), m["term_key"].iat[0],
                     genes(m["term_key"], hits)))
    rows = sorted(rows)[:n_clusters]
    th = themes.head(n_themes)
    cons_hits = c["hit_genes"].to_dict()
    th_genes = [genes(str(r.member_keys).split("|"), cons_hits) for r in th.itertuples()]
    J = np.array([[len(g & t) / max(len(g | t), 1) for t in th_genes] for *_, g in rows])
    fig, ax, lax = _figure(FULL_W, 0.2 * len(rows) + 2.1, legend_w=1.3)
    ax.imshow(J, aspect="auto", cmap=NAVY, vmin=0, vmax=1)
    for i in range(J.shape[0]):
        for j in range(J.shape[1]):
            if J[i, j] >= 0.2:
                ax.text(j, i, f"{J[i, j]:.1f}", ha="center", va="center", fontsize=7.5,
                        color="white" if J[i, j] > 0.55 else INK)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([_short(f"M{i + 1}  {disp(r[2])}", 46) for i, r in enumerate(rows)])
    ax.set_xticks(range(len(th)))
    ax.set_xticklabels([_short(f"T{r.theme_rank} {disp(r.theme_name)}", 36) for r in th.itertuples()],
                       rotation=90, ha="center", va="top", fontsize=8.5)
    ax.set_xticks(np.arange(-0.5, len(th)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows)), minor=True)
    ax.grid(which="minor", color="white", lw=0.8)
    ax.tick_params(which="minor", length=0)
    for i, r in enumerate(rows):
        tier = int(c["tier"].get(r[3], 0) or 0)
        ax.annotate(TIER_NAMES.get(tier, "not significant").lower(), xy=(1, i), xycoords=("axes fraction", "data"),
                    xytext=(4, 0), textcoords="offset points", fontsize=8, va="center", ha="left",
                    color=INK if tier in (1, 2) else INK2, fontweight="bold" if tier == 1 else "normal")
    _box(ax)
    ax.set_title(title)
    L = LegendColumn(fig, lax, x0=0.62)
    L.colorbar("Shared hit genes\n(Jaccard index)", Normalize(0, 1), NAVY, height_in=0.9)
    L.note("M#: Metascape cluster,\nnamed by its best term\nRight of map: consensus\ntier of that term")
    return _save(fig, path)


def term_landscape(landscape: pd.DataFrame, cons: pd.DataFrame, db: str, title: str, path: Path,
                   n_label: int = 8) -> str | None:
    """Map of every term in one database (grey); significant consensus terms coloured by map cluster."""
    if landscape is None or landscape.empty:
        return None
    sig = cons[(cons["tier"] > 0) & (cons["database"] == db)].drop_duplicates("term_key").set_index("term_key")
    L_ = landscape.copy()
    L_["sig"] = L_["term_key"].isin(sig.index)
    if not L_["sig"].any():
        return None
    fg = np.where(L_["sig"].values)[0]
    # a map cluster keeps its colour + marker in every run: both depend only on the cluster id
    clusters = sorted(int(c) for c in L_["cluster"].unique())
    n_sig = L_.iloc[fg]["cluster"].value_counts()
    fig, ax, lax = _figure(FULL_W, max(5.2, 0.18 * len(clusters) + 1.3), legend_w=2.45)
    ax.set_box_aspect(1)                    # square map
    ax.set_anchor("N")
    bg = ~L_["sig"].values
    ax.scatter(L_.loc[bg, "x"], L_.loc[bg, "y"], s=3, c="#D9D9D9", edgecolors="none", rasterized=True)
    scores = np.array([sig.loc[L_["term_key"].iat[i], "consensus_score"] for i in fg])
    fg_cl = L_["cluster"].values[fg].astype(int)
    for m in sorted({cluster_style(c)[1] for c in fg_cl}):
        k = np.array([cluster_style(c)[1] == m for c in fg_cl])
        ax.scatter(L_["x"].values[fg][k], L_["y"].values[fg][k], s=14 + 0.9 * scores[k], marker=m,
                   c=[cluster_style(c)[0] for c in fg_cl[k]], edgecolors=INK, linewidths=0.6, zorder=3)
    ax.set_xticks([])
    ax.set_yticks([])
    _box(ax)
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    ax.set_title(title)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_bb = ax.get_window_extent(renderer)
    pts_disp = ax.transData.transform(np.c_[L_["x"].values[fg], L_["y"].values[fg]])
    placed, leaders = [], []
    cx = np.median(L_["x"])
    halo = [pe.withStroke(linewidth=2.2, foreground="white")]

    def _text_bbox(a):                      # text only: Annotation.get_window_extent also spans the leader line
        a.update_positions(renderer)
        return Text.get_window_extent(a, renderer)

    def _leader(p0, bb, n=16):              # points along the leader line, from the marker to the label
        p1 = np.array([np.clip(p0[0], bb.x0, bb.x1), np.clip(p0[1], bb.y0, bb.y1)])
        return [p0 + (p1 - p0) * t for t in np.linspace(0.1, 0.95, n)]

    cands = [(r * np.cos(t), r * np.sin(t)) for r in (22, 34, 50, 68, 90, 115, 145)
             for t in np.deg2rad([30, -30, 60, -60, 0, 90, -90, 135, -135, 180, 150, -150])]
    pad = np.sqrt(14 + 0.9 * scores.max()) / 2 * fig.dpi / 72      # largest marker radius, in pixels
    ranked, per_cluster = [], {}            # top-scoring terms, at most two per map cluster
    for i in fg[np.argsort(-scores)]:
        c_ = L_["cluster"].iat[i]
        if per_cluster.get(c_, 0) < 2:
            per_cluster[c_] = per_cluster.get(c_, 0) + 1
            ranked.append(i)
    for i in ranked[:n_label]:
        x0, y0 = L_["x"].iat[i], L_["y"].iat[i]
        side = -1 if x0 > cx else 1
        best = None
        for dx, dy in cands:
            dx *= side
            a = ax.annotate(_wrap(disp(L_["name"].iat[i]), 26, 3), (x0, y0), xytext=(dx, dy),
                            textcoords="offset points", ha="left" if dx > 0 else "right", va="center",
                            fontsize=8, zorder=5, path_effects=halo,
                            arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2, shrinkA=0, shrinkB=2))
            bb = _text_bbox(a)
            inside = ax_bb.x0 <= bb.x0 and bb.x1 <= ax_bb.x1 and ax_bb.y0 <= bb.y0 and bb.y1 <= ax_bb.y1
            n_pts = int(((pts_disp[:, 0] > bb.x0 - pad) & (pts_disp[:, 0] < bb.x1 + pad) &
                         (pts_disp[:, 1] > bb.y0 - pad) & (pts_disp[:, 1] < bb.y1 + pad)).sum())
            n_lab = sum(bb.padded(4).overlaps(o) for o in placed)
            # leader lines: this one must not run through a placed label, nor placed ones through this label
            lead = _leader(np.array(ax.transData.transform((x0, y0))), bb)
            n_cross = sum(o.contains(*q) for o in placed for q in lead) + \
                sum(bb.contains(*q) for q in leaders)
            cost = (0 if inside else 1e6) + (1e5 if n_lab else 0) + 40 * n_pts + 60 * n_cross + 0.05 * np.hypot(dx, dy)
            a.remove()
            if best is None or cost < best[0]:
                best = (cost, dx, dy)
            if cost == 0:
                break
        if best and best[0] < 1e5:          # inside the axes and not on another label
            a = ax.annotate(_wrap(disp(L_["name"].iat[i]), 26, 3), (x0, y0), xytext=(best[1], best[2]),
                            textcoords="offset points", ha="left" if best[1] > 0 else "right", va="center",
                            fontsize=8, zorder=5, path_effects=halo,
                            arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2, shrinkA=0, shrinkB=2))
            placed.append(_text_bbox(a))
            leaders.extend(_leader(np.array(ax.transData.transform((x0, y0))), placed[-1]))
    L = LegendColumn(fig, lax, x0=0.04)
    L.title("Map cluster: keywords\n(number of significant terms)")
    kw = L_.drop_duplicates("cluster").set_index("cluster")["cluster_name"].to_dict() \
        if "cluster_name" in L_.columns else {}
    for c in clusters:                      # every cluster, same order in every run; empty ones greyed
        n = int(n_sig.get(c, 0))
        col, m = cluster_style(c)
        yc = L.y - L._dy(0.09)
        lax.scatter([L.x0 + 0.04], [yc], s=26, marker=m, color=col if n else "white",
                    edgecolor=INK if n else NA_GREY, linewidths=0.5, clip_on=False)
        words = " · ".join(str(kw.get(c, "")).split(" · ")[:2])
        lax.text(L.x0 + 0.1, yc, f"C{c + 1} {words} ({n})", va="center", fontsize=8.5, color=INK if n else NA_GREY)
        L.y -= L._dy(0.18)
    L.y -= L._dy(0.1)
    L.note(f"{len(L_):,} {DB_LABEL.get(db, db)} terms\n(grey dots: not significant)")
    return _save(fig, path)
