"""
Term landscape (Enrichr-appyter style): a 2-D map of every gene set in a database, so a
run's significant terms can be shown in the context of the whole library.

TF-IDF of the set x gene matrix -> TruncatedSVD(50) -> t-SNE(2D); clusters = k-means on the SVD
space. Computed once per database at library build: library/<org>/landscape_<DB>.parquet
"""
from __future__ import annotations

import math
import re
import time
from collections import Counter

import numpy as np
import pandas as pd

from ..config import library_path

LANDSCAPE_DBS = ("GO_BP", "REACTOME", "KEGG", "WIKIPATHWAYS", "HALLMARK")
_GENERIC = {"regulation", "positive", "negative", "process", "response", "pathway", "activity", "cell", "cellular",
            "involved", "signaling", "other", "events", "mediated", "the", "and", "via", "from", "for", "with",
            "into"}


def cluster_keywords(df: pd.DataFrame, k: int = 3) -> dict[int, str]:
    """Name each map cluster by the k words most specific to its term names (TF-IDF over clusters)."""
    words = {c: Counter(w for n in g for w in re.findall(r"[a-z][a-z0-9\-]{2,}", str(n).lower())
                        if w not in _GENERIC)
             for c, g in df.groupby("cluster")["name"]}
    n_docs = Counter(w for cnt in words.values() for w in cnt)
    out = {}
    for c, cnt in words.items():
        score = {w: (1 + math.log(f)) * (math.log((1 + len(words)) / (1 + n_docs[w])) + 1) for w, f in cnt.items()}
        out[int(c)] = " · ".join(sorted(score, key=lambda w: (-score[w], w))[:k])
    return out


def build_landscape(organism: str, db: str, n_clusters: int | None = None, seed: int = 0) -> pd.DataFrame:
    from sklearn.cluster import KMeans
    from sklearn.decomposition import TruncatedSVD
    from sklearn.feature_extraction.text import TfidfTransformer
    from sklearn.manifold import TSNE

    from ..core.library import get_library

    lib = get_library(organism)
    rows = lib.db_rows([db])
    X = lib.matrix[rows]
    X = X[:, np.asarray(X.sum(axis=0)).ravel() > 0]
    Xt = TfidfTransformer().fit_transform(X)
    k = min(50, X.shape[0] - 1, X.shape[1] - 1)
    Z = TruncatedSVD(n_components=k, random_state=seed).fit_transform(Xt)
    Z = Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-12)
    perp = max(5, min(30, (len(rows) - 1) // 3))
    xy = TSNE(n_components=2, perplexity=perp, init="pca", random_state=seed, metric="cosine").fit_transform(Z)
    n_clusters = n_clusters or int(np.clip(round(np.sqrt(len(rows)) / 3.5), 4, 24))
    cl = KMeans(n_clusters=n_clusters, n_init=4, random_state=seed).fit_predict(Z)
    df = pd.DataFrame({"set_id": lib.sets["set_id"].values[rows], "term_key": lib.sets["term_key"].values[rows],
                       "name": lib.sets["name"].values[rows], "x": xy[:, 0], "y": xy[:, 1], "cluster": cl})
    df["cluster_name"] = df["cluster"].map(cluster_keywords(df))
    df.to_parquet(library_path(organism) / f"landscape_{db}.parquet", index=False)
    return df


def build_all(organism: str, dbs=LANDSCAPE_DBS):
    for db in dbs:
        t0 = time.time()
        df = build_landscape(organism, db)
        print(f"  landscape {organism} {db}: {len(df)} terms, {df['cluster'].nunique()} clusters "
              f"({time.time() - t0:.0f}s)")


def load_landscape(organism: str, db: str) -> pd.DataFrame | None:
    p = library_path(organism) / f"landscape_{db}.parquet"
    if not p.exists():
        return None
    df = pd.read_parquet(p)
    if "cluster_name" not in df.columns:            # maps built before clusters were named
        df["cluster_name"] = df["cluster"].map(cluster_keywords(df))
    return df
