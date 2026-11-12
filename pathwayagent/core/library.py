"""In-memory gene-set library (design §5, §7.2). Loaded once per process (<1 s)."""
from __future__ import annotations

import json
import re
import threading
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy import sparse

from ..config import library_path

_lock = threading.Lock()

_GREEK = [("alpha", "a"), ("beta", "b"), ("gamma", "g"), ("delta", "d"), ("kappa", "k")]


_HM_TOKENS = {"TNFA": "TNF-α", "NFKB": "NF-κB", "IL2": "IL-2", "IL6": "IL-6", "STAT5": "STAT5", "STAT3": "STAT3",
              "JAK": "JAK", "KRAS": "KRAS", "MTORC1": "mTORC1", "E2F": "E2F", "G2M": "G2/M", "MYC": "MYC",
              "V1": "v1", "V2": "v2", "P53": "p53", "UV": "UV", "DN": "down", "DNA": "DNA", "PI3K": "PI3K",
              "AKT": "AKT", "MTOR": "mTOR", "TGF": "TGF", "BETA": "β", "ALPHA": "α", "GAMMA": "γ", "WNT": "WNT",
              "NOTCH": "NOTCH", "IFN": "IFN", "ROS": "ROS", "HEME": "heme", "APICAL": "apical"}


def hallmark_name(set_id: str) -> str:
    """'HALLMARK:HALLMARK_TNFA_SIGNALING_VIA_NFKB' -> 'TNF-α signaling via NF-κB'."""
    toks = str(set_id).split(":", 1)[-1].replace("HALLMARK_", "").split("_")
    out = [_HM_TOKENS.get(t, t.lower()) for t in toks]
    s = " ".join(out).replace("TGF β", "TGF-β").replace("WNT β", "WNT/β").replace("α response", "α response")
    return s[:1].upper() + s[1:]


def compact_name(name: str) -> str:
    """Punctuation/spacing/greek-insensitive key, e.g. 'TNF-alpha Signaling via NF-kB'
    and 'TNFA_SIGNALING_VIA_NFKB' both -> 'tnfasignalingvianfkb'."""
    s = re.sub(r"[^0-9a-z]+", " ", str(name).lower())
    for a, b in _GREEK:
        s = re.sub(rf"\b{a}\b|(?<=[a-z0-9]){a}\b", b, s)
    return s.replace(" ", "")


class GeneSetLibrary:
    def __init__(self, organism: str):
        path = library_path(organism)
        if not (path / "sets.parquet").exists():
            raise FileNotFoundError(
                f"No compiled library for {organism} at {path}. Run: python -m pathwayagent db build "
                f"--organism {organism}")
        self.organism = organism
        self.path = path
        self.sets = pd.read_parquet(path / "sets.parquet")
        hm = self.sets["database"] == "HALLMARK"
        self.sets.loc[hm, "name"] = [hallmark_name(x) for x in self.sets.loc[hm, "set_id"]]
        self.matrix = sparse.load_npz(path / "membership.npz").tocsr().astype(np.float32)
        self.genes = np.array((path / "gene_index.txt").read_text().split("\n"))
        self.gene_pos = {g: i for i, g in enumerate(self.genes)}
        self.manifest = json.loads((path / "manifest.json").read_text())
        self.set_pos = {s: i for i, s in enumerate(self.sets["set_id"])}
        self.key_pos = {}
        for i, k in enumerate(self.sets["term_key"]):
            self.key_pos.setdefault(k, i)
        self.name_pos = {}
        self.compact_pos = {}
        for i, (db, nm) in enumerate(zip(self.sets["database"], self.sets["name_norm"])):
            self.name_pos.setdefault((db, nm), i)
            self.compact_pos.setdefault((db, compact_name(nm)), i)
        self._annotated = np.asarray(self.matrix.sum(axis=0)).ravel() > 0

    # ── selection helpers ──────────────────────────────────────────────
    def db_rows(self, databases) -> np.ndarray:
        return np.where(self.sets["database"].isin(list(databases)).values)[0]

    def default_universe(self, databases) -> np.ndarray:
        """Gene columns annotated in any set of the selected databases (Enrichr-style)."""
        rows = self.db_rows(databases)
        return np.where(np.asarray(self.matrix[rows].sum(axis=0)).ravel() > 0)[0]

    def gene_indices(self, symbols) -> np.ndarray:
        return np.array([self.gene_pos[g] for g in symbols if g in self.gene_pos], dtype=int)

    def set_genes(self, row: int) -> list[str]:
        return self.genes[self.matrix[row].indices].tolist()

    def set_gene_idx(self, row: int) -> np.ndarray:
        return self.matrix[row].indices

    def gmt(self, db: str) -> str:
        return str(self.path / "gmt" / f"{db}.gmt")


@lru_cache(maxsize=4)
def _load(organism: str) -> GeneSetLibrary:
    return GeneSetLibrary(organism)


def get_library(organism: str) -> GeneSetLibrary:
    with _lock:
        return _load(organism)
