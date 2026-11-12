"""
Local gene ID mapping (design §3.2): symbols (any case), synonyms/previous symbols,
Entrez and Ensembl IDs -> official symbols; mouse <-> human orthologs.
Unmapped genes are reported, never silently dropped.
"""
from __future__ import annotations

import re
from functools import lru_cache

import pandas as pd

from ..config import library_path


class GeneMapper:
    def __init__(self, organism: str):
        self.organism = organism
        g = pd.read_parquet(library_path(organism) / "genes.parquet")
        self.table = g
        self.official = set(g["symbol"])
        self.upper = {}
        for s in g["symbol"]:
            self.upper.setdefault(s.upper(), s)
        self.syn: dict[str, set] = {}
        for s, syns in zip(g["symbol"], g["synonyms"].fillna("")):
            for a in str(syns).split("|"):
                a = a.strip().upper()
                if a and a not in self.upper:
                    self.syn.setdefault(a, set()).add(s)
        self.entrez = {str(e): s for e, s in zip(g["entrez"], g["symbol"]) if isinstance(e, str) and e}
        self.ensembl = {str(e): s for e, s in zip(g["ensembl"], g["symbol"]) if isinstance(e, str) and e}
        self.ortholog = {s: o.split("|") for s, o in zip(g["symbol"], g["ortholog"].fillna("")) if o}
        self.entrez_of = {s: e for e, s in self.entrez.items()}

    def map_one(self, x: str) -> tuple[str | None, str]:
        x = str(x).strip()
        if x in self.official:
            return x, "exact"
        if re.fullmatch(r"ENS[A-Z]*G\d+(\.\d+)?", x):
            s = self.ensembl.get(x.split(".")[0])
            return (s, "ensembl") if s else (None, "unmapped")
        if x.isdigit():
            s = self.entrez.get(x)
            return (s, "entrez") if s else (None, "unmapped")
        u = x.upper()
        if u in self.upper:
            return self.upper[u], "case"
        cands = self.syn.get(u)
        if cands:
            if len(cands) == 1:
                return next(iter(cands)), "synonym"
            return None, "ambiguous:" + "|".join(sorted(cands))
        return None, "unmapped"

    def map(self, genes: list[str]) -> tuple[list[str], pd.DataFrame]:
        rows, out, seen = [], [], set()
        for g in genes:
            s, how = self.map_one(g)
            rows.append({"input": g, "symbol": s, "status": how})
            if s and s not in seen:
                seen.add(s)
                out.append(s)
        return out, pd.DataFrame(rows)

    def to_orthologs(self, symbols: list[str]) -> list[str]:
        out, seen = [], set()
        for s in symbols:
            for o in self.ortholog.get(s, []):
                if o and o not in seen:
                    seen.add(o)
                    out.append(o)
        return out

    def to_entrez(self, symbols: list[str]) -> list[str]:
        return [self.entrez_of[s] for s in symbols if s in self.entrez_of]


@lru_cache(maxsize=4)
def get_mapper(organism: str) -> GeneMapper:
    return GeneMapper(organism)


def guess_organism(genes: list[str]) -> str:
    """Mouse symbols are mostly Capitalised-lowercase (Tox, Pdcd1); human are UPPER (TOX)."""
    alpha = [g for g in genes if re.search(r"[A-Za-z]", g) and not g.startswith("ENS")]
    if not alpha:
        ens = [g for g in genes if g.startswith("ENSMUSG")]
        return "mmusculus" if ens else "hsapiens"
    mouse_like = sum(bool(re.search(r"[a-z]", g)) for g in alpha)
    return "mmusculus" if mouse_like / len(alpha) > 0.5 else "hsapiens"
