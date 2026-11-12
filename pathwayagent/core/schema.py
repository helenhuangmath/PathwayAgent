"""Common data structures (design §8)."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import pandas as pd

from ..config import DEFAULT_DATABASES, DEFAULTS, PRESETS

# One row = one term tested by one tool with one method on one query list.
RECORD_COLUMNS = [
    "query", "tool", "method", "family", "mode", "database",
    "term_id", "term_key", "term_name", "term_size",
    "query_size", "universe_size", "overlap",
    "pval", "padj", "padj_method", "effect", "effect_type", "direction",
    "hit_genes", "rank", "n_tested",
]


def empty_records() -> pd.DataFrame:
    return pd.DataFrame(columns=RECORD_COLUMNS)


def finalize_records(df: pd.DataFrame) -> pd.DataFrame:
    """Fill missing columns, coerce types and add within-(tool, database) ranks."""
    if df is None or df.empty:
        return empty_records()
    df = df.copy()
    for c in RECORD_COLUMNS:
        if c not in df.columns:
            df[c] = None
    for c in ("pval", "padj"):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(1.0).clip(1e-300, 1.0)
    for c in ("overlap", "term_size", "query_size", "universe_size"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["effect"] = pd.to_numeric(df["effect"], errors="coerce")
    df["direction"] = pd.to_numeric(df["direction"], errors="coerce").fillna(0).astype(int)
    df["hit_genes"] = df["hit_genes"].apply(
        lambda x: ";".join(x) if isinstance(x, (list, tuple, set)) else ("" if x is None or x != x else str(x)))
    grp = df.groupby(["query", "tool", "method", "database"], sort=False)
    df["rank"] = grp["pval"].rank(method="min")
    df["n_tested"] = grp["pval"].transform("size")
    return df[RECORD_COLUMNS].reset_index(drop=True)


@dataclass
class QueryInput:
    """One gene list (optionally with scores) to analyse."""
    name: str
    genes: list[str]
    scores: Optional[dict[str, float]] = None      # gene -> score (higher = more "up")
    raw_genes: list[str] = field(default_factory=list)
    mapping: Optional[pd.DataFrame] = None           # ID-mapping report

    @property
    def is_ranked(self) -> bool:
        return bool(self.scores)

    def ranked_genes(self) -> list[str]:
        if not self.scores:
            return []
        return [g for g, _ in sorted(self.scores.items(), key=lambda kv: -kv[1])]


@dataclass
class AnalysisRequest:
    queries: list[QueryInput]
    organism: str = "mmusculus"
    universe: Optional[list[str]] = None
    databases: list[str] = field(default_factory=lambda: list(DEFAULT_DATABASES))
    tools: list[str] = field(default_factory=lambda: list(PRESETS["standard"]))
    min_size: int = DEFAULTS["min_size"]
    max_size: int = DEFAULTS["max_size"]
    alpha: float = DEFAULTS["alpha"]
    deadline_s: float = DEFAULTS["deadline_s"]
    local_only: bool = False
    min_ranked_for_gsea: int = 1000
    profile: Optional[dict] = None
    metascape: Optional[dict[str, str]] = None      # list name (or "*") -> Metascape result package
    title: str = "PathwayAgent run"


@dataclass
class AdapterResult:
    tool: str
    query: str
    status: str                     # ok | empty | error | timeout | skipped
    message: str = ""
    runtime_s: float = 0.0
    records: pd.DataFrame = field(default_factory=empty_records)
    databases: list[str] = field(default_factory=list)   # databases this tool covered
    family: str = ""
    version: str = ""
    method: str = ""
    extra: dict = field(default_factory=dict)             # e.g. native tables for figures

    def summary(self) -> dict:
        n_sig = int((self.records["padj"] < 0.05).sum()) if len(self.records) else 0
        return {"tool": self.tool, "query": self.query, "status": self.status, "message": self.message,
                "runtime_s": round(self.runtime_s, 2), "n_tested": int(len(self.records)),
                "n_significant": n_sig, "databases": self.databases, "family": self.family,
                "method": self.method, "version": self.version}
