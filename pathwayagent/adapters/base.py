"""Adapter contract (design §6) and the per-run context shared by all adapters."""
from __future__ import annotations

import re
import time
import traceback
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..core.idmap import GeneMapper
from ..core.library import GeneSetLibrary, compact_name
from ..core.schema import AdapterResult, AnalysisRequest, QueryInput, finalize_records
from ..db.build_library import norm_name


@dataclass
class RunContext:
    request: AnalysisRequest
    library: GeneSetLibrary
    mapper: GeneMapper
    universe_idx: np.ndarray                 # gene columns of the universe
    universe_symbols: set = field(default_factory=set)
    user_universe: bool = False              # True when the user supplied a background

    @property
    def organism(self) -> str:
        return self.request.organism


class Adapter:
    name = "base"
    family = ""
    method = ""
    mode = "native"               # native | shared_library
    web = False                   # sends gene lists to an external service
    needs_ranked = False          # requires a genome-wide ranked list
    needs_scores = False          # requires scores (any length)
    databases: tuple = ()         # database types this tool can cover
    version = ""

    def covered(self, ctx: RunContext) -> list[str]:
        return [d for d in ctx.request.databases if d in self.databases]

    def applicable(self, q: QueryInput, ctx: RunContext) -> str | None:
        """Return a reason string when the tool cannot run on this query (-> 'skipped')."""
        if self.web and ctx.request.local_only:
            return "local-only mode: external services disabled"
        if not self.covered(ctx):
            return "none of the selected databases are covered by this tool"
        if self.needs_scores and not q.is_ranked:
            return "needs a ranked list (gene + score)"
        if self.needs_ranked:
            if not q.is_ranked:
                return "needs a ranked list (gene + score)"
            n = len([g for g in q.genes if g in ctx.universe_symbols])
            if n < ctx.request.min_ranked_for_gsea:
                return (f"ranked list covers only {n} genes; GSEA-type tests need a genome-wide ranking "
                        f"(>= {ctx.request.min_ranked_for_gsea} genes). Rank-based ORA (mHG) is used instead.")
        return None

    def run(self, q: QueryInput, ctx: RunContext) -> AdapterResult:
        res = self._run_safe(q, ctx)
        res.method = self.method
        return res

    def _run_safe(self, q: QueryInput, ctx: RunContext) -> AdapterResult:
        t0 = time.time()
        reason = self.applicable(q, ctx)
        if reason:
            return AdapterResult(self.name, q.name, "skipped", reason, 0.0, family=self.family)
        try:
            df, extra = self._run(q, ctx)
            df = finalize_records(df)
            status = "ok" if len(df) else "empty"
            return AdapterResult(self.name, q.name, status, "", time.time() - t0, df,
                                 databases=self.covered(ctx), family=self.family,
                                 version=self.version, extra=extra or {})
        except Exception as exc:  # every failure is surfaced, never swallowed (bug B13)
            msg = f"{type(exc).__name__}: {exc}"
            return AdapterResult(self.name, q.name, "error", msg, time.time() - t0, family=self.family,
                                 extra={"traceback": traceback.format_exc(limit=3)})

    def _run(self, q: QueryInput, ctx: RunContext) -> tuple[pd.DataFrame, dict]:
        raise NotImplementedError

    # ── helpers ────────────────────────────────────────────────────────
    def base_record(self, q: QueryInput, **kw) -> dict:
        rec = {"query": q.name, "tool": self.name, "method": self.method, "family": self.family,
               "mode": self.mode, "padj_method": "BH"}
        rec.update(kw)
        return rec


# ── term keys: map any tool's term id / name onto the library's canonical term_key ──

def key_from_id(term_id: str) -> str | None:
    t = str(term_id).strip()
    if m := re.search(r"(GO:\d{7})", t):
        return m.group(1)
    if m := re.search(r"R-[A-Z]{3}-(\d+)", t) or re.fullmatch(r"(?:REAC:)?[A-Z]{3}-(\d+)", t):
        return f"REACT:{m.group(1)}"
    if m := re.search(r"\bWP(\d+)\b", t):
        return f"WP:WP{m.group(1)}"
    if m := re.fullmatch(r"(?:KEGG:|path:)?(?:[a-z]{3,4})?(\d{5})", t):
        return f"KEGG:{m.group(1)}"
    return None


def resolve_key(lib: GeneSetLibrary, database: str, term_id: str = "", name: str = "") -> tuple[str, int | None]:
    """Return (term_key, library row or None). Exact id first, then name within the database."""
    k = key_from_id(term_id) if term_id else None
    if k is None and name:
        k = key_from_id(name) if database in ("GO_BP", "GO_MF", "GO_CC", "REACTOME", "WIKIPATHWAYS") else None
    if k and k in lib.key_pos:
        return k, lib.key_pos[k]
    nm = norm_name(re.sub(r"\s*\((GO:\d+)\)\s*$|\s+R-[A-Z]{3}-\d+$|\s+WP\d+$", "", name))
    if nm:
        row = lib.name_pos.get((database, nm))
        if row is None:
            row = lib.compact_pos.get((database, compact_name(nm)))
        if row is not None:
            return lib.sets["term_key"].iat[row], row
    if k:
        return k, None
    return f"{database}:{nm}", None
