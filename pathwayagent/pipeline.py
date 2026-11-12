"""
Orchestrator (design §13): parse -> map IDs -> run tools in parallel under a deadline
-> harmonise -> consensus -> themes -> key genes -> figures & reports.
"""
from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import pandas as pd

from .adapters import PLANNED, REGISTRY
from .adapters.base import RunContext
from .config import RESULTS_DIR
from .consensus.keygenes import key_genes
from .consensus.score import consensus, stability
from .consensus.themes import build_themes
from .core.idmap import get_mapper, guess_organism
from .core.library import get_library
from .core.schema import AdapterResult, AnalysisRequest, QueryInput, empty_records


# tool -> tool it exactly reproduces (shared-library mode); reported per tool, counted once in consensus
REPLICATES = {"clusterprofiler": "local_ora"}


def _replicate_note(records: pd.DataFrame, rep: str, ref: str) -> str:
    a = records[records["tool"] == ref].set_index(["database", "term_key"])["pval"]
    b = records[records["tool"] == rep].set_index(["database", "term_key"])["pval"]
    common = a.index.intersection(b.index)
    if len(common):
        diff = float(np.max(np.abs(np.log10(a.loc[common].values) - np.log10(b.loc[common].values))))
        agree = f"p-values agree on {len(common)} shared terms (max |Δlog10 p| = {diff:.1e})"
    else:
        agree = "no shared terms to compare"
    return (f"{rep} runs the same hypergeometric test on the same gene sets as {ref}: {agree}. "
            f"It is reported per tool as a cross-check but counted once in the consensus.")


@dataclass
class QueryResult:
    query: QueryInput
    tool_results: list[AdapterResult] = field(default_factory=list)
    records: pd.DataFrame = field(default_factory=empty_records)
    consensus: pd.DataFrame = field(default_factory=pd.DataFrame)
    themes: pd.DataFrame = field(default_factory=pd.DataFrame)
    key_genes: pd.DataFrame = field(default_factory=pd.DataFrame)
    figures: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


@dataclass
class RunResult:
    request: AnalysisRequest
    run_id: str
    outdir: Path
    queries: list[QueryResult]
    timings: dict
    universe_size: int
    universe_source: str
    library_manifest: dict
    warnings: list[str] = field(default_factory=list)


def _emit(cb, **ev):
    if cb:
        try:
            cb(ev)
        except Exception:
            pass


def prepare(request: AnalysisRequest) -> tuple[RunContext, list[str]]:
    warnings = []
    if request.organism in (None, "", "auto"):
        request.organism = guess_organism([g for q in request.queries for g in q.genes])
    lib = get_library(request.organism)
    mapper = get_mapper(request.organism)
    for q in request.queries:
        mapped, report = mapper.map(q.raw_genes or q.genes)
        q.mapping = report
        if q.scores:
            sym = dict(zip(report["input"], report["symbol"]))
            sc = {}
            for g, s in q.scores.items():
                m = sym.get(g)
                if m and m not in sc:
                    sc[m] = s
            q.scores = sc
            mapped = [g for g, _ in sorted(sc.items(), key=lambda kv: -kv[1])]
        q.genes = mapped
        n_un = int((report["symbol"].isna()).sum())
        if n_un:
            warnings.append(f"[{q.name}] {n_un}/{len(report)} input genes could not be mapped "
                            f"(see id_mapping.tsv)")
    if request.universe:
        u_sym, _ = mapper.map(request.universe)
        uni = lib.gene_indices(u_sym)
        user_uni = True
    else:
        uni = lib.default_universe(request.databases)
        user_uni = False
        warnings.append("No background universe supplied: using all genes annotated in the selected "
                        "databases (Enrichr-style). Supplying the detected/expressed genes is recommended.")
    ctx = RunContext(request, lib, mapper, uni, set(lib.genes[uni].tolist()), user_uni)
    return ctx, warnings


def run_analysis(request: AnalysisRequest, outdir: str | Path | None = None,
                 on_event: Optional[Callable[[dict], None]] = None, make_outputs: bool = True) -> RunResult:
    t0 = time.time()
    timings = {}
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(outdir) if outdir else RESULTS_DIR / run_id
    if make_outputs:
        outdir.mkdir(parents=True, exist_ok=True)

    ctx, warnings = prepare(request)
    timings["prepare_s"] = round(time.time() - t0, 2)
    _emit(on_event, type="stage", stage="prepared", organism=request.organism,
          universe=len(ctx.universe_idx), queries=[q.name for q in request.queries])

    unknown = [t for t in request.tools if t not in REGISTRY]
    for t in unknown:
        warnings.append(f"Tool '{t}' not available: {PLANNED.get(t, 'unknown tool')}")
    tools = [t for t in request.tools if t in REGISTRY]

    results: dict[tuple, AdapterResult] = {}
    t_tools = time.time()
    pool = ThreadPoolExecutor(max_workers=12)
    futs = {pool.submit(REGISTRY[t]().run, q, ctx): (q.name, t) for q in request.queries for t in tools}
    _emit(on_event, type="stage", stage="tools_started", n_tasks=len(futs))
    pending = set(futs)
    deadline = t_tools + request.deadline_s
    while pending:
        remaining = deadline - time.time()
        if remaining <= 0:
            break
        done, pending = wait(pending, timeout=remaining, return_when=FIRST_COMPLETED)
        for f in done:
            qn, t = futs[f]
            res = f.result()
            results[(qn, t)] = res
            _emit(on_event, type="tool_done", **res.summary())
    for f in pending:
        qn, t = futs[f]
        res = AdapterResult(t, qn, "timeout", f"no answer within the {request.deadline_s:.0f}s deadline",
                            request.deadline_s)
        results[(qn, t)] = res
        _emit(on_event, type="tool_done", **res.summary())
    pool.shutdown(wait=False, cancel_futures=True)
    timings["tools_s"] = round(time.time() - t_tools, 2)

    t_cons = time.time()
    qresults = []
    interest = set((request.profile or {}).get("genes_of_interest", []) or [])
    for q in request.queries:
        qr = QueryResult(q, [results[(q.name, t)] for t in tools])
        frames = [r.records for r in qr.tool_results if r.status == "ok" and len(r.records)]
        qr.records = pd.concat(frames, ignore_index=True) if frames else empty_records()
        tool_info = {r.tool: {"family": r.family, "databases": r.databases}
                     for r in qr.tool_results if r.status in ("ok", "empty")}
        # Replicates: same test on the same gene sets must not be counted twice (design §1.2).
        for rep, ref in REPLICATES.items():
            if rep in tool_info and ref in tool_info:
                tool_info.pop(rep)
                qr.warnings.append(_replicate_note(qr.records, rep, ref))
        if len(qr.records):
            c = consensus(qr.records, tool_info, ctx.library, request.alpha)
            c = stability(qr.records, tool_info, ctx.library, c, request.alpha)
            c, th = build_themes(c)
            qr.consensus, qr.themes = c, th
            qr.key_genes = key_genes(qr.records, c, th, q, request.alpha, interest)
            if qr.key_genes is not None and len(qr.key_genes) and qr.key_genes.attrs.get("fallback"):
                qr.warnings.append("No ROBUST/SUPPORTED themes: key genes are drawn from single-tool themes "
                                   "and should be treated as exploratory.")
        n_ran = len(tool_info)
        if n_ran < 2:
            qr.warnings.append("Fewer than two tools returned results: no cross-tool consensus possible.")
        qresults.append(qr)
        _emit(on_event, type="query_done", query=q.name,
              n_robust=int((qr.consensus.get("tier", pd.Series(dtype=int)) == 1).sum()) if len(qr.consensus) else 0)
    timings["consensus_s"] = round(time.time() - t_cons, 2)

    run = RunResult(request, run_id, outdir, qresults, timings, len(ctx.universe_idx),
                    "user-supplied" if ctx.user_universe else "annotated genes of selected databases",
                    ctx.library.manifest, warnings)
    if make_outputs:
        t_out = time.time()
        from .report.writer import write_outputs
        write_outputs(run, ctx, on_event)
        timings["outputs_s"] = round(time.time() - t_out, 2)
    timings["total_s"] = round(time.time() - t0, 2)
    if make_outputs:
        from .report.writer import write_manifest
        write_manifest(run)
    _emit(on_event, type="done", run_id=run_id, outdir=str(outdir), timings=timings)
    return run
