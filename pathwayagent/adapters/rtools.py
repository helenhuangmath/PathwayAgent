"""
R-based tools through one persistent R worker (design §7.3):
clusterProfiler (shared library + native enrichGO), ReactomePA (native), fgsea (shared library).
"""
from __future__ import annotations

import json
import queue
import subprocess
import threading
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import RSCRIPT
from ..core.schema import QueryInput
from .base import Adapter, RunContext, resolve_key

WORKER_R = Path(__file__).resolve().parent.parent / "rworker" / "worker.R"


class RWorker:
    """A long-lived Rscript process; requests are serialised with a lock."""

    def __init__(self, rscript: str = RSCRIPT, startup_timeout: float = 90):
        self.rscript = rscript
        self.startup_timeout = startup_timeout
        self.proc = None
        self.versions = {}
        self._lock = threading.Lock()
        self._out: queue.Queue = queue.Queue()
        self._stderr: list[str] = []
        self._n = 0

    def _reader(self, stream, sink):
        for line in iter(stream.readline, ""):
            sink(line)

    def start(self):
        if self.proc and self.proc.poll() is None:
            return
        self._out = queue.Queue()
        self.proc = subprocess.Popen([self.rscript, "--vanilla", str(WORKER_R)], stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
        threading.Thread(target=self._reader, args=(self.proc.stdout, self._on_stdout), daemon=True).start()
        threading.Thread(target=self._reader, args=(self.proc.stderr, self._stderr.append), daemon=True).start()
        msg = self._wait(self.startup_timeout)
        self.versions = msg.get("versions", {})

    def _on_stdout(self, line: str):
        if line.startswith("@@JSON@@"):
            self._out.put(json.loads(line[len("@@JSON@@"):]))

    def _wait(self, timeout: float) -> dict:
        try:
            return self._out.get(timeout=timeout)
        except queue.Empty:
            self.stop()
            raise TimeoutError(f"R worker did not answer within {timeout:.0f}s; stderr: "
                               f"{''.join(self._stderr[-5:])[-500:]}")

    def call(self, task: str, args: dict, timeout: float = 120) -> dict:
        with self._lock:
            self.start()
            self._n += 1
            req = {"id": self._n, "task": task, "args": args}
            self.proc.stdin.write(json.dumps(req) + "\n")
            self.proc.stdin.flush()
            msg = self._wait(timeout)
        if not msg.get("ok"):
            raise RuntimeError(f"R task {task} failed: {msg.get('error')}")
        return msg["result"]

    def stop(self):
        if self.proc and self.proc.poll() is None:
            self.proc.kill()
        self.proc = None


N_WORKERS = 2
_pool: queue.Queue | None = None
_all_workers: list[RWorker] = []
_wlock = threading.Lock()


class _Lease:
    """Borrow a free worker from the pool for one call (so slow tasks don't block fast ones)."""

    def __init__(self):
        global _pool
        with _wlock:
            if _pool is None:
                _pool = queue.Queue()
                for _ in range(N_WORKERS):
                    w = RWorker()
                    _all_workers.append(w)
                    _pool.put(w)

    @property
    def versions(self) -> dict:
        return next((w.versions for w in _all_workers if w.versions), {})

    def call(self, task: str, args: dict, timeout: float = 120) -> dict:
        w = _pool.get()
        try:
            return w.call(task, args, timeout)
        finally:
            _pool.put(w)

    def start(self):
        for w in _all_workers:
            threading.Thread(target=w.start, daemon=True).start()


def get_worker() -> _Lease:
    return _Lease()


def r_available() -> bool:
    return Path(RSCRIPT).exists()


def _ratio(s) -> float:
    try:
        a, b = str(s).split("/")
        return int(a) / max(int(b), 1)
    except Exception:
        return np.nan


def _enrich_table_to_records(adapter: Adapter, q: QueryInput, ctx: RunContext, tables: dict,
                             ids_are_library: bool) -> pd.DataFrame:
    lib = ctx.library
    recs = []
    for db, rows in tables.items():
        if not rows:
            continue
        df = pd.DataFrame(rows).rename(columns={"p.adjust": "p_adjust"})
        for r in df.itertuples(index=False):
            if ids_are_library and r.ID in lib.set_pos:
                row = lib.set_pos[r.ID]
                key, name = lib.sets["term_key"].iat[row], lib.sets["name"].iat[row]
            else:
                key, _ = resolve_key(lib, db, r.ID, r.Description)
                name = r.Description
            gr, br = _ratio(r.GeneRatio), _ratio(r.BgRatio)
            n_q = int(str(r.GeneRatio).split("/")[1])
            n_u = int(str(r.BgRatio).split("/")[1])
            recs.append(adapter.base_record(
                q, database=db, term_id=r.ID, term_key=key, term_name=name,
                term_size=int(str(r.BgRatio).split("/")[0]), query_size=n_q, universe_size=n_u,
                overlap=int(r.Count), pval=r.pvalue, padj=r.p_adjust, effect=gr / br if br else np.nan, effect_type="fold_enrichment", direction=0,
                hit_genes=str(r.geneID).split("/") if isinstance(r.geneID, str) else []))
    return pd.DataFrame(recs)


class RAdapter(Adapter):
    r_task = ""
    timeout = 120

    def applicable(self, q, ctx):
        if not r_available():
            return f"Rscript not found at {RSCRIPT}"
        return super().applicable(q, ctx)

    def _universe(self, ctx):
        return sorted(ctx.universe_symbols) if ctx.user_universe else []


class ClusterProfiler(RAdapter):
    """clusterProfiler::enricher on the SAME gene sets as the local engine (method-only comparison)."""
    name, family, method, mode = "clusterprofiler", "F1_ORA", "hypergeometric (clusterProfiler::enricher)", "shared_library"
    databases = ("GO_BP", "GO_MF", "GO_CC", "KEGG", "REACTOME", "WIKIPATHWAYS", "HALLMARK", "DISGENET")

    def _run(self, q: QueryInput, ctx: RunContext):
        w = get_worker()
        args = {"genes": q.genes, "universe": sorted(ctx.universe_symbols),
                "gmts": {db: ctx.library.gmt(db) for db in self.covered(ctx)},
                "min_size": ctx.request.min_size, "max_size": ctx.request.max_size}
        res = w.call("enricher", args, self.timeout)
        self.version = f"clusterProfiler {w.versions.get('clusterProfiler', '')}"
        return _enrich_table_to_records(self, q, ctx, res, ids_are_library=True), {}


class ClusterProfilerGO(RAdapter):
    """clusterProfiler with its own org.*.eg.db GO annotation (GOALL, as enrichGO uses; native mode).
    The propagated table is built once and cached to library/<org>/rcache (enrichGO rebuilds it per call)."""
    name, family, method, mode = "clusterprofiler_go", "F1_ORA", "hypergeometric (clusterProfiler, org.db GO)", "native"
    timeout = 300
    databases = ("GO_BP", "GO_MF", "GO_CC")

    def _run(self, q: QueryInput, ctx: RunContext):
        w = get_worker()
        args = {"genes": q.genes, "universe": self._universe(ctx), "organism": ctx.organism,
                "cache_dir": str(ctx.library.path / "rcache"),
                "onts": [d.split("_")[1] for d in self.covered(ctx)],
                "min_size": ctx.request.min_size, "max_size": ctx.request.max_size}
        res = w.call("enrichgo", args, self.timeout)
        self.version = f"clusterProfiler {w.versions.get('clusterProfiler', '')} + org db"
        return _enrich_table_to_records(self, q, ctx, res, ids_are_library=False), {}


class ReactomePA(RAdapter):
    name, family, method, mode = "reactomepa", "F1_ORA", "hypergeometric (ReactomePA::enrichPathway)", "native"
    databases = ("REACTOME",)

    def _run(self, q: QueryInput, ctx: RunContext):
        w = get_worker()
        args = {"genes": q.genes, "universe": self._universe(ctx), "organism": ctx.organism,
                "min_size": ctx.request.min_size, "max_size": ctx.request.max_size}
        res = w.call("reactomepa", args, self.timeout)
        self.version = "ReactomePA + reactome.db"
        return _enrich_table_to_records(self, q, ctx, res, ids_are_library=False), {}


class FGSEA(RAdapter):
    name, family, method, mode = "fgsea", "F2_GSEA", "fgseaMultilevel", "shared_library"
    databases = ("KEGG", "REACTOME", "WIKIPATHWAYS", "HALLMARK", "GO_BP", "GO_MF", "GO_CC")
    needs_ranked = True

    def _run(self, q: QueryInput, ctx: RunContext):
        w = get_worker()
        genes = [g for g in q.genes if g in ctx.universe_symbols]
        args = {"genes": genes, "scores": [q.scores[g] for g in genes],
                "gmts": {db: ctx.library.gmt(db) for db in self.covered(ctx)},
                "min_size": ctx.request.min_size, "max_size": ctx.request.max_size}
        res = w.call("fgsea", args, self.timeout)
        self.version = f"fgsea {w.versions.get('fgsea', '')}"
        lib = ctx.library
        recs = []
        for db, rows in res.items():
            for r in rows:
                row = lib.set_pos.get(r["pathway"])
                recs.append(self.base_record(
                    q, database=db, term_id=r["pathway"],
                    term_key=lib.sets["term_key"].iat[row] if row is not None else r["pathway"],
                    term_name=lib.sets["name"].iat[row] if row is not None else r["pathway"],
                    term_size=r["size"], overlap=r["size"], query_size=len(genes),
                    pval=r["pval"], padj=r["padj"], effect=r["NES"], effect_type="NES",
                    direction=int(np.sign(r["NES"] or 0)), hit_genes=str(r["leadingEdge"]).split(";")))
        return pd.DataFrame(recs), {}
