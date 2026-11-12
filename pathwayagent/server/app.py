"""
Local web app + REST API (design §13.1-13.2).

  pathwayagent serve            -> http://127.0.0.1:8765
  POST /api/v1/jobs             JSON AnalysisRequest-like body -> {job_id}
  POST /api/v1/jobs/upload      multipart: file + params (JSON string)
  GET  /api/v1/jobs/{id}        status, per-tool progress, results summary when done
  GET  /api/v1/jobs/{id}/events Server-Sent Events stream of progress
  GET  /api/v1/options          organisms, databases, tools, presets
  GET  /api/v1/health           also readable cross-origin (used by the docs site's launch page)
  GET  /results/<run>/...       every output file (report.html, TSV, figures)
"""
from __future__ import annotations

import asyncio
import json
import threading
import time
import uuid
import webbrowser
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Optional

import pandas as pd
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from .. import __version__
from ..adapters import PLANNED, REGISTRY
from ..config import DATABASES, DEFAULT_DATABASES, DEFAULTS, FAMILIES, ORGANISMS, PRESETS, RESULTS_DIR
from ..core.inputs import parse_text
from ..core.schema import AnalysisRequest, QueryInput
from ..pipeline import run_analysis

STATIC = Path(__file__).resolve().parent / "static"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

app = FastAPI(title="PathwayAgent", version=__version__)
app.mount("/results", StaticFiles(directory=str(RESULTS_DIR)), name="results")
app.mount("/examples", StaticFiles(directory=str(STATIC / "examples")), name="examples")
EXEC = ThreadPoolExecutor(max_workers=2)
JOBS: dict[str, dict] = {}


class ListIn(BaseModel):
    name: str
    genes: list[str]
    scores: Optional[list[float]] = None


class JobIn(BaseModel):
    text: Optional[str] = None            # pasted genes / table (auto-detected)
    lists: Optional[list[ListIn]] = None  # or explicit lists
    name: str = "query"
    organism: str = "auto"
    preset: Optional[str] = "standard"
    tools: Optional[list[str]] = None
    databases: list[str] = list(DEFAULT_DATABASES)
    universe: Optional[str] = None
    min_size: int = DEFAULTS["min_size"]
    max_size: int = DEFAULTS["max_size"]
    alpha: float = DEFAULTS["alpha"]
    deadline_s: float = DEFAULTS["deadline_s"]
    local_only: bool = False
    genes_of_interest: Optional[list[str]] = None
    background_note: Optional[str] = None


def _build_request(job: JobIn) -> AnalysisRequest:
    queries: list[QueryInput] = []
    if job.lists:
        for l in job.lists:
            sc = dict(zip(l.genes, l.scores)) if l.scores else None
            queries.append(QueryInput(name=l.name, genes=list(l.genes), scores=sc, raw_genes=list(l.genes)))
    if job.text and job.text.strip():
        queries.extend(parse_text(job.text, name=job.name or "query"))
    if not queries:
        raise HTTPException(400, "No genes supplied.")
    tools = job.tools or PRESETS.get(job.preset or "standard", PRESETS["standard"])
    universe = None
    if job.universe and job.universe.strip():
        universe = [x for x in job.universe.replace(",", "\n").split() if x.strip()]
    return AnalysisRequest(queries=queries, organism=job.organism, universe=universe,
                           databases=job.databases, tools=tools, min_size=job.min_size, max_size=job.max_size,
                           alpha=job.alpha, deadline_s=job.deadline_s, local_only=job.local_only,
                           profile={"genes_of_interest": job.genes_of_interest or [],
                                    "background": job.background_note or ""})


def _records(df: pd.DataFrame, cols: list[str], n: int) -> list[dict]:
    if df is None or df.empty:
        return []
    cols = [c for c in cols if c in df.columns]
    return json.loads(df[cols].head(n).to_json(orient="records"))


def _summarise(run) -> dict:
    out = {"run_id": run.run_id, "report": f"/results/{run.outdir.name}/report.html",
           "summary_txt": f"/results/{run.outdir.name}/summary.txt", "timings": run.timings,
           "warnings": run.warnings, "organism": run.request.organism, "universe_size": run.universe_size,
           "queries": []}
    root = RESULTS_DIR.resolve()
    for qr in run.queries:
        folder = "".join(c if c.isalnum() or c in "-_." else "_" for c in qr.query.name)
        base = f"/results/{run.outdir.name}/{folder}"
        figs = {k: "/results/" + Path(v).resolve().relative_to(root).as_posix()
                for k, v in qr.figures.items() if isinstance(v, str)}
        out["queries"].append({
            "name": qr.query.name, "n_genes": len(qr.query.genes), "warnings": qr.warnings,
            "tools": [r.summary() for r in qr.tool_results],
            "consensus": _records(qr.consensus[qr.consensus["tier"] > 0] if len(qr.consensus) else qr.consensus,
                                  ["rank", "tier_label", "term_name", "database", "consensus_score", "tool_support",
                                   "best_padj", "hit_genes", "url"], 60),
            "themes": _records(qr.themes, ["theme_rank", "tier_label", "theme_name", "theme_score", "n_terms",
                                           "databases", "core_genes"], 20),
            "key_genes": _records(qr.key_genes, ["rank", "gene", "gene_score", "n_themes", "n_tools",
                                                 "driver_of", "why"], 30),
            "figures": figs, "network": qr.figures.get("_cy"),
            "folder": base,
        })
    return out


def _run_job(job_id: str, req: AnalysisRequest):
    job = JOBS[job_id]

    def on_event(ev):
        ev["t"] = round(time.time() - job["started"], 2)
        job["events"].append(ev)

    try:
        job["status"] = "running"
        run = run_analysis(req, outdir=RESULTS_DIR / f"web_{job_id}", on_event=on_event)
        job["result"] = _summarise(run)
        job["status"] = "done"
    except Exception as exc:
        job["status"] = "error"
        job["error"] = f"{type(exc).__name__}: {exc}"
        job["events"].append({"type": "error", "message": job["error"]})


def _submit(req: AnalysisRequest) -> str:
    job_id = uuid.uuid4().hex[:10]
    JOBS[job_id] = {"id": job_id, "status": "queued", "events": [], "started": time.time(), "result": None,
                    "queries": [q.name for q in req.queries], "tools": req.tools}
    EXEC.submit(_run_job, job_id, req)
    return job_id


@app.get("/")
def index():
    return FileResponse(STATIC / "index.html")


# The documentation site's "Launch web app" page asks this endpoint whether PathwayAgent is running
# on the visitor's machine. Only the health check is readable from other origins; it reveals the
# version and nothing else. The private-network header answers the browser's preflight for a public
# page contacting localhost.
_HEALTH_CORS = {"Access-Control-Allow-Origin": "*", "Access-Control-Allow-Methods": "GET, OPTIONS",
                "Access-Control-Allow-Headers": "*", "Access-Control-Allow-Private-Network": "true"}


@app.get("/api/v1/health")
def health():
    return JSONResponse({"ok": True, "version": __version__}, headers=_HEALTH_CORS)


@app.options("/api/v1/health")
def health_preflight():
    return Response(status_code=204, headers=_HEALTH_CORS)


@app.get("/api/v1/options")
def options():
    return {"organisms": {k: v["label"] for k, v in ORGANISMS.items()}, "databases": DATABASES,
            "default_databases": DEFAULT_DATABASES, "families": FAMILIES, "presets": PRESETS,
            "tools": {n: {"family": c.family, "method": c.method, "web": c.web, "databases": list(c.databases),
                          "needs_ranked": c.needs_ranked or c.needs_scores} for n, c in REGISTRY.items()},
            "planned": PLANNED, "defaults": DEFAULTS}


@app.post("/api/v1/jobs")
def create_job(job: JobIn):
    return {"job_id": _submit(_build_request(job))}


@app.post("/api/v1/jobs/upload")
async def create_job_upload(file: UploadFile = File(...), params: str = Form("{}")):
    p = JobIn(**json.loads(params or "{}"))
    p.text = (await file.read()).decode("utf-8-sig", errors="replace")
    p.name = Path(file.filename or "upload").stem
    return {"job_id": _submit(_build_request(p))}


@app.get("/api/v1/jobs/{job_id}")
def get_job(job_id: str):
    job = JOBS.get(job_id)
    if not job:
        raise HTTPException(404, "unknown job")
    return JSONResponse({k: v for k, v in job.items() if k != "started"} |
                        {"elapsed_s": round(time.time() - job["started"], 1)})


@app.get("/api/v1/jobs/{job_id}/events")
async def job_events(job_id: str):
    if job_id not in JOBS:
        raise HTTPException(404, "unknown job")

    async def stream():
        i = 0
        while True:
            job = JOBS[job_id]
            evs = job["events"]
            while i < len(evs):
                yield f"data: {json.dumps(evs[i], default=str)}\n\n"
                i += 1
            if job["status"] in ("done", "error") and i >= len(job["events"]):
                yield f"data: {json.dumps({'type': 'finished', 'status': job['status']})}\n\n"
                return
            await asyncio.sleep(0.25)

    return StreamingResponse(stream(), media_type="text/event-stream")


def warm_up():
    """Load libraries and start the R worker in the background so the first job is fast."""
    def _w():
        try:
            from ..core.idmap import get_mapper
            from ..core.library import get_library
            for org in ORGANISMS:
                try:
                    get_library(org)
                    get_mapper(org)
                except FileNotFoundError:
                    pass
            from ..adapters.rtools import get_worker, r_available
            if r_available():
                get_worker().start()
        except Exception as exc:
            print(f"[PathwayAgent] warm-up warning: {exc}")
    threading.Thread(target=_w, daemon=True).start()


def serve(port: int = 8765, open_browser: bool = True, host: str = "127.0.0.1"):
    import uvicorn

    warm_up()
    url = f"http://{host}:{port}"
    print(f"[PathwayAgent] serving on {url}  (Ctrl+C to stop)")
    if open_browser:
        threading.Timer(1.5, lambda: webbrowser.open(url)).start()
    uvicorn.run(app, host=host, port=port, log_level="warning")
