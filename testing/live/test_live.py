"""
Live tests against external services + end-to-end runs on the Example gene lists.
Regression tests for the bugs found in the CS_pathway_agent prototype (design §2).
Run with:  pytest -m live
"""
import json
import time

import pytest

from pathwayagent.adapters import REGISTRY
from pathwayagent.core.inputs import parse_file

from conftest import make_ctx

pytestmark = pytest.mark.live


@pytest.fixture(scope="module")
def module4(examples):
    q, = parse_file(examples / "genels_ex_module4.txt")
    req, ctx = make_ctx([q])
    return q, ctx


def test_gprofiler_hits_are_gene_symbols(module4):          # bug B2
    q, ctx = module4
    r = REGISTRY["gprofiler"]().run(q, ctx)
    assert r.status == "ok", r.message
    genes = set(q.genes)
    hits = {g for h in r.records["hit_genes"] for g in h.split(";") if g}
    assert hits and hits <= genes
    assert not hits & {"IEA", "ISS", "IDA", "IBA", "TAS"}
    # overlap equals number of hit genes
    ok = r.records["overlap"] == r.records["hit_genes"].str.split(";").apply(len)
    assert ok.mean() > 0.95


def test_gprofiler_returns_non_significant_terms(module4):  # bug B8
    q, ctx = module4
    r = REGISTRY["gprofiler"]().run(q, ctx)
    assert (r.records["padj"] >= 0.05).sum() > 100


def test_reactome_returns_found_genes(module4):             # bug B5
    q, ctx = module4
    r = REGISTRY["reactome"]().run(q, ctx)
    assert r.status == "ok", r.message
    assert (r.records["hit_genes"].str.len() > 0).mean() > 0.9


def test_enrichr_terms_map_to_library(module4):             # bug B9
    q, ctx = module4
    r = REGISTRY["enrichr"]().run(q, ctx)
    assert r.status == "ok", r.message
    rec = r.records
    lib_keys = set(ctx.library.sets["term_key"])
    for db in ("GO_BP", "HALLMARK", "REACTOME", "KEGG"):
        sub = rec[rec["database"] == db]
        frac = sub["term_key"].isin(lib_keys).mean()
        assert frac > 0.6, (db, frac)


def test_string_adapter(module4):
    q, ctx = module4
    r = REGISTRY["string"]().run(q, ctx)
    assert r.status == "ok", r.message
    assert r.records["term_key"].str.startswith("GO:").any()


@pytest.mark.slow
@pytest.mark.parametrize("fname", ["genels_ex_module123.txt", "genels_ex_module4.txt",
                                   "genels_ex_module4_wStat.csv", "genels_ex_module6.csv"])
def test_examples_end_to_end(examples, tmp_path, fname):
    from pathwayagent.config import PRESETS
    from pathwayagent.core.schema import AnalysisRequest
    from pathwayagent.pipeline import run_analysis
    qs = parse_file(examples / fname)
    t0 = time.time()
    run = run_analysis(AnalysisRequest(queries=qs, organism="auto", tools=PRESETS["standard"]), outdir=tmp_path)
    assert time.time() - t0 < 120
    assert (tmp_path / "report.html").exists() and (tmp_path / "summary.txt").exists()
    man = json.loads((tmp_path / "manifest.json").read_text())
    assert man["organism"] == "mmusculus"
    for qr in run.queries:
        st = {r.tool: r.status for r in qr.tool_results}
        assert st["local_ora"] == "ok"
        assert sum(s == "ok" for s in st.values()) >= 4, st
        assert not any(s == "error" for s in st.values()), [(r.tool, r.message) for r in qr.tool_results]
        assert len(qr.consensus) and (qr.consensus["tier"] > 0).any()


def test_server_job_roundtrip(examples):
    from fastapi.testclient import TestClient
    from pathwayagent.server.app import app
    c = TestClient(app)
    assert c.get("/api/v1/health").json()["ok"]
    opts = c.get("/api/v1/options").json()
    assert "local_ora" in opts["tools"]
    text = (examples / "genels_ex_module4.txt").read_text()
    jid = c.post("/api/v1/jobs", json={"text": text, "name": "m4", "tools": ["local_ora", "enrichr"],
                                       "databases": ["HALLMARK", "GO_BP"]}).json()["job_id"]
    for _ in range(240):
        j = c.get(f"/api/v1/jobs/{jid}").json()
        if j["status"] in ("done", "error"):
            break
        time.sleep(0.5)
    assert j["status"] == "done", j.get("error")
    q = j["result"]["queries"][0]
    assert q["consensus"] and q["figures"].get("consensus_bar")
    assert c.get(q["figures"]["consensus_bar"]).status_code == 200
    assert c.get(j["result"]["report"]).status_code == 200
