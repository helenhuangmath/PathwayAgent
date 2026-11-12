"""
Validation tests (design §15):
  1. planted signal is recovered and ranked top
  2. null calibration: random lists rarely give significant terms
  4. local ORA == clusterProfiler::enricher on the same gene sets (R)
  6. timing: the local "quick" preset answers in seconds
"""
import time

import numpy as np
import pytest

from pathwayagent.core import stats
from pathwayagent.core.schema import QueryInput
from pathwayagent.pipeline import run_analysis
from pathwayagent.core.schema import AnalysisRequest

from conftest import make_ctx


def _planted_query(lib, set_name="HALLMARK:HALLMARK_INTERFERON_ALPHA_RESPONSE", frac=0.6, n_noise=200, seed=0):
    rng = np.random.default_rng(seed)
    row = lib.set_pos[set_name]
    members = lib.set_genes(row)
    planted = list(rng.choice(members, int(frac * len(members)), replace=False))
    uni = lib.genes[lib.default_universe(["HALLMARK", "GO_BP", "KEGG", "REACTOME"])]
    noise = list(rng.choice(np.setdiff1d(uni, members), n_noise, replace=False))
    genes = planted + noise
    scores = {g: 10.0 - i * 0.01 for i, g in enumerate(genes)}   # planted genes ranked on top
    return QueryInput("planted", genes=genes, scores=scores, raw_genes=genes), lib.sets["term_key"].iat[row]


def test_planted_signal_recovered(lib):
    q, key = _planted_query(lib)
    req = AnalysisRequest(queries=[q], organism="mmusculus", tools=["local_ora", "local_mhg"],
                          databases=["HALLMARK", "GO_BP", "KEGG", "REACTOME"])
    run = run_analysis(req, make_outputs=False)
    c = run.queries[0].consensus
    top3 = c.head(3)["term_key"].tolist()
    assert key in top3, top3
    assert c.set_index("term_key").loc[key, "tier_label"] == "ROBUST"
    themes = run.queries[0].themes
    assert key in themes.iloc[0]["member_keys"].split("|")


@pytest.mark.slow
def test_null_calibration(lib):
    """Under the global null, BH controls the chance of >=1 false positive per database at ~alpha."""
    rng = np.random.default_rng(42)
    dbs = ["HALLMARK", "KEGG", "REACTOME", "GO_MF"]
    uni = lib.default_universe(dbs)
    hits = {d: 0 for d in dbs}
    n_runs = 60
    for _ in range(n_runs):
        q = rng.choice(uni, 200, replace=False)
        for d in dbs:
            df = stats.ora(lib.matrix, lib.db_rows([d]), uni, q, lib.genes, 10, 500)
            if len(df) and (stats.bh(df["pval"].values) < 0.05).any():
                hits[d] += 1
    rates = {d: h / n_runs for d, h in hits.items()}
    assert all(r <= 0.12 for r in rates.values()), rates


@pytest.mark.r
def test_local_ora_matches_clusterprofiler(examples):
    from pathwayagent.adapters import REGISTRY
    from pathwayagent.adapters.rtools import r_available
    from pathwayagent.core.inputs import parse_file
    if not r_available():
        pytest.skip("Rscript not available")
    q, = parse_file(examples / "genels_ex_module4.txt")
    req, ctx = make_ctx([q], databases=["HALLMARK", "KEGG", "GO_MF"])
    a = REGISTRY["local_ora"]().run(q, ctx).records.set_index("term_key")
    b = REGISTRY["clusterprofiler"]().run(q, ctx).records.set_index("term_key")
    common = a.index.intersection(b.index)
    assert len(common) > 100
    assert np.allclose(a.loc[common, "pval"], b.loc[common, "pval"], rtol=1e-5)
    assert np.allclose(a.loc[common, "padj"], b.loc[common, "padj"], rtol=1e-5)
    assert (a.loc[common, "overlap"] == b.loc[common, "overlap"]).all()


def test_quick_preset_is_fast(examples, lib):
    from pathwayagent.core.inputs import parse_file
    qs = parse_file(examples / "genels_ex_module123.txt")
    req = AnalysisRequest(queries=qs, organism="mmusculus", tools=["local_ora", "local_mhg"])
    run_analysis(req, make_outputs=False)          # warm caches (library load)
    t0 = time.time()
    run = run_analysis(AnalysisRequest(queries=parse_file(examples / "genels_ex_module123.txt"),
                                       organism="mmusculus", tools=["local_ora"]), make_outputs=False)
    dt = time.time() - t0
    assert dt < 5, f"local ORA on 3 lists took {dt:.1f}s"
    assert all(r.status == "ok" for qr in run.queries for r in qr.tool_results)
