"""Offline unit tests: input parsing, ID mapping, statistics, consensus logic, themes."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats as sps

from pathwayagent.core.inputs import parse_file, parse_text
from pathwayagent.core import stats
from pathwayagent.consensus.score import consensus, rra_pvalue
from pathwayagent.consensus.themes import kappa_matrix


# ───────────────────────────── inputs ─────────────────────────────

def test_parse_multicolumn_modules(examples):
    qs = parse_file(examples / "genels_ex_module123.txt")
    assert [q.name for q in qs] == ["Module1", "Module2", "Module3"]
    assert all(len(q.genes) > 20 for q in qs)
    assert all(not q.is_ranked for q in qs)
    assert all("" not in q.genes for q in qs)


def test_parse_multicolumn_stem_header():
    qs = parse_text("TEX_G1,TEX_G2,TEX_G5G6\r\nCcl5,Rnf213,Chsy1\r\nPrex1,,Nsmaf\r\n")
    assert [q.name for q in qs] == ["TEX_G1", "TEX_G2", "TEX_G5G6"]
    assert [q.genes for q in qs] == [["Ccl5", "Prex1"], ["Rnf213"], ["Chsy1", "Nsmaf"]]
    qs = parse_text("Gimap3,Gimap7\nGimap4,Gimap8\n")            # gene family, no header
    assert [q.genes for q in qs] == [["Gimap3", "Gimap4"], ["Gimap7", "Gimap8"]]


def test_parse_single_list_crlf(examples):
    q, = parse_file(examples / "genels_ex_module4.txt")
    assert len(q.genes) == 248 and not q.is_ranked
    assert not any(g.endswith("\r") for g in q.genes)


def test_parse_ranked_with_header(examples):
    q, = parse_file(examples / "genels_ex_module4_wStat.csv")
    assert q.is_ranked and len(q.scores) == 248
    assert q.ranked_genes()[0] == max(q.scores, key=q.scores.get)


def test_parse_module6(examples):
    q, = parse_file(examples / "genels_ex_module6.csv")
    assert len(q.genes) >= 390


def test_parse_de_table():
    txt = "gene,log2FoldChange,pvalue\nTox,2.0,1e-5\nTcf7,-1.5,1e-3\nPdcd1,0.1,0.5\n"
    q, = parse_text(txt)
    assert q.is_ranked
    assert q.scores["Tox"] > q.scores["Pdcd1"] > q.scores["Tcf7"]


# ───────────────────────────── Metascape import ───────────────────

def test_metascape_package(tmp_path):
    from conftest import make_ctx
    from pathwayagent.adapters.metascape import Metascape
    from pathwayagent.core.schema import QueryInput
    go = tmp_path / "pkg" / "Enrichment_GO"
    go.mkdir(parents=True)
    head = ("Category,CategoryID,GO,Description,PARENT_GO,LogP,Enrichment,Z-score,#TotalGeneInLibrary,#GeneInGO,"
            "#GeneInHitList,#GeneInGOAndHitList,%InGO,STDV %InGO,GeneID,Hits,Log(q-value),GeneList\n")
    (go / "GO_AllLists.csv").write_text(
        head + "GO Biological Processes,19,GO:0007059,chromosome segregation,,-6,9,8,20000,300,3,3,100,0,1|2|3,"
               "Cdk1|Plk1|Ccnb1,-3,listA\n"
               "KEGG Pathway,24,mmu04110,Cell cycle,,-4,8,7,20000,120,3,3,100,0,1|2|3,Cdk1|Plk1|Ccnb1,-0.5,listA\n"
               "CORUM,33,CORUM:1,Some complex,,-4,8,7,20000,12,3,3,100,0,1|2|3,Cdk1|Plk1|Ccnb1,-2,listA\n"
               "GO Biological Processes,19,GO:0006955,immune response,,-5,4,6,20000,900,3,3,100,0,4|5|6,"
               "Ifng|Tox|Batf,-2,listB\n")
    (go / "_FINAL_GO.csv").write_text("GO,GROUP_ID\nGO:0007059,1\nmmu04110,1\nGO:0006955,2\n")
    qa = QueryInput("listA", ["Cdk1", "Plk1", "Ccnb1"], raw_genes=["Cdk1", "Plk1", "Ccnb1"])
    _, ctx = make_ctx([qa], tools=("metascape",), metascape={"*": str(tmp_path / "pkg")})
    res = Metascape().run(qa, ctx)
    d = res.records.set_index("term_key")
    assert res.status == "ok" and set(d.index) == {"GO:0007059", "KEGG:04110"}      # own list, known databases
    assert d.loc["GO:0007059", "padj"] == pytest.approx(1e-3) and d.loc["GO:0007059", "pval"] == pytest.approx(1e-6)
    assert d.loc["KEGG:04110", "hit_genes"] == "Ccnb1;Cdk1;Plk1"
    assert res.databases == ["GO_BP", "KEGG"] and res.extra["categories_not_used"] == ["CORUM"]
    assert {c["GROUP_ID"] for c in res.extra["clusters"]} == {1}
    _, ctx = make_ctx([qa], tools=("metascape",))
    res = Metascape().run(qa, ctx)
    assert res.status == "skipped" and "--metascape listA=PATH" in res.message


# ───────────────────────────── ID mapping ─────────────────────────

def test_idmap(mapper):
    assert mapper.map_one("Tox") == ("Tox", "exact")
    assert mapper.map_one("TOX") == ("Tox", "case")
    sym, how = mapper.map_one("21859")          # Entrez id of Timp3? -> any mapped symbol
    assert how in ("entrez", "unmapped")
    s, how = mapper.map_one("ENSMUSG00000041872")
    assert how in ("ensembl", "unmapped")
    assert mapper.map_one("NotAGene123")[1] == "unmapped"


def test_guess_organism():
    from pathwayagent.core.idmap import guess_organism
    assert guess_organism(["Tox", "Tcf7", "Pdcd1"]) == "mmusculus"
    assert guess_organism(["TOX", "TCF7", "PDCD1"]) == "hsapiens"


# ───────────────────────────── statistics ─────────────────────────

def test_bh_matches_reference():
    p = np.array([0.01, 0.04, 0.03, 0.2, 0.001])
    ref = np.array([0.025, 0.05, 0.05, 0.2, 0.005])
    assert np.allclose(stats.bh(p), ref)


def test_ora_equals_fisher_exact(lib):
    rows = lib.db_rows(["HALLMARK"])
    uni = lib.default_universe(["HALLMARK"])
    rng = np.random.default_rng(0)
    target = rows[0]
    members = lib.set_gene_idx(target)
    query = np.unique(np.concatenate([members[:20], rng.choice(uni, 180, replace=False)]))
    df = stats.ora(lib.matrix, rows, uni, query, lib.genes, 5, 5000)
    r = df[df["row"] == target].iloc[0]
    k, K, n, N = r.overlap, r.term_size, r.query_size, r.universe_size
    _, p = sps.fisher_exact([[k, n - k], [K - k, N - K - n + k]], alternative="greater")
    assert np.isclose(r.pval, p, rtol=1e-6)
    assert set(r.hit_genes) <= set(lib.genes[query])      # hits are query genes only


def test_mhg_prefers_top_ranked(lib):
    rows = lib.db_rows(["HALLMARK"])
    uni = lib.default_universe(["HALLMARK"])
    members = lib.set_gene_idx(rows[3])
    others = np.setdiff1d(uni, members)[:300]
    ranked_top = np.concatenate([members[:25], others])
    ranked_bottom = np.concatenate([others, members[:25]])
    p_top = stats.mhg(lib.matrix, rows, uni, ranked_top, lib.genes, 5, 5000).set_index("row").loc[rows[3], "pval"]
    p_bot = stats.mhg(lib.matrix, rows, uni, ranked_bottom, lib.genes, 5, 5000).set_index("row").loc[rows[3], "pval"]
    assert p_top < 1e-10 and p_top < p_bot


def test_rank_test_direction(lib):
    rows = lib.db_rows(["HALLMARK"])
    uni = lib.default_universe(["HALLMARK"])
    rng = np.random.default_rng(1)
    scored = rng.choice(uni, 3000, replace=False)
    scores = rng.normal(size=len(scored))
    members = set(lib.set_gene_idx(rows[5]).tolist())
    scores[[i for i, g in enumerate(scored) if g in members]] += 3
    df = stats.rank_test(lib.matrix, rows, scored, scores, lib.genes, 5, 5000).set_index("row")
    assert df.loc[rows[5], "pval"] < 1e-6 and df.loc[rows[5], "direction"] == 1


def test_kappa_and_rra():
    X = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1]])
    K = kappa_matrix(X)
    assert np.allclose(K, K.T) and np.isclose(K[0, 1], 1.0) and K[0, 2] < 0
    assert rra_pvalue(np.array([0.001, 0.002, 0.003])) < 1e-5
    assert rra_pvalue(np.array([0.5, 0.6, 0.9])) > 0.5


# ───────────────────────────── consensus logic ────────────────────

def _rec(tool, family, key, padj, db="GO_BP", genes="A;B"):
    return {"query": "q", "tool": tool, "method": "m", "family": family, "mode": "x", "database": db,
            "term_id": key, "term_key": key, "term_name": key, "term_size": 50, "query_size": 100,
            "universe_size": 10000, "overlap": 5, "pval": padj, "padj": padj, "padj_method": "BH",
            "effect": 3.0, "effect_type": "odds_ratio", "direction": 0, "hit_genes": genes,
            "rank": 1, "n_tested": 100}


def test_eligibility_uncovered_tool_does_not_penalise(lib):
    # 'reactome' only covers REACTOME, so it must not count against a GO term
    recs = pd.DataFrame([_rec("t1", "F1_ORA", "GO:X", 1e-6), _rec("t2", "F1_ORA", "GO:X", 1e-5)])
    info = {"t1": {"family": "F1_ORA", "databases": ["GO_BP"]},
            "t2": {"family": "F1_ORA", "databases": ["GO_BP"]},
            "reactome": {"family": "F1_ORA", "databases": ["REACTOME"]}}
    c = consensus(recs, info, lib).set_index("term_key")
    assert c.loc["GO:X", "tool_support"] == "2/2"
    assert c.loc["GO:X", "tier_label"] == "ROBUST"


def test_family_balancing(lib):
    # 6 ORA tools agree on A; 1 GSEA tool disagrees. B is supported by both families.
    recs, info = [], {}
    for i in range(6):
        info[f"ora{i}"] = {"family": "F1_ORA", "databases": ["GO_BP"]}
        recs.append(_rec(f"ora{i}", "F1_ORA", "GO:A", 1e-6))
        recs.append(_rec(f"ora{i}", "F1_ORA", "GO:B", 1e-4))
    info["gsea"] = {"family": "F2_GSEA", "databases": ["GO_BP"]}
    recs.append(_rec("gsea", "F2_GSEA", "GO:B", 1e-4))
    recs.append(_rec("gsea", "F2_GSEA", "GO:A", 0.9))
    c = consensus(pd.DataFrame(recs), info, lib).set_index("term_key")
    assert c.loc["GO:B", "n_families_sig"] == 2 and c.loc["GO:A", "n_families_sig"] == 1
    assert c.loc["GO:B", "tier"] == 1          # both families agree
    assert c.loc["GO:A", "tier"] != 1          # 6 correlated ORA votes are not enough
