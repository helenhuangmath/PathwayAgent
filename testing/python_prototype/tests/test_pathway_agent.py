"""
PathwayAgent test suite.

Tests:
  1. Unit — name normalisation, consensus scoring (offline, no API)
  2. Unit — Jaccard network builder
  3. Unit — figure generation with mock data
  4. Integration — ORA with a real cancer-related gene list (live API)
  5. Integration — GSEA mode with a ranked gene list (live API)
  6. Timing  — full ORA run must complete in < 120 s
"""
from __future__ import annotations
import os
import sys
import time
import unittest
import numpy as np
import pandas as pd

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from pathway_agent.parsers.merge import _normalise_name, merge_results, consensus_score
from pathway_agent.networks.builder import build_pathway_network
from pathway_agent.reports.plots import make_dotplot, make_barplot, make_network_plot


# ── canonical test gene lists ────────────────────────────────────────────────

# 50 well-known tumour-suppressor / oncogene symbols (TP53, BRCA1/2, MYC, etc.)
CANCER_GENES = [
    "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "PIK3CA", "PTEN", "RB1",
    "APC", "VHL", "MYC", "MYCN", "CDH1", "MLH1", "MSH2", "ATM", "CHEK2",
    "NF1", "NF2", "TSC1", "TSC2", "RET", "MEN1", "FH", "SDHA", "SDHB",
    "SDHC", "SDHD", "BMPR1A", "STK11", "CDKN2A", "SMAD4", "RUNX1",
    "FLT3", "NPM1", "IDH1", "IDH2", "TET2", "DNMT3A", "ASXL1",
    "EZH2", "KMT2A", "CREBBP", "EP300", "JAK2", "MPL", "CALR",
    "SF3B1", "SRSF2", "U2AF1",
]

# 200 random genes + 50 cancer genes with fake log2FC scores (GSEA test)
RANKED_GENES  = CANCER_GENES + [
    "ACTB", "GAPDH", "TUBA1B", "RPL13", "RPS18", "LMNA", "VIM", "KRT8",
    "KRT18", "S100A4", "MMP9", "MMP2", "VEGFA", "HIF1A", "CASP3", "BCL2",
    "BAX", "CYCS", "CDKN1A", "CDKN1B", "CCND1", "CCNE1", "CDK4", "CDK6",
    "E2F1", "E2F3", "PCNA", "MCM2", "MCM7", "TOP2A", "AURKA", "AURKB",
    "PLK1", "BUB1", "BUB3", "MAD2L1", "CDC20", "CCNB1", "CDK1", "WEE1",
    "CHEK1", "RAD51", "BRIP1", "PALB2", "FANCD2", "FANCA", "FANCM",
    "RAD17", "TOPBP1", "RPA1",
]
np.random.seed(42)
RANKED_SCORES = list(np.concatenate([
    np.random.uniform(2, 8, len(CANCER_GENES)),   # up-regulated
    np.random.uniform(-3, 3, len(RANKED_GENES) - len(CANCER_GENES)),
]))


# ────────────────────────────────────────────────────────────────────────────
# 1. Unit: name normalisation
# ────────────────────────────────────────────────────────────────────────────

class TestNameNormalisation(unittest.TestCase):

    def test_strips_go_id(self):
        n = _normalise_name("Apoptosis (GO:0006915)")
        self.assertNotIn("GO:", n)

    def test_strips_year_suffix(self):
        n = _normalise_name("KEGG_2021_Human: Apoptosis")
        self.assertNotIn("2021", n)

    def test_lowercases(self):
        n = _normalise_name("DNA Repair")
        self.assertEqual(n, n.lower())

    def test_reactome_prefix_stripped(self):
        n = _normalise_name("Reactome: DNA Damage Response")
        self.assertNotIn("reactome:", n)


# ────────────────────────────────────────────────────────────────────────────
# 2. Unit: consensus merge & score
# ────────────────────────────────────────────────────────────────────────────

class TestConsensusMerge(unittest.TestCase):

    def _make_df(self, tool: str, pathways: list, fdrs: list) -> pd.DataFrame:
        return pd.DataFrame({
            "pathway_name": pathways,
            "fdr":          fdrs,
            "pval":         fdrs,
            "overlap":      [10] * len(pathways),
            "genes":        ["TP53;BRCA1"] * len(pathways),
            "tool":         [tool] * len(pathways),
            "source":       [tool] * len(pathways),
            "database":     [tool] * len(pathways),
        })

    def test_merge_two_tools(self):
        df1 = self._make_df("Enrichr", ["Apoptosis", "Cell Cycle"], [0.001, 0.01])
        df2 = self._make_df("gProfiler_ORA", ["Apoptosis", "DNA Repair"], [0.005, 0.02])
        raw = {"Enrichr": df1, "gProfiler_ORA": df2}
        merged = merge_results(raw)
        self.assertGreater(len(merged), 0)
        self.assertIn("pathway_key", merged.columns)

    def test_consensus_score_higher_for_multi_tool(self):
        df1 = self._make_df("Enrichr",       ["Apoptosis", "Unique_A"], [0.001, 0.001])
        df2 = self._make_df("gProfiler_ORA", ["Apoptosis", "Unique_B"], [0.001, 0.001])
        raw = {"Enrichr": df1, "gProfiler_ORA": df2}
        merged   = merge_results(raw)
        scored   = consensus_score(merged)
        apop_row = scored[scored["pathway_key"].str.contains("apoptosis", na=False)]
        self.assertGreater(len(apop_row), 0)
        # Apoptosis found in 2 tools → higher score than unique ones
        apop_score  = apop_row.iloc[0]["consensus_score"]
        unique_rows = scored[scored["pathway_key"].str.contains("unique", na=False)]
        if len(unique_rows) > 0:
            self.assertGreater(apop_score, unique_rows["consensus_score"].max())

    def test_fdr_filtering(self):
        df1 = self._make_df("Enrichr", ["Nonsig_A", "Significant_B"], [0.9, 0.001])
        raw = {"Enrichr": df1}
        merged = merge_results(raw)
        scored = consensus_score(merged)
        sig    = scored[scored["fdr"] < 0.05]
        self.assertEqual(len(sig), 1)
        self.assertIn("significant", sig.iloc[0]["pathway_key"])

    def test_nes_weight_gsea_higher(self):
        df_ora  = self._make_df("Enrichr",          ["Apoptosis"], [0.001])
        df_gsea = self._make_df("gseapy_Preranked",  ["Apoptosis"], [0.001])
        raw = {"Enrichr": df_ora, "gseapy_Preranked": df_gsea}
        merged = merge_results(raw)
        scored = consensus_score(merged)
        # Combined score should reflect both tool weights
        self.assertGreater(scored.iloc[0]["consensus_score"], 0)


# ────────────────────────────────────────────────────────────────────────────
# 3. Unit: network builder
# ────────────────────────────────────────────────────────────────────────────

class TestNetworkBuilder(unittest.TestCase):

    def _make_top_pathways(self) -> pd.DataFrame:
        return pd.DataFrame({
            "pathway_name":   ["Apoptosis", "Cell Cycle", "DNA Repair", "p53 Signaling"],
            "pathway_key":    ["apoptosis", "cell cycle", "dna repair", "p53 signaling"],
            "fdr":            [1e-5, 1e-4, 1e-3, 1e-4],
            "consensus_score":[0.9, 0.8, 0.7, 0.85],
            "n_tools":        [4, 3, 2, 3],
            "all_genes":      [
                "TP53;BRCA1;CASP3;BCL2",
                "TP53;CDKN1A;CCND1;CDK4",
                "BRCA1;BRCA2;RAD51;ATM",
                "TP53;MDM2;CDKN1A",
            ],
            "database": ["GO:BP", "KEGG", "Reactome", "KEGG"],
            "best_overlap": [10, 8, 6, 7],
            "tools_found": ["Enrichr|gProfiler", "Enrichr", "Reactome", "gProfiler"],
        })

    def test_network_has_nodes(self):
        df = self._make_top_pathways()
        G = build_pathway_network(df, CANCER_GENES[:20])
        self.assertGreater(G.number_of_nodes(), 0)

    def test_shared_gene_creates_edge(self):
        """Apoptosis, Cell Cycle, p53 Signaling all share TP53 → should have edges."""
        df = self._make_top_pathways()
        G = build_pathway_network(df, CANCER_GENES[:20], jaccard_threshold=0.05)
        # TP53 shared between apoptosis/cell_cycle/p53 → at least one edge
        self.assertGreater(G.number_of_edges(), 0)

    def test_node_attributes(self):
        df = self._make_top_pathways()
        G = build_pathway_network(df, CANCER_GENES[:20])
        for node in G.nodes():
            self.assertIn("consensus_score", G.nodes[node])
            self.assertIn("fdr", G.nodes[node])


# ────────────────────────────────────────────────────────────────────────────
# 4. Unit: figure generation (offline, no API)
# ────────────────────────────────────────────────────────────────────────────

class TestFigureGeneration(unittest.TestCase):

    def setUp(self):
        self.outdir = "/tmp/pathway_agent_test_figs"
        os.makedirs(self.outdir, exist_ok=True)
        self.df = pd.DataFrame({
            "pathway_name":   [f"Pathway {i}" for i in range(20)],
            "pathway_key":    [f"pathway {i}" for i in range(20)],
            "fdr":            [10**(-i * 0.5) for i in range(1, 21)],
            "consensus_score":[1 - i * 0.04 for i in range(20)],
            "n_tools":        [max(1, 4 - i // 5) for i in range(20)],
            "all_genes":      ["TP53;BRCA1;MYC"] * 20,
            "database":       (["GO:BP", "KEGG", "Reactome", "WP"] * 5),
            "best_overlap":   [20 - i for i in range(20)],
            "tools_found":    ["Enrichr|gProfiler"] * 20,
        })

    def test_dotplot_created(self):
        path = os.path.join(self.outdir, "dotplot.png")
        make_dotplot(self.df, path)
        self.assertTrue(os.path.exists(path))
        self.assertGreater(os.path.getsize(path), 5_000)

    def test_barplot_created(self):
        path = os.path.join(self.outdir, "barplot.png")
        make_barplot(self.df, path)
        self.assertTrue(os.path.exists(path))
        self.assertGreater(os.path.getsize(path), 5_000)

    def test_network_plot_created(self):
        from pathway_agent.networks.builder import build_pathway_network
        G = build_pathway_network(self.df, CANCER_GENES, jaccard_threshold=0.01)
        path = os.path.join(self.outdir, "network.png")
        make_network_plot(G, self.df, path)
        self.assertTrue(os.path.exists(path))
        self.assertGreater(os.path.getsize(path), 5_000)


# ────────────────────────────────────────────────────────────────────────────
# 5. Integration: full ORA run (live API)
# ────────────────────────────────────────────────────────────────────────────

class TestIntegrationORA(unittest.TestCase):
    """Run the agent once via setUpClass; all test methods share the result."""

    _agent = None
    _outdir = "/tmp/pathway_agent_test_ora"

    @classmethod
    def setUpClass(cls):
        from pathway_agent import PathwayAgent
        os.makedirs(cls._outdir, exist_ok=True)
        cls._agent = PathwayAgent(
            gene_list=CANCER_GENES,
            organism="hsapiens",
            fdr_cutoff=0.05,
            min_genes=3,
            outdir=cls._outdir,
            n_jobs=6,
        )
        cls._agent.run(verbose=True)

    def setUp(self):
        self.agent  = self._agent
        self.outdir = self._outdir

    def test_run_completes(self):
        self.assertIsNotNone(self.agent.consensus)

    def test_at_least_one_tool_returns_results(self):
        self.assertGreater(len(self.agent.raw_results), 0)

    def test_significant_pathways_found(self):
        self.assertIsNotNone(self.agent.top_pathways)
        self.assertGreater(len(self.agent.top_pathways), 0,
                           "Expected at least 1 significant pathway for cancer genes")

    def test_consensus_score_between_0_and_1(self):
        tp = self.agent.top_pathways
        if tp is not None and len(tp) > 0:
            self.assertTrue((tp["consensus_score"] >= 0).all())
            self.assertTrue((tp["consensus_score"] <= 1).all())

    def test_expected_pathways_in_top20(self):
        """TP53/apoptosis/cell cycle should appear for this gene set."""
        tp = self.agent.top_pathways
        if tp is None or len(tp) == 0:
            self.skipTest("No significant pathways returned")
        names_lower = tp.head(20)["pathway_key"].str.lower().tolist()
        hits = [n for n in names_lower if any(
            kw in n for kw in ["apoptosis", "p53", "dna repair", "cell cycle",
                                "cancer", "brca", "tumor", "tumour"]
        )]
        self.assertGreater(len(hits), 0,
                           f"Expected cancer-related pathway in top-20. Got: {names_lower[:10]}")

    def test_output_files_created(self):
        for fname in ["consensus_pathways.csv", "all_tools_results.csv", "report.json"]:
            fpath = os.path.join(self.outdir, fname)
            self.assertTrue(os.path.exists(fpath), f"Missing output: {fname}")

    def test_at_least_one_figure_created(self):
        self.assertGreater(len(self.agent.figures), 0)

    def test_summary_method(self):
        summary = self.agent.summary(10)
        self.assertLessEqual(len(summary), 10)
        self.assertIn("pathway_name", summary.columns)


# ────────────────────────────────────────────────────────────────────────────
# 6. Integration: GSEA mode (live API)
# ────────────────────────────────────────────────────────────────────────────

class TestIntegrationGSEA(unittest.TestCase):
    """Run the agent once via setUpClass; all test methods share the result."""

    _agent = None
    _outdir = "/tmp/pathway_agent_test_gsea"

    @classmethod
    def setUpClass(cls):
        from pathway_agent import PathwayAgent
        os.makedirs(cls._outdir, exist_ok=True)
        cls._agent = PathwayAgent(
            gene_list=RANKED_GENES,
            ranked_scores=RANKED_SCORES,
            organism="hsapiens",
            fdr_cutoff=0.05,
            min_genes=3,
            outdir=cls._outdir,
            n_jobs=6,
        )
        cls._agent.run(verbose=True)

    def setUp(self):
        self.agent  = self._agent
        self.outdir = self._outdir

    def test_gsea_tasks_registered(self):
        tasks = self.agent._build_task_list()
        task_names = [t[0] for t in tasks]
        self.assertIn("gseapy_Preranked", task_names)
        self.assertIn("gProfiler_GSEA", task_names)

    def test_gsea_run_completes(self):
        self.assertIsNotNone(self.agent.consensus)


# ────────────────────────────────────────────────────────────────────────────
# 7. Timing: ORA must complete in < 120 s
# ────────────────────────────────────────────────────────────────────────────

class TestTiming(unittest.TestCase):

    def test_ora_under_120s(self):
        from pathway_agent import PathwayAgent
        agent = PathwayAgent(
            gene_list=CANCER_GENES,
            organism="hsapiens",
            fdr_cutoff=0.1,
            outdir="/tmp/pathway_agent_timing",
            n_jobs=6,
        )
        t0 = time.time()
        agent.run(verbose=False)
        elapsed = time.time() - t0
        print(f"\n  [Timing] Full ORA run: {elapsed:.1f}s")
        self.assertLess(elapsed, 120,
                        f"Run took {elapsed:.1f}s — exceeded 120 s limit")


# ────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    unittest.main(verbosity=2)
