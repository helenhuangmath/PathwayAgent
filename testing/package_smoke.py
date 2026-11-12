"""Exercise an installed distribution with synthetic data, never a downloaded reference library.

Run from outside the checkout with the installed wheel's Python interpreter.
"""
from __future__ import annotations

import importlib.metadata
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
from scipy import sparse


def main():
    with tempfile.TemporaryDirectory(prefix="pathwayagent-package-") as tmp:
        root = Path(tmp)
        os.environ["PATHWAYAGENT_LIBRARY"] = str(root / "library")
        os.environ["PATHWAYAGENT_RESULTS"] = str(root / "results")
        os.environ["MPLCONFIGDIR"] = str(root / "matplotlib")
        library = root / "library" / "mmusculus"
        library.mkdir(parents=True)
        genes = [f"Gene{i:03}" for i in range(200)]
        (library / "gene_index.txt").write_text("\n".join(genes))
        pd.DataFrame({
            "symbol": genes, "synonyms": "", "entrez": "", "ensembl": "", "ortholog": "",
        }).to_parquet(library / "genes.parquet", index=False)
        sets = []
        matrix = np.zeros((4, len(genes)), dtype=np.float32)
        for i in range(4):
            key = f"HALLMARK:SYNTHETIC_{i}"
            sets.append({
                "set_id": key, "term_key": key, "database": "HALLMARK",
                "name": f"Synthetic {i}", "name_norm": f"synthetic {i}", "url": "",
                "size": 50, "source": "Synthetic package fixture", "go_depth": 0,
            })
            matrix[i, i * 50:(i + 1) * 50] = 1
        pd.DataFrame(sets).to_parquet(library / "sets.parquet", index=False)
        sparse.save_npz(library / "membership.npz", sparse.csr_matrix(matrix))
        (library / "manifest.json").write_text(json.dumps({"built": "synthetic", "n_sets": {"HALLMARK": 4}}))

        import pathwayagent
        from pathwayagent.core.schema import AnalysisRequest, QueryInput
        from pathwayagent.pipeline import run_analysis
        from pathwayagent.examples.interferon import GENES

        assert pathwayagent.__version__ == importlib.metadata.version("pathwayagent")
        assert len(GENES) == 20
        query = QueryInput("synthetic", genes=genes[:20], raw_genes=genes[:20])
        request = AnalysisRequest(queries=[query], tools=["local_ora"], databases=["HALLMARK"], local_only=True)
        run = run_analysis(request, outdir=root / "results" / "smoke")
        result = run.queries[0]
        assert result.tool_results[0].status == "ok", result.tool_results[0].message
        top = result.consensus.iloc[0]
        assert top["term_key"] == "HALLMARK:SYNTHETIC_0"
        assert top["best_padj"] < 0.05
        assert top["tier_label"] == "SINGLE-TOOL"
        for name in ["report.html", "manifest.json", "methods.txt", "summary.txt", "synthetic/results.xlsx",
                     "synthetic/consensus_terms.tsv", "synthetic/figures/01_consensus_bar.pdf"]:
            assert (run.outdir / name).stat().st_size > 0, name
        from pathwayagent.server.app import app
        from fastapi.testclient import TestClient

        with TestClient(app) as client:
            assert client.get("/api/v1/health").status_code == 200
            assert client.get("/").status_code == 200
            assert client.get("/examples/examples.json").status_code == 200
        from importlib.resources import files
        assert files("pathwayagent").joinpath("rworker/worker.R").is_file()
        print(f"Package smoke passed: v{pathwayagent.__version__}, imported from {pathwayagent.__file__}")


if __name__ == "__main__":
    main()
