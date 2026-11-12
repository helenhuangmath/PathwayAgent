"""Run a small, local mouse interferon example: python -m pathwayagent.examples.interferon."""
from __future__ import annotations

import argparse
from pathlib import Path

from pathwayagent.core.schema import AnalysisRequest, QueryInput
from pathwayagent.pipeline import run_analysis

# Illustrative genes selected for interferon biology, not an experimental dataset.
GENES = [
    "Stat1", "Stat2", "Irf7", "Irf9", "Isg15", "Ifit1", "Ifit2", "Ifit3",
    "Mx1", "Mx2", "Oas1a", "Oas2", "Rsad2", "Usp18", "Ddx58", "Ifih1",
    "Ifitm3", "Bst2", "Gbp2", "Cxcl10",
]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/interferon"))
    args = parser.parse_args(argv)
    request = AnalysisRequest(
        queries=[QueryInput("interferon", genes=list(GENES), raw_genes=list(GENES))],
        organism="mmusculus", databases=["HALLMARK"], tools=["local_ora"],
        local_only=True, title="Small mouse interferon example",
    )
    run = run_analysis(request, outdir=args.out)
    result = run.queries[0]
    failed = [r for r in result.tool_results if r.status != "ok"]
    if failed:
        raise RuntimeError("; ".join(f"{r.tool}: {r.status}: {r.message}" for r in failed))
    print(result.consensus[["term_name", "tier_label", "best_padj"]].head(5).to_string(index=False))
    print(f"Report: {(run.outdir / 'report.html').resolve()}")
    return run


if __name__ == "__main__":
    main()
