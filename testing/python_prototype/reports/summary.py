"""Write CSV and JSON summary reports."""
from __future__ import annotations
import json
import pandas as pd
from datetime import datetime


def write_summary_table(raw_results: dict, path: str):
    """Write all-tools results to a single CSV."""
    frames = []
    for tool, df in raw_results.items():
        if df is not None and not df.empty:
            frames.append(df)
    if not frames:
        return
    pd.concat(frames, ignore_index=True).to_csv(path, index=False)


def write_json_report(top_pathways: pd.DataFrame,
                      raw_results: dict,
                      timings: dict,
                      path: str):
    """Write a structured JSON report for downstream parsing."""
    report = {
        "generated_at": datetime.utcnow().isoformat() + "Z",
        "n_tools_run":  len(raw_results),
        "tools":        list(raw_results.keys()),
        "timings_sec":  {k: round(v, 2) for k, v in timings.items()},
        "n_significant": len(top_pathways),
        "top_pathways": [],
    }

    cols = ["pathway_name", "fdr", "consensus_score", "n_tools",
            "tools_found", "best_overlap", "all_genes", "database"]
    for _, row in top_pathways.head(100).iterrows():
        entry = {}
        for c in cols:
            val = row.get(c, None)
            if pd.isna(val) if not isinstance(val, str) else False:
                val = None
            elif isinstance(val, float):
                val = round(val, 6)
            entry[c] = val
        report["top_pathways"].append(entry)

    with open(path, "w") as f:
        json.dump(report, f, indent=2, default=str)
