# Small worked example: mouse interferon genes

This example uses 20 illustrative mouse genes selected for interferon biology. It is a demonstration,
not an experimental dataset or an unbiased scientific benchmark. The script is included in the installed
package at `pathwayagent.examples.interferon`.

## Install and build the library

Download the wheel from the [v0.1.0 GitHub release](https://github.com/helenhuangmath/PathwayAgent/releases/tag/v0.1.0), then:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install ./pathwayagent-0.1.0-py3-none-any.whl
export PATHWAYAGENT_LIBRARY="$PWD/library"
export PATHWAYAGENT_RESULTS="$PWD/results"
export MPLCONFIGDIR="$PWD/.mplconfig"
.venv/bin/python -m pathwayagent db build --organism mmusculus
```

Library construction downloads reference data. R is optional for this example. If you already have a
built mouse library, point `PATHWAYAGENT_LIBRARY` at its parent directory (the one containing `mmusculus/`)
and skip the build. Reference releases evolve; numerical results can change when rebuilding databases.

## Run

```bash
.venv/bin/python -m pathwayagent.examples.interferon --out results/interferon
```

The analysis runs only `local_ora` against Hallmark sets, with no web-service calls or R requirement.
The example uses the annotated Hallmark genes as background; real experiments should supply their
measured/expressed gene universe through the Python API or CLI.

## Inspect

- `results/interferon/report.html`: open in a browser for tables and figures.
- `results/interferon/interferon/tools/local_ora.tsv`: per-term test results.
- `results/interferon/interferon/consensus_terms.tsv`: significant exploratory terms.
- `results/interferon/interferon/results.xlsx`: workbook with results and tool status.
- `results/interferon/manifest.json`: version, parameters, database information, and timings.

Expect interferon response terms near the top, with low adjusted p-values. Every significant term
has SINGLE-TOOL evidence: this example demonstrates enrichment and reporting, not cross-tool consensus.
The script prints the top five terms and report path and raises an error if the adapter does not succeed.

With the existing MSigDB 2026.1.Mm library used for release verification, interferon γ and α response
were the top two terms, both with BH-adjusted p = 3.64 × 10⁻¹⁹. IL-6 JAK STAT3 signaling followed
(adjusted p = 8.71 × 10⁻⁴). Exact values and order depend on the reference library version.

The [Python API reference](PYTHON_API.md) shows how to supply a background, use additional tools,
inspect partial failures, and work directly with returned DataFrames.
