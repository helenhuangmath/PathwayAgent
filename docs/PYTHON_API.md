# Python API (v1.0.1)

This reference documents the supported Python API for PathwayAgent v1.0.1.
Python 3.11 or newer is required.

## Minimal analysis

Build a database before the first analysis (`python -m pathwayagent db build --organism mmusculus`).
Set `PATHWAYAGENT_LIBRARY` and `PATHWAYAGENT_RESULTS` **before importing** PathwayAgent when using
custom paths. For an installed wheel, explicitly choose writable paths; see [INSTALL.md](INSTALL.md).

```python
from pathwayagent.core.inputs import parse_text
from pathwayagent.core.schema import AnalysisRequest
from pathwayagent.pipeline import run_analysis

queries = parse_text("Stat1\nStat2\nIrf7\nIrf9\nIsg15\nIfit1\nIfit2\nMx1", name="interferon")
request = AnalysisRequest(
    queries=queries,
    organism="mmusculus",
    databases=["HALLMARK"],
    tools=["local_ora"],
    local_only=True,
)
run = run_analysis(request, outdir="results/python-example")
query = run.queries[0]
print(query.consensus[["term_name", "tier_label", "consensus_score", "best_padj"]].head())
print(run.outdir / "report.html")
```

A single tool produces exploratory SINGLE-TOOL evidence. Cross-tool agreement requires multiple
eligible tools; additional algorithms using the same data do not constitute independent biological validation.
See the [worked example](WORKED_EXAMPLE.md) for an executable version and expected results.

## Inputs

`parse_text(text: str, name: str = "query") -> list[QueryInput]` and
`parse_file(path) -> list[QueryInput]` detect gene lists, delimited ranked lists,
DE tables, and multiple gene columns. `parse_file` accepts text/CSV/TSV and Excel files.
DE scores use a statistic column when available, otherwise signed −log10(p) when fold-change
and p-value columns exist. Check parsed lists before analysis when using unusual headers.

Construct `QueryInput(name, genes, scores=None, raw_genes=[], mapping=None)` directly if needed:

```python
from pathwayagent.core.schema import QueryInput

scores = {"Stat1": 8.0, "Irf7": 6.0, "Tcf7": -2.0}
query = QueryInput("ranked", genes=list(scores), scores=scores, raw_genes=list(scores))
```

Higher scores mean more positive regulation. `is_ranked` reports whether scores exist;
`ranked_genes()` returns genes sorted by decreasing score. `mapping` is populated by analysis.
Use at least 1,000 mapped, background-covered genes for GSEA; partial rankings can use `local_mhg`.

## AnalysisRequest

| Field | Default | Meaning |
|---|---|---|
| `queries` | Required | List of `QueryInput` objects; use distinct names |
| `organism` | `"mmusculus"` | `"mmusculus"`, `"hsapiens"`, or `"auto"`; explicit selection is recommended |
| `universe` | `None` | Background gene symbols/IDs; supply measured/expressed genes when available |
| `databases` | All configured databases | List of database keys, such as `"HALLMARK"`, `"GO_BP"`, `"REACTOME"` |
| `tools` | Standard preset | List of adapter names; use `list(PRESETS["quick"])` for local tools |
| `min_size`, `max_size` | `10`, `500` | Gene-set size bounds |
| `alpha` | `0.05` | Adjusted p-value threshold |
| `deadline_s` | `60.0` | Time budget for the tool execution stage; excludes preparation and reporting |
| `local_only` | `False` | Skip adapters that send genes to external services |
| `min_ranked_for_gsea` | `1000` | Minimum background-covered ranking size for GSEA/rank-sum adapters |
| `profile` | `None` | Optional dictionary; `genes_of_interest` affects gene annotations |
| `metascape` | `None` | Mapping of query names (or `"*"`) to imported Metascape packages |
| `title` | `"PathwayAgent run"` | Report title |

Get presets from `pathwayagent.config.PRESETS`; list adapters with `python -m pathwayagent tools`.
Unranked lists skip adapters requiring scores. Missing R installations and tool errors appear in
per-tool status rather than necessarily failing the entire run.

## Execution and returned results

`run_analysis(request, outdir=None, on_event=None, make_outputs=True) -> RunResult`

- `outdir`: string or `Path`; defaults to the configured results directory plus a timestamp.
  Choose a fresh directory per run to avoid overwriting earlier files.
- `make_outputs=False`: compute in memory without writing reports or figures. Library loading
  and adapter caches can still access disk.
- `on_event`: optional callable receiving event dictionaries for stages, tool completion,
  query completion, and run completion. Callback exceptions are suppressed.

The request and its queries are modified during organism selection and ID mapping. Construct
fresh objects when rerunning an analysis with different settings. Deadline expiry marks unfinished
tools as timed out; already running worker threads are not forcibly terminated.

| Object | Useful attributes |
|---|---|
| `RunResult` | `request`, `run_id`, `outdir`, `queries`, `timings`, `universe_size`, `universe_source`, `library_manifest`, `warnings` |
| `QueryResult` | `query`, `tool_results`, `records`, `consensus`, `themes`, `key_genes`, `figures`, `warnings` |
| `AdapterResult` | `tool`, `status`, `message`, `runtime_s`, `records`, `databases`, `family`, `version`, `method`, `extra` |

Tables are pandas DataFrames. `records` contains harmonized per-tool tests. `consensus` includes
`term_key`, `term_name`, `database`, `tier_label`, `consensus_score`, `best_padj`, and `hit_genes`.
`themes` groups related terms; `key_genes` ranks genes by theme membership. Tables can be empty,
and empty tables may lack the usual populated-result columns. Check `.empty` before selecting columns.

## Progress and errors

```python
def progress(event):
    if event["type"] == "tool_done":
        print(event["tool"], event["status"], event.get("message", ""))

run = run_analysis(request, outdir="results/progress-example", on_event=progress)
for result in run.queries:
    for tool in result.tool_results:
        if tool.status in {"error", "timeout", "skipped"}:
            print(tool.tool, tool.status, tool.message)
```

Statuses are `ok`, `empty`, `error`, `timeout`, and `skipped`. Inspect them and both run-level and
query-level warnings before interpreting results. A completed run can have partial or no tool evidence.
Missing compiled libraries raise `FileNotFoundError`; input parsing, preparation, and output errors
can propagate to the caller. `local_only=True` applies to analysis adapters; initial database building
still needs network access.
