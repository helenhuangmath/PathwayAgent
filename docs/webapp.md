# Web app and REST API

## Web app

The [Launch web app](launch.md) page opens the app when it is running, and on macOS
`PathwayAgent.command` starts it with a double-click. From a terminal:

```bash
pathwayagent serve                      # opens http://127.0.0.1:8765
pathwayagent serve --port 9000 --no-browser
```

Paste genes or a table, or upload a file; choose the organism, preset and databases; and run. The
page shows progress per tool and then the same report that the command line writes. The
**Or load test data** menu offers bundled example inputs:

| Example | What it shows |
|---|---|
| TEX modules G1–G6 | Five gene lists, one per column. The tutorial data. |
| Modules 1–3 | Three columns under a header row. |
| Module 4 with Wald statistic | A ranked list, which enables the rank-based tools. |
| Module 6 | A single gene list; over-representation only. |

The server runs on your machine and listens on `127.0.0.1` by default. Web-service tools still send
genes to their providers unless **local only** is selected.

## REST API

The web app is a client of this API, so anything it does can be scripted.

| Method | Path | Purpose |
|---|---|---|
| `GET` | `/api/v1/health` | Liveness and version. The only endpoint readable from other websites; the launch page uses it. |
| `GET` | `/api/v1/options` | Organisms, databases, tools, presets and defaults. |
| `POST` | `/api/v1/jobs` | Start a job from JSON. Returns `{"job_id": ...}`. |
| `POST` | `/api/v1/jobs/upload` | Start a job from an uploaded file (`file`, plus `params` as a JSON string). |
| `GET` | `/api/v1/jobs/{job_id}` | Status and results of a job. |
| `GET` | `/api/v1/jobs/{job_id}/events` | Progress as a server-sent event stream. |

### Start a job

```bash
curl -X POST localhost:8765/api/v1/jobs -H 'Content-Type: application/json' \
     -d '{"text":"Fos\nJun\nEgr1\nNr4a1","preset":"standard"}'
curl localhost:8765/api/v1/jobs/<job_id>
```

Upload a file instead:

```bash
curl -X POST localhost:8765/api/v1/jobs/upload \
     -F file=@lists.csv -F 'params={"preset":"quick","local_only":true}'
```

### Job parameters

| Field | Default | Meaning |
|---|---|---|
| `text` | | Pasted genes or table; the format is detected as for files. |
| `lists` | | Alternative to `text`: `[{"name": ..., "genes": [...], "scores": [...]}]`. |
| `name` | `"query"` | Name of the list when `text` holds a single list. |
| `organism` | `"auto"` | `auto`, `mmusculus` or `hsapiens`. |
| `preset` | `"standard"` | `quick`, `standard` or `everything`. |
| `tools` | | List of tool names; overrides `preset`. |
| `databases` | all eight | List of database keys. |
| `universe` | | Background genes as text. |
| `min_size`, `max_size` | `10`, `500` | Gene-set size bounds. |
| `alpha` | `0.05` | Adjusted P threshold. |
| `deadline_s` | `60` | Time budget for the tool stage. |
| `local_only` | `false` | Skip tools that send genes to external services. |
| `genes_of_interest` | | Genes to highlight in the key-gene table. |

For analyses inside Python, the [Python API](PYTHON_API.md) is more direct than the REST API.
