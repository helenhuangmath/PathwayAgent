# Tests and examples

All development tests, example inputs, and archived prototypes are grouped here.

| Folder | Purpose |
|---|---|
| `unit/` | Offline unit tests |
| `validation/` | Planted signals, null calibration, stability, and R agreement |
| `live/` | External-service integration and end-to-end tests |
| `gene_lists/` | Small gene lists and ranked inputs for smoke tests |
| `tcell_modules/` | CD8 T cell module tables used in the tutorial |
| `metascape_reference/` | Metascape reference results used in the tutorial |
| `browser_prototype/` | Standalone browser prototype and its browser tests; its `data/` and `data_source/` (23 MB) are symlinks into `../large_data/browser_prototype/` and are not committed |
| `python_prototype/` | Archived Python prototype with its original tests |

Run the active offline suite from the project root:

```bash
python -m pytest -m "not live"
```

Library-dependent tests skip explicitly when a compiled mouse library is absent. A clean CI checkout
therefore exercises the database-independent tests; full scientific validation requires building or
providing the reference library. Use `-ra` to inspect skipped tests. The v1.0.1 release was checked
locally with the mouse reference library, including the R agreement test.

`package_smoke.py` independently creates synthetic genes and sets in a temporary folder and checks
enrichment, HTML/PDF/Excel generation, API health, bundled examples, and the packaged R worker. CI runs
it against a freshly installed wheel from outside the checkout. It requires `httpx`, installed by the
development extra, but no reference downloads, web-service calls, or R installation:

```bash
python /path/to/checkout/testing/package_smoke.py
```

Run the live integration tests (requires network access and a built library):

```bash
python -m pytest -m live
```

Pytest's default collection targets `unit/`, `validation/`, and `live/`. Archived prototype tests are retained for reference.
The shared fixtures in `conftest.py` resolve example paths from the project root.

Try a small local analysis after building the library:

```bash
python -m pathwayagent run testing/gene_lists/genels_ex_module123.txt --preset quick --local-only
```

The browser prototype can be opened at `testing/browser_prototype/index.html`; its checks are in
`testing/browser_prototype/test.html`. Keep the browser folder's internal structure intact so its
relative script and data links continue to work.

The installed web app's bundled demonstration inputs remain in
`pathwayagent/server/static/examples/` as runtime package assets.
Generated analysis outputs go to `../large_data/results/` (or the ignored `results/` folder when no `large_data/` sibling exists).
