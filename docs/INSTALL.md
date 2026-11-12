# Installation

## Quick Local Install

```bash
python -m venv .venv
.venv/bin/pip install -e ".[dev]"
.venv/bin/python -m pathwayagent db build
```

The database build creates the gene-set library (several hundred MB). It is generated state and is never committed to Git.

Where it goes: if a folder named `large_data/` exists next to the repository folder (the development layout, see below), the library is built in `large_data/library/` and run outputs go to `large_data/results/`. Otherwise they go to `library/` and `results/` inside the repository, both ignored by Git.

```
PathwayAgent/
├── github_repo/     # this repository: code, tests, docs, small example inputs
└── large_data/      # library/, results/, dist/, browser_prototype/ databases
```

## Runtime Paths

Environment variables override the locations above:

```bash
export PATHWAYAGENT_DATA=/path/to/large_data              # parent of library/ and results/
export PATHWAYAGENT_LIBRARY=/path/to/pathwayagent-library
export PATHWAYAGENT_RESULTS=/path/to/pathwayagent-results
export PATHWAYAGENT_CACHE=/path/to/pathwayagent-cache.sqlite
export PATHWAYAGENT_RSCRIPT=/path/to/Rscript
export MPLCONFIGDIR=/tmp/pathwayagent-matplotlib
```

For wheel installations, set explicit writable library and results paths before building the database
or importing the Python API. Without overrides, the fallback paths are relative to the installed
package location, which may not be writable in a shared Python environment.

To install a versioned wheel, see the [small worked example](WORKED_EXAMPLE.md).

`MPLCONFIGDIR` is useful on servers or sandboxes where the default Matplotlib cache directory is not writable.

## Verify

```bash
.venv/bin/python -m pathwayagent db list
.venv/bin/python -m pytest -m "not live"
.venv/bin/python -m pathwayagent run testing/gene_lists/genels_ex_module123.txt --preset quick --local-only --out /tmp/pathwayagent-smoke
```
