# Contributing

Thanks for helping improve PathwayAgent. 

## Development Setup

```bash
python -m venv .venv
.venv/bin/pip install -e ".[dev]"
.venv/bin/python -m pathwayagent db build
```

The local gene-set library is generated under `../large_data/library/` when that sibling folder exists, else under `library/`. It is intentionally ignored by Git because it is large and reproducible. Set `PATHWAYAGENT_LIBRARY=/path/to/library` if you want to share one library across clones.

## Tests

Run the offline suite before opening a pull request:

```bash
.venv/bin/python -m pytest -m "not live"
```

Live tests call external services and are best run manually:

```bash
.venv/bin/python -m pytest -m live
```

R-dependent tests require `Rscript` plus the Bioconductor packages documented in the README. Set `PATHWAYAGENT_RSCRIPT` if R is not at `/usr/local/bin/Rscript`.

## Pull Request Checklist

- Add or update tests for behavior changes.
- Keep generated outputs, local libraries, and virtual environments out of Git.
- Update `README.md` or `TUTORIAL.md` when user-facing commands or outputs change.
- Run `pytest -m "not live"` locally.
- Note whether live web-service tests were run.

## Release Checklist

- Bump `version` in `pyproject.toml`, `pathwayagent/__init__.py` and `CITATION.cff`; date the entry in `CHANGELOG.md`.
- Run the offline test suite.
- Run live tests when web services are reachable.
- Build the package with `python -m build`.
- Check the distribution with `twine check dist/*`.
- Tag the release and attach any large prebuilt library or example-result artifacts separately (kept locally in `../large_data/`).

See [docs/RELEASING.md](docs/RELEASING.md) for artifact verification and GitHub release commands.
