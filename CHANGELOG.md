# Changelog

Published release version: `v1.0.1`. Package metadata, citation, Python API documentation,
and installation examples use the same version.

### Added

- Multi-tool pathway analysis pipeline.
- Local and web enrichment adapters.
- Cross-tool consensus scoring, theme grouping, key-gene ranking, figures, HTML report, CLI, REST API, and local web app.
- Offline unit and validation tests.
- Bundled example inputs in the web app ("Or load test data" menu), served from `/examples`.
- Example figures and captions in README and TUTORIAL (`docs/figures/`, 300-dpi PNG).
- Documentation site (MkDocs Material) with introduction, quick start, tutorial, user guide and Python API,
  deployed to GitHub Pages by `.github/workflows/docs.yml`.
- "Launch web app" page on the documentation site that detects a locally running PathwayAgent and opens it,
  and `PathwayAgent.command`, a double-click launcher (macOS) that installs on first use and starts the web app.
  The health endpoint is readable cross-origin for this purpose; no other endpoint is.

### Changed
- Repository split into `github_repo/` (code, tests, docs, small example data) and a sibling `large_data/` folder (gene-set library, results, build artifacts, browser-prototype databases); `pathwayagent.config` finds the sibling automatically, or use `PATHWAYAGENT_DATA`.
- Figure text enlarged by one point (9–10 pt Arial) for readability.
- Web app and HTML report restyled with a flat Bootstrap-3 "Flatly" look (fixed dark navbar, wordmark banner, pill tabs, panels, dashed footer), with Help and About panels and a citation/contact footer.
- MIT license, author and citation metadata (`LICENSE`, `CITATION.cff`, `pyproject.toml`).

### Known Limitations

- Some adapters depend on external web services and can change behavior as those services change.
- R-based tools require a local R/Bioconductor environment.
- Large generated libraries and example outputs should be distributed separately from the source repository.
