# Repository Health Roadmap

This file tracks the work needed to make PathwayAgent feel like a mature, reusable GitHub package.

## Current State

- Core package: implemented under `pathwayagent/`.
- CLI and local server: implemented.
- Offline tests: implemented and passing locally.
- Documentation: README and tutorial are substantial.
- Generated artifacts: currently present locally, but should stay out of Git.

## Before Public GitHub Release

- Initialize or restore a top-level Git repository.
- ~~Choose and add a top-level open-source license.~~ Done: MIT (`LICENSE`), author and citation metadata added 2026-11-09.
- ~~Replace placeholder GitHub URLs in `pyproject.toml`.~~ Done for v0.1.0: verified `helenhuangmath/PathwayAgent`.
- ~~Keep `.venv/`, `library/`, and `results/` untracked.~~ Done 2026-11-09: large data moved to the sibling `large_data/` folder; the code finds it automatically.
- Test code, inputs, and archived prototypes are grouped under `testing/`; the active Python tests are in `testing/unit/`, `testing/validation/`, and `testing/live/`.
- ~~Decide whether `PeakAnnoAgent/` should move to a separate repository.~~ Done: it now lives beside this project.
- ~~Add a small, stable test fixture for docs and CI.~~ Synthetic installed-package smoke test plus a packaged worked example added for v0.1.0.
- Verify installation on additional clean operating systems; v0.1.0 wheel and source distribution are checked in fresh Python environments on macOS, with Linux wheel checks configured in CI.

## Maturity Milestones

### Alpha

- Package installs with `pip install -e .`.
- Offline tests pass.
- README and tutorial explain normal use.
- Generated data is not committed.

### Beta

- CI runs offline tests on pull requests.
- Live tests are documented and optionally scheduled.
- R dependencies are documented or containerized.
- Large libraries are downloadable or generated reproducibly.
- Example outputs are available as release artifacts.

### Stable

- Public API is documented in `docs/PYTHON_API.md` (alpha contract; may change before 1.0).
- Versioning follows semantic versioning.
- Releases include changelog entries.
- Citation metadata is complete (`CITATION.cff` added; add a DOI once archived on Zenodo).
- Reproducibility story is clear for Python, R, and database versions.
