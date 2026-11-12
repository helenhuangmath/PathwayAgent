# Releasing PathwayAgent

The current release version is v1.0.1. Publish it as a regular GitHub release. A GitHub release does not change
repository visibility and does not publish to PyPI.

## Prepare and verify

1. Update the version in `pyproject.toml`, `pathwayagent/__init__.py`, and `CITATION.cff` together;
   update the release date, changelog, and versioned example commands.
2. Install the development extra in a supported Python environment. Run `ruff check pathwayagent
   testing/conftest.py testing/unit testing/validation testing/live testing/package_smoke.py`.
3. Run `python -m pytest -m "not live" -ra` with a compiled mouse reference library.
   Without one, library-dependent tests skip. R agreement also requires the documented R environment.
4. Run `python -m build` and `python -m twine check dist/*`.
5. Create fresh virtual environments; install the wheel in one and the source archive in the other.
   From outside the source directory, run `pathwayagent --help` and
   `python /absolute/path/to/testing/package_smoke.py` with each environment's interpreter.
   Confirm imports resolve to the fresh environment, not an editable checkout.
6. Run `python -m pathwayagent.examples.interferon --out /tmp/interferon-example` from the installed
   wheel with `PATHWAYAGENT_LIBRARY` pointing to the reference library. Inspect the report, adapter
   status, and example results. Bundle the generated example outputs separately from source.
7. Check that archives contain the R worker, web assets, and packaged example, and exclude local
   environments, generated libraries, caches, and result folders. Generate SHA-256 checksums for assets.

CI tests database-independent code and installs the built wheel for the synthetic smoke test on Linux.
It does not substitute for reference-library scientific validation or external-service integration tests.

## Publish to GitHub

Commit the verified source, push it to the repository, and create an annotated version tag at that
commit. Do not overwrite an existing release/tag. Prepare release notes describing changes, example
commands, validation results, and limitations. From the repository checkout:

```bash
git tag -a v1.0.1 -m "PathwayAgent v1.0.1"
git push origin main v1.0.1
gh release create v1.0.1 --verify-tag --title "PathwayAgent v1.0.1" \
  --notes-file /path/to/release-notes.md \
  /path/to/pathwayagent-1.0.1-py3-none-any.whl \
  /path/to/pathwayagent-1.0.1.tar.gz \
  /path/to/interferon-example-v1.0.1.zip \
  /path/to/SHA256SUMS
```

Confirm the tag, attached assets, and their checksums after publishing. Live web-service tests and
a fresh reference-library download should be reported separately if not exercised.
