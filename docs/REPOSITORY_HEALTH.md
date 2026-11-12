# Repository Health and Maintenance

PathwayAgent [v1.0.1](https://github.com/helenhuangmath/PathwayAgent/releases/tag/v1.0.1)
is the current formal release. It is published as a regular GitHub release, not a draft or prerelease.
The tasks below track ongoing maintenance and do not change that release status.

## Completed release work

- Package, CLI, local web app, and REST API are published in the GitHub repository.
- Package metadata, runtime version, citation, and current API documentation use version 1.0.1.
- MIT license and citation metadata are included; repository links point to the actual project.
- The release includes a wheel, source archive, worked-example outputs, and SHA-256 checksums.
- The Python API is documented in [PYTHON_API.md](PYTHON_API.md).
- A small executable example is documented in [WORKED_EXAMPLE.md](WORKED_EXAMPLE.md).
- All 20 offline tests passed locally with the mouse reference library, including R agreement.
- GitHub CI passed on Python 3.11 and 3.12 and verified the installed wheel on Linux.
- A synthetic package smoke test checks enrichment, reports, API health, and bundled assets
  without downloading reference data.
- Generated libraries, result folders, virtual environments, and caches are excluded from Git.

## Ongoing maintenance

- Verify installation on additional operating systems and Python versions.
- Exercise live external-service integrations and fresh reference-library downloads separately
  from the offline suite; document the versions and dates used for each release.
- Keep R/Bioconductor setup instructions current and consider providing a container image.
- Extend reproducibility records for reference data and dependency versions.
- Add a DOI to `CITATION.cff` if a release is archived on Zenodo.
- Evaluate persistent job storage before offering the local web app as a shared hosted service.

Reference-library-dependent tests skip explicitly when the compiled library is absent. CI's
synthetic smoke test checks package operation; the full scientific validation suite still requires
the reference library. See [testing/README.md](https://github.com/helenhuangmath/PathwayAgent/blob/main/testing/README.md) and
[RELEASING.md](RELEASING.md) for verification instructions.
