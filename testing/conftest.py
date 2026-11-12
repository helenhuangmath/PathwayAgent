import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
EXAMPLES = ROOT / "testing" / "gene_lists"


def require_library(organism="mmusculus"):
    from pathwayagent.config import library_path
    path = library_path(organism)
    if not (path / "sets.parquet").exists() or not (path / "genes.parquet").exists():
        pytest.skip(f"Reference library absent; build with: pathwayagent db build --organism {organism}")


@pytest.fixture(scope="session")
def lib():
    require_library()
    from pathwayagent.core.library import get_library
    return get_library("mmusculus")


@pytest.fixture(scope="session")
def mapper():
    require_library()
    from pathwayagent.core.idmap import get_mapper
    return get_mapper("mmusculus")


@pytest.fixture(scope="session")
def examples():
    return EXAMPLES


def make_ctx(queries, tools=("local_ora",), organism="mmusculus", **kw):
    require_library(organism)
    from pathwayagent.core.schema import AnalysisRequest
    from pathwayagent.pipeline import prepare
    req = AnalysisRequest(queries=queries, organism=organism, tools=list(tools), **kw)
    ctx, warnings = prepare(req)
    return req, ctx
