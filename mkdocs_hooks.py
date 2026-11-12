"""MkDocs hook: include files from the repository root in the documentation site.

A page in ``docs/`` that contains a line of the form

    <!-- include: TUTORIAL.md -->

has that line replaced by the named file (path relative to the repository root). Image
and link targets written for the repository root (``docs/figures/x.png``) are rewritten
to be relative to ``docs/``. This keeps TUTORIAL.md and CHANGELOG.md as the single
source for both GitHub and the site.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INCLUDE = re.compile(r"^<!--\s*include:\s*(\S+)\s*-->\s*$", re.MULTILINE)


def _load(match: re.Match) -> str:
    path = (ROOT / match.group(1)).resolve()
    if ROOT not in path.parents or not path.is_file():
        raise FileNotFoundError(f"mkdocs include target not found in repository: {match.group(1)}")
    text = path.read_text(encoding="utf-8")
    return re.sub(r"\]\(docs/", "](", text)


def on_page_markdown(markdown: str, **_kwargs) -> str:
    return INCLUDE.sub(_load, markdown)
