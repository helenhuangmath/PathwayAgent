"""
Project profile (design §11): background notes, genes/pathways of interest.
Used for highlighting and LLM context only — never to change statistics (design §1.6).

Example profile.yaml:
    project: Tex_vs_Tmem
    organism: mmusculus
    background: >
      CD8 T cells from chronic vs acute infection ...
    genes_of_interest: [Tox, Tcf7, Pdcd1]
    pathways_of_interest: ["type I interferon", "TCR signaling"]
"""
from __future__ import annotations

from pathlib import Path

import yaml


def load_profile(path: str | Path | None) -> dict:
    if not path:
        return {}
    d = yaml.safe_load(Path(path).read_text()) or {}
    d.setdefault("genes_of_interest", [])
    d.setdefault("pathways_of_interest", [])
    d.setdefault("background", "")
    return d


def interest_hits(profile: dict, consensus) -> list[str]:
    """Consensus terms matching the profile's pathways of interest (case-insensitive substring)."""
    pats = [p.lower() for p in profile.get("pathways_of_interest", [])]
    if not pats or consensus is None or consensus.empty:
        return []
    names = consensus["term_name"].str.lower()
    return consensus.loc[names.apply(lambda n: any(p in n for p in pats)), "term_key"].tolist()
