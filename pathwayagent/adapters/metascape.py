"""
Metascape adapter (design §6).

Metascape has no public REST API, so its results reach PathwayAgent in one of two ways:

  import  a result package from metascape.org ("All in One Zip File", zipped or unpacked) or an MSBio
          output folder, passed as AnalysisRequest.metascape = {list name or "*": path}
          (CLI: --metascape [LIST=]PATH). A multi-list package serves every list it contains.
  MSBio   the licensed Docker distribution of Metascape, used when PATHWAYAGENT_MSBIO points at an
          MSBio installation whose container is running (bin/up.sh) and no package was supplied.

Both give the same files. Enrichment_GO/GO_AllLists.csv holds every enriched term of every list;
Enrichment_GO/_FINAL_GO.csv adds Metascape's own term clusters (GROUP_ID). Metascape reports enriched
terms only (its defaults: P < 0.01, >= 3 hit genes, enrichment > 1.5), like STRING.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import re
import subprocess
import zipfile
from pathlib import Path

import pandas as pd

from ..core.schema import AdapterResult, QueryInput
from .base import Adapter, RunContext, resolve_key

_CATEGORY = {"GO Biological Processes": "GO_BP", "GO Molecular Functions": "GO_MF",
             "GO Cellular Components": "GO_CC", "KEGG Pathway": "KEGG", "Reactome Gene Sets": "REACTOME",
             "WikiPathways": "WIKIPATHWAYS", "Hallmark Gene Sets": "HALLMARK"}
_GO_ALL, _GO_FINAL = "Enrichment_GO/GO_AllLists.csv", "Enrichment_GO/_FINAL_GO.csv"
MSBIO_DIR = os.environ.get("PATHWAYAGENT_MSBIO", "")
MSBIO_CONTAINER = os.environ.get("PATHWAYAGENT_MSBIO_CONTAINER", "msbio1")
MSBIO_TIMEOUT = 900


# ───────────────────────────── result packages ─────────────────────────────

def _read_member(package: Path, member: str) -> pd.DataFrame | None:
    """Read one CSV of a Metascape package (folder or zip; the zip may have a top-level folder)."""
    if package.is_dir():
        hits = [package / member] if (package / member).exists() else sorted(package.glob(f"*/{member}"))
        return pd.read_csv(hits[0]) if hits and hits[0].exists() else None
    with zipfile.ZipFile(package) as z:
        names = [n for n in z.namelist() if n == member or n.endswith("/" + member)]
        return pd.read_csv(io.BytesIO(z.read(names[0]))) if names else None


def read_package(package: str | Path) -> pd.DataFrame:
    """All enriched terms of a Metascape package: one row per (GeneList, term), with GROUP_ID when known."""
    package = Path(package).expanduser()
    if not package.exists():
        raise FileNotFoundError(f"Metascape result not found: {package}")
    if package.is_file() and package.suffix.lower() != ".zip":
        package = package.parent                      # e.g. .../metascape_result.xlsx -> its folder
    terms = _read_member(package, _GO_ALL)
    if terms is None:
        raise ValueError(f"{package} has no {_GO_ALL}: pass the complete Metascape package "
                         f"('All in One Zip File'), not only the spreadsheet")
    final = _read_member(package, _GO_FINAL)
    if final is not None and "GROUP_ID" in final.columns:
        terms = terms.merge(final[["GO", "GROUP_ID"]].drop_duplicates("GO"), on="GO", how="left")
    else:
        terms["GROUP_ID"] = pd.NA
    return terms


def _pick_list(terms: pd.DataFrame, name: str) -> pd.DataFrame | None:
    """Rows of the list called `name`; a single-list package is used whatever its list is called."""
    lists = terms["GeneList"].astype(str)
    if lists.nunique() == 1:
        return terms
    wanted = re.sub(r"[^0-9a-z]+", "", name.lower())
    match = lists.map(lambda s: re.sub(r"[^0-9a-z]+", "", s.lower())) == wanted
    return terms[match] if match.any() else None


# ───────────────────────────── MSBio (Docker) ──────────────────────────────

def msbio_ready() -> str | None:
    """None when MSBio can be run, else the reason it cannot."""
    if not MSBIO_DIR:
        return "PATHWAYAGENT_MSBIO is not set"
    if not (Path(MSBIO_DIR) / "data").is_dir():
        return f"{MSBIO_DIR} is not an MSBio installation (no data/ folder)"
    try:
        out = subprocess.run(["docker", "ps", "--filter", f"name=^{MSBIO_CONTAINER}$", "--format", "{{.Names}}"],
                             capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"docker is not available ({exc})"
    if MSBIO_CONTAINER not in out.stdout.split():
        return f"MSBio container '{MSBIO_CONTAINER}' is not running (start it with bin/up.sh)"
    return None


def run_msbio(q: QueryInput, organism: str) -> Path:
    """Run one list through MSBio and return the output folder. Results are reused for an identical list."""
    tag = hashlib.sha1(("\n".join(sorted(q.genes)) + organism).encode()).hexdigest()[:12]
    work = Path(MSBIO_DIR) / "data" / "pathwayagent" / tag
    out = work / "output"
    if (out / _GO_ALL).exists():
        return out
    work.mkdir(parents=True, exist_ok=True)
    (work / "genes.csv").write_text("Gene\n" + "\n".join(q.genes) + "\n")
    inside = f"/data/pathwayagent/{tag}"
    job = {"input": f"{inside}/genes.csv", "output": f"{inside}/output", "single": True,
           "source_tax_id": "mouse" if organism == "mmusculus" else "human"}
    (work / "run.job").write_text(json.dumps(job) + "\n")
    cmd = ["docker", "exec", "-w", "/home/meta_test", MSBIO_CONTAINER, "/home/meta_test/mylib/ms/msbio.py",
           "--license=/license", f"{inside}/run.job"]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=MSBIO_TIMEOUT)
    if not (out / _GO_ALL).exists():
        tail = (p.stderr or p.stdout or "").strip().splitlines()[-3:]
        raise RuntimeError("MSBio produced no enrichment table: " + " | ".join(tail))
    return out


# ───────────────────────────── adapter ─────────────────────────────────────

class Metascape(Adapter):
    """Metascape enrichment (hypergeometric, BH q-values), from an imported package or a local MSBio run."""
    name, family, method = "metascape", "F1_ORA", "hypergeometric (Metascape)"
    databases = tuple(_CATEGORY.values())
    version = "Metascape"

    @staticmethod
    def _package(q: QueryInput, ctx: RunContext) -> str | None:
        given = ctx.request.metascape or {}
        return given.get(q.name) or given.get("*")

    def applicable(self, q: QueryInput, ctx: RunContext) -> str | None:
        reason = super().applicable(q, ctx)
        if reason or self._package(q, ctx):
            return reason
        why = msbio_ready()
        if why:
            return ("no Metascape result for this list: analyse it at metascape.org, download 'All in One "
                    f"Zip File' and pass --metascape {q.name}=PATH (MSBio not used: {why})")
        return None

    def run(self, q: QueryInput, ctx: RunContext) -> AdapterResult:
        res = super().run(q, ctx)
        if res.status == "ok":      # only the ontologies present in the package count as covered
            res.databases = [d for d in res.databases if d in set(res.records["database"])]
        return res

    def _run(self, q: QueryInput, ctx: RunContext):
        package = self._package(q, ctx)
        source = "imported package"
        if not package:
            package, source = run_msbio(q, ctx.organism), "MSBio"
        terms = _pick_list(read_package(package), q.name)
        if terms is None:
            raise ValueError(f"{package} has no gene list named '{q.name}'")

        # hit genes come back as the analysed species' symbols: map them onto the user's symbols
        back = {g.upper(): g for g in q.genes}
        for g in q.genes:
            for o in ctx.mapper.ortholog.get(g, []):
                back.setdefault(o.upper(), g)
        recs, skipped, clusters = [], set(), []
        for d in terms.to_dict("records"):
            db = _CATEGORY.get(str(d["Category"]))
            if db not in self.covered(ctx):
                skipped.add(str(d["Category"]))
                continue
            name = re.sub(r"^HALLMARK[_ ]", "", str(d["Description"]))
            key, _ = resolve_key(ctx.library, db, str(d["GO"]), name)
            hits = sorted({back.get(h.upper(), h) for h in str(d["Hits"]).split("|") if h})
            recs.append(self.base_record(
                q, database=db, term_id=str(d["GO"]), term_key=key, term_name=str(d["Description"]),
                term_size=d.get("#GeneInGO"), query_size=d.get("#GeneInHitList"),
                universe_size=d.get("#TotalGeneInLibrary"), overlap=d.get("#GeneInGOAndHitList"),
                pval=10 ** float(d["LogP"]), padj=10 ** min(float(d["Log(q-value)"]), 0.0),
                padj_method="Metascape q (BH)",
                effect=d["Enrichment"], effect_type="fold_enrichment", direction=0, hit_genes=hits))
            if pd.notna(d.get("GROUP_ID")):         # Metascape's own clustering, kept for the comparison figure
                clusters.append({"GROUP_ID": int(d["GROUP_ID"]), "term_key": key, "Description": str(d["Description"]),
                                 "LogP": float(d["LogP"]), "hit_genes": ";".join(hits)})
        return pd.DataFrame(recs), {"source": source, "package": str(package),
                                    "categories_not_used": sorted(skipped), "clusters": clusters,
                                    "note": "Metascape returns enriched terms only"}
