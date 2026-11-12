"""
External web-service adapters (design §6): Enrichr, g:Profiler, Reactome, STRING.

Each adapter:
  * requests ALL tested terms where the service allows it (bug B8), so that
    "tested but not significant" is distinguishable from "not tested";
  * maps hit genes back to the user's symbols;
  * tags every row with a canonical term_key (bug B9) and a database type.
"""
from __future__ import annotations

import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

from ..core import cache
from ..core.schema import QueryInput
from .base import Adapter, RunContext, resolve_key

UA = {"User-Agent": "PathwayAgent/0.1 (academic research tool)"}
TIMEOUT = 45
# polite per-service concurrency limits (avoid HTTP 429 when many lists run at once)
_SEM = {"enrichr": threading.Semaphore(2), "gprofiler": threading.Semaphore(3),
        "reactome": threading.Semaphore(3), "string": threading.Semaphore(2)}


def http(service: str, method: str, url: str, retries: int = 4, **kw) -> requests.Response:
    """HTTP call with a per-service semaphore and retry/backoff on 429, 5xx and connection errors."""
    kw.setdefault("timeout", TIMEOUT)
    kw.setdefault("headers", UA)
    delay = 1.0
    for attempt in range(retries + 1):
        try:
            with _SEM[service]:
                r = requests.request(method, url, **kw)
            if r.status_code == 429 or r.status_code >= 500:
                if attempt == retries:
                    r.raise_for_status()
                wait = float(r.headers.get("Retry-After", delay)) if r.headers.get("Retry-After", "").isdigit() else delay
                time.sleep(min(wait, 15))
                delay *= 2
                continue
            r.raise_for_status()
            return r
        except (requests.ConnectionError, requests.Timeout):
            if attempt == retries:
                raise
            time.sleep(delay)
            delay *= 2
    raise RuntimeError("unreachable")


# ─────────────────────────────── Enrichr ─────────────────────────────────────

ENRICHR = "https://maayanlab.cloud/Enrichr/"
_ENRICHR_PATTERNS = {
    # database: (regex for library name, species of library symbols)
    "GO_BP": (r"^GO_Biological_Process_(\d{4})$", "human"),
    "GO_MF": (r"^GO_Molecular_Function_(\d{4})$", "human"),
    "GO_CC": (r"^GO_Cellular_Component_(\d{4})$", "human"),
    "KEGG": (r"^KEGG_(\d{4})(?:_Human)?$", "human"),
    "REACTOME": (r"^Reactome_(?:Pathways_)?(\d{4})$", "human"),
    "WIKIPATHWAYS": (r"^WikiPathways?_(\d{4})_{species}$", "native"),
    "HALLMARK": (r"^MSigDB_Hallmark_(\d{4})$", "human"),
    "DISGENET": (r"^DisGeNET$", "human"),
}


def _enrichr_libraries() -> list[str]:
    def fetch():
        d = http("enrichr", "GET", ENRICHR + "datasetStatistics").json()
        return [x["libraryName"] for x in d["statistics"]]
    libs, _ = cache.cached("enrichr_libs", 1, fetch)
    return libs


def pick_enrichr_library(db: str, organism: str, libs: list[str]) -> tuple[str | None, str]:
    pat, sp = _ENRICHR_PATTERNS[db]
    species = "Mouse" if organism == "mmusculus" else "Human"
    pat = pat.replace("{species}", species)
    best, best_year = None, -1
    for name in libs:
        m = re.match(pat, name)
        if m:
            year = int(m.group(1)) if m.groups() else 0
            if year > best_year:
                best, best_year = name, year
    return best, ("native" if sp == "native" else "human")


class Enrichr(Adapter):
    name, family, method = "enrichr", "F1_ORA", "Fisher exact (Enrichr)"
    web = True
    databases = tuple(_ENRICHR_PATTERNS)
    version = "Enrichr API"

    def _run(self, q: QueryInput, ctx: RunContext):
        libs = _enrichr_libraries()
        human_mode = ctx.organism != "hsapiens"
        # symbol lists: human libraries get orthologs, species libraries get upper-cased symbols
        back: dict[str, dict[str, str]] = {"human": {}, "native": {}}
        lists = {}
        if human_mode:
            orth = []
            for g in q.genes:
                for o in ctx.mapper.ortholog.get(g, []):
                    back["human"].setdefault(o.upper(), g)
                    orth.append(o)
            lists["human"] = sorted(set(orth))
            lists["native"] = sorted({g.upper() for g in q.genes})
            back["native"] = {g.upper(): g for g in q.genes}
        else:
            lists["human"] = lists["native"] = list(q.genes)
            back["human"] = back["native"] = {g.upper(): g for g in q.genes}

        chosen = {}
        for db in self.covered(ctx):
            lib, sp = pick_enrichr_library(db, ctx.organism, libs)
            if lib:
                chosen[db] = (lib, sp)
        list_ids: dict[str, int] = {}
        id_lock = threading.Lock()

        def list_id(sp: str) -> int:
            """Upload each list once (bug B14), lazily — skipped entirely on a full cache hit."""
            with id_lock:
                if sp not in list_ids:
                    r = http("enrichr", "POST", ENRICHR + "addList",
                             files={"list": (None, "\n".join(lists[sp])), "description": (None, "PathwayAgent")})
                    list_ids[sp] = r.json()["userListId"]
                return list_ids[sp]

        def enrich(db):
            lib, sp = chosen[db]
            payload = {"lib": lib, "genes": lists[sp]}
            def fetch():
                r = http("enrichr", "GET", ENRICHR + "enrich",
                         params={"userListId": list_id(sp), "backgroundType": lib})
                return r.json().get(lib, [])
            rows, _ = cache.cached("enrichr", payload, fetch)
            return db, lib, sp, rows

        recs = []
        used = {}
        with ThreadPoolExecutor(max_workers=8) as pool:
            for db, lib, sp, rows in pool.map(enrich, list(chosen)):
                used[db] = lib
                n_list = len(lists[sp])
                for row in rows:
                    _, term, p, odds, comb, genes, padj = row[:7]
                    key, lrow = resolve_key(ctx.library, db, "", term)
                    hits = sorted({back[sp].get(g.upper(), g) for g in genes})
                    recs.append(self.base_record(
                        q, database=db, term_id=term, term_key=key, term_name=re.sub(
                            r"\s*\(GO:\d+\)$|\s+R-HSA-\d+$|\s+WP\d+$", "", term),
                        term_size=ctx.library.sets["size"].iat[lrow] if lrow is not None else None,
                        query_size=n_list, overlap=len(genes), pval=p, padj=padj,
                        effect=odds, effect_type="odds_ratio", direction=0, hit_genes=hits))
        extra = {"libraries": used,
                 "note": "mouse genes mapped to human orthologs for human-symbol libraries" if human_mode else ""}
        return pd.DataFrame(recs), extra


# ─────────────────────────────── g:Profiler ──────────────────────────────────

GPROFILER = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
_GP_SOURCES = {"GO_BP": "GO:BP", "GO_MF": "GO:MF", "GO_CC": "GO:CC", "KEGG": "KEGG",
               "REACTOME": "REAC", "WIKIPATHWAYS": "WP"}
_GP_BACK = {v: k for k, v in _GP_SOURCES.items()}


class GProfiler(Adapter):
    name, family, method = "gprofiler", "F1_ORA", "hypergeometric (g:Profiler g:GOSt)"
    web = True
    databases = tuple(_GP_SOURCES)
    version = "g:Profiler API"
    ordered = False

    def _run(self, q: QueryInput, ctx: RunContext):
        genes = q.ranked_genes() if self.ordered else list(q.genes)
        payload = {"organism": ctx.organism, "query": genes,
                   "sources": [_GP_SOURCES[d] for d in self.covered(ctx)],
                   "user_threshold": 1.0, "all_results": True, "no_evidences": False,
                   "significance_threshold_method": "fdr", "ordered": self.ordered,
                   "domain_scope": "annotated"}
        if ctx.user_universe:
            payload["domain_scope"] = "custom"
            payload["background"] = sorted(ctx.universe_symbols)

        def fetch():
            r = http("gprofiler", "POST", GPROFILER, json=payload, timeout=TIMEOUT * 2)
            return r.json()
        data, _ = cache.cached("gprofiler", payload, fetch)
        md = data["meta"]["genes_metadata"]["query"]["query_1"]
        ensgs = md["ensgs"]
        ensg2in = {}
        for inp, ids in md["mapping"].items():
            for e in ids:
                ensg2in.setdefault(e, inp)
        recs = []
        for r in data.get("result", []):
            db = _GP_BACK.get(r["source"])
            if db is None:
                continue
            # intersections: one entry per query ENSG (evidence codes), NOT gene names (bug B2)
            hits = sorted({ensg2in.get(ensgs[i], ensgs[i]) for i, ev in enumerate(r.get("intersections", []))
                           if ev})
            key, _ = resolve_key(ctx.library, db, r["native"], r["name"])
            recs.append(self.base_record(
                q, database=db, term_id=r["native"], term_key=key, term_name=r["name"],
                term_size=r["term_size"], query_size=r["query_size"],
                universe_size=r.get("effective_domain_size"), overlap=r["intersection_size"],
                pval=None, padj=r["p_value"],   # g:Profiler reports ADJUSTED p only (bug B3)
                padj_method="g:Profiler FDR", effect=(r["intersection_size"] / max(r["query_size"], 1)) /
                (r["term_size"] / max(r.get("effective_domain_size", 1), 1)),
                effect_type="fold_enrichment", direction=0, hit_genes=hits))
        df = pd.DataFrame(recs)
        if len(df):
            df["pval"] = df["padj"]  # raw p not returned; rank on adjusted value
        return df, {"version": data["meta"].get("version", "")}


class GProfilerOrdered(GProfiler):
    """g:Profiler ordered query = minimum-hypergeometric family (bug B4)."""
    name, family, method = "gprofiler_ordered", "F3_RANK", "ordered query / mHG (g:Profiler)"
    needs_scores = True
    ordered = True


# ─────────────────────────────── Reactome ────────────────────────────────────

REACTOME = "https://reactome.org/AnalysisService/"


class Reactome(Adapter):
    name, family, method = "reactome", "F1_ORA", "hypergeometric (Reactome AnalysisService)"
    web = True
    databases = ("REACTOME",)
    version = "Reactome AnalysisService"

    def _run(self, q: QueryInput, ctx: RunContext):
        genes = list(q.genes)
        payload = {"genes": genes, "organism": ctx.organism}

        def fetch():
            url = (REACTOME + "identifiers/projection?interactors=false&pageSize=5000&page=1"
                   "&sortBy=ENTITIES_PVALUE&order=ASC&resource=TOTAL&includeDisease=true")
            r = http("reactome", "POST", url, data="\n".join(genes), headers={"Content-Type": "text/plain", **UA})
            res = r.json()
            token = res["summary"]["token"]
            ids = [p["stId"] for p in res.get("pathways", [])]
            found = {}
            for i in range(0, len(ids), 200):  # found genes per pathway (bug B5)
                f = http("reactome", "POST", f"{REACTOME}token/{token}/found/all?resource=TOTAL",
                         data=",".join(ids[i:i + 200]), headers={"Content-Type": "text/plain", **UA})
                if f.ok:
                    for p in f.json():
                        found[p["pathway"]] = sorted({e["id"] for e in p.get("entities", [])})
            return {"pathways": res.get("pathways", []), "found": found,
                    "not_found": res.get("identifiersNotFound", 0), "token": token}
        data, _ = cache.cached("reactome", payload, fetch, store_if=lambda d: bool(d.get("pathways")))
        n_q = len(genes) - int(data.get("not_found", 0))
        recs = []
        for p in data["pathways"]:
            e = p["entities"]
            key, lrow = resolve_key(ctx.library, "REACTOME", p["stId"], p["name"])
            hits = data["found"].get(p["stId"], [])
            recs.append(self.base_record(
                q, database="REACTOME", term_id=p["stId"], term_key=key, term_name=p["name"],
                term_size=e.get("total"), query_size=n_q, overlap=e.get("found"),
                pval=e.get("pValue"), padj=e.get("fdr"), padj_method="Reactome FDR (BH)",
                effect=(e.get("found", 0) / max(n_q, 1)) / max(e.get("ratio", 1e-9), 1e-9),
                effect_type="fold_enrichment", direction=0, hit_genes=hits))
        return pd.DataFrame(recs), {"token": data.get("token"),
                                    "note": "projected to human pathways (Reactome projection)"}


# ─────────────────────────────── STRING ──────────────────────────────────────

STRING = "https://string-db.org/api/json/"
_STRING_CAT = {"Process": "GO_BP", "Function": "GO_MF", "Component": "GO_CC", "KEGG": "KEGG",
               "RCTM": "REACTOME", "WikiPathways": "WIKIPATHWAYS"}


class StringDB(Adapter):
    """STRING functional enrichment (hypergeometric vs STRING background). Returns significant terms only."""
    name, family, method = "string", "F1_ORA", "hypergeometric (STRING)"
    web = True
    databases = tuple(_STRING_CAT.values())
    version = "STRING API"

    def _run(self, q: QueryInput, ctx: RunContext):
        taxon = 10090 if ctx.organism == "mmusculus" else 9606
        payload = {"identifiers": "\r".join(q.genes), "species": taxon, "caller_identity": "pathwayagent"}

        def fetch():
            enr = http("string", "POST", STRING + "enrichment", data=payload).json()
            try:
                ppi = http("string", "POST", STRING + "ppi_enrichment", data=payload, retries=1).json()
            except requests.RequestException:
                ppi = []
            return {"enrichment": enr, "ppi": ppi}
        data, _ = cache.cached("string", payload, fetch)
        recs = []
        for r in data["enrichment"]:
            db = _STRING_CAT.get(r["category"])
            if db not in self.covered(ctx):
                continue
            key, lrow = resolve_key(ctx.library, db, r["term"], r["description"])
            recs.append(self.base_record(
                q, database=db, term_id=r["term"], term_key=key, term_name=r["description"],
                term_size=r["number_of_genes_in_background"], overlap=r["number_of_genes"],
                query_size=len(q.genes), pval=r["p_value"], padj=r["fdr"], padj_method="STRING FDR",
                effect=None, effect_type="", direction=0, hit_genes=sorted(r.get("inputGenes", []))))
        ppi = data.get("ppi") or [{}]
        return pd.DataFrame(recs), {"ppi_enrichment": ppi[0] if ppi else {},
                                    "note": "STRING returns significant terms only"}
