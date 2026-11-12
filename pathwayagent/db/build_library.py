"""
Build the local gene-set library (design §5).

Downloads every source into library/raw/, then compiles, per organism:
  library/<organism>/genes.parquet      gene table: symbol, entrez, ensembl, synonyms, orthologs
  library/<organism>/sets.parquet       one row per gene set (id, term_key, database, name, size, ...)
  library/<organism>/membership.npz     sparse sets x genes boolean matrix
  library/<organism>/gene_index.txt     column order of the matrix
  library/<organism>/gmt/<DB>.gmt       per-database GMT (set_id, name, genes) for the R worker
  library/<organism>/manifest.json      source URLs, versions, download dates, counts

Run:  python -m pathwayagent db build [--organism mmusculus hsapiens] [--force]
"""
from __future__ import annotations

import json
import re
import time
from collections import defaultdict
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import sparse

from ..config import LIBRARY_DIR, ORGANISMS, library_path

RAW = LIBRARY_DIR / "raw"
MSIGDB_BASE = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
WP_BASE = "https://data.wikipathways.org/current/gmt/"
MGI = "https://www.informatics.jax.org/downloads/reports/"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
GO_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"
REACTOME_NAMES_URL = "https://reactome.org/download/current/ReactomePathways.txt"
KEGG_REST = "https://rest.kegg.jp/"
ENRICHR_LIB = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="

MIN_SET_SIZE = 3


# ───────────────────────────── download helpers ─────────────────────────────

def fetch(url: str, name: str, force: bool = False, timeout: int = 300) -> Path:
    """Download url to library/raw/name once (re-used unless force=True)."""
    RAW.mkdir(parents=True, exist_ok=True)
    dest = RAW / name
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return dest
    t0 = time.time()
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
        tmp.rename(dest)
    print(f"  downloaded {name} ({dest.stat().st_size / 1e6:.1f} MB, {time.time() - t0:.1f}s)")
    return dest


def latest_msigdb_release(species: str) -> str:
    html = requests.get(MSIGDB_BASE, timeout=60).text
    rel = sorted(set(re.findall(rf'href="(\d{{4}}\.\d+\.{species})/"', html)))
    if not rel:
        raise RuntimeError(f"No MSigDB release found for {species}")
    return rel[-1]


def latest_wikipathways(species_name: str) -> str:
    html = requests.get(WP_BASE, timeout=60).text
    files = re.findall(r'href="(wikipathways-\d+-gmt-' + species_name.replace(" ", "_") + r'\.gmt)"', html)
    if not files:
        raise RuntimeError(f"No WikiPathways GMT for {species_name}")
    return sorted(files)[-1]


# ───────────────────────────── naming helpers ───────────────────────────────

def norm_name(s: str) -> str:
    """Lower-case alphanumeric name used for cross-tool name matching."""
    s = re.sub(r"[^0-9a-z]+", " ", str(s).lower())
    return re.sub(r"\s+", " ", s).strip()


def pretty_msigdb(name: str, prefix: str) -> str:
    s = name[len(prefix):] if name.startswith(prefix) else name
    s = s.replace("_", " ").lower()
    return s[:1].upper() + s[1:]


# ───────────────────────────── gene tables ──────────────────────────────────

def mouse_gene_table(force: bool) -> tuple[pd.DataFrame, dict]:
    mrk = pd.read_csv(fetch(MGI + "MRK_List2.rpt", "MRK_List2.rpt", force), sep="\t", dtype=str)
    mrk = mrk[(mrk["Status"] == "O") & (mrk["Marker Type"] == "Gene")]
    genes = pd.DataFrame({
        "symbol": mrk["Marker Symbol"].values,
        "mgi": mrk["MGI Accession ID"].values,
        "synonyms": mrk["Marker Synonyms (pipe-separated)"].fillna("").values,
        "feature_type": mrk["Feature Type"].fillna("").values,
    })
    ez = pd.read_csv(fetch(MGI + "MGI_EntrezGene.rpt", "MGI_EntrezGene.rpt", force), sep="\t",
                     header=None, dtype=str)
    ez = ez[ez[2] == "O"][[0, 8]].rename(columns={0: "mgi", 8: "entrez"}).dropna()
    ens = pd.read_csv(fetch(MGI + "MRK_ENSEMBL.rpt", "MRK_ENSEMBL.rpt", force), sep="\t",
                      header=None, dtype=str)[[0, 5]].rename(columns={0: "mgi", 5: "ensembl"})
    genes = genes.merge(ez, on="mgi", how="left").merge(ens, on="mgi", how="left")

    hom = pd.read_csv(fetch(MGI + "HOM_MouseHumanSequence.rpt", "HOM_MouseHumanSequence.rpt", force),
                      sep="\t", dtype=str)
    m = hom[hom["NCBI Taxon ID"] == "10090"][["DB Class Key", "Symbol"]]
    h = hom[hom["NCBI Taxon ID"] == "9606"][["DB Class Key", "Symbol"]]
    pairs = m.merge(h, on="DB Class Key", suffixes=("_mouse", "_human"))
    m2h = pairs.groupby("Symbol_mouse")["Symbol_human"].apply(lambda x: "|".join(sorted(set(x))))
    genes["ortholog"] = genes["symbol"].map(m2h).fillna("")
    genes = genes.drop_duplicates("symbol").reset_index(drop=True)
    orth = {"mouse_to_human": m2h.to_dict(),
            "human_to_mouse": pairs.groupby("Symbol_human")["Symbol_mouse"]
            .apply(lambda x: "|".join(sorted(set(x)))).to_dict()}
    return genes, orth


def human_gene_table(force: bool, orth: dict) -> pd.DataFrame:
    hg = pd.read_csv(fetch(HGNC_URL, "hgnc_complete_set.txt", force), sep="\t", dtype=str,
                     low_memory=False)
    hg = hg[hg["status"] == "Approved"]
    syn = (hg["alias_symbol"].fillna("").str.strip('"') + "|" + hg["prev_symbol"].fillna("").str.strip('"'))
    genes = pd.DataFrame({
        "symbol": hg["symbol"].values,
        "hgnc": hg["hgnc_id"].values,
        "synonyms": syn.str.strip("|").values,
        "feature_type": hg["locus_type"].fillna("").values,
        "entrez": hg["entrez_id"].values,
        "ensembl": hg["ensembl_gene_id"].values,
    })
    genes["ortholog"] = genes["symbol"].map(orth.get("human_to_mouse", {})).fillna("")
    return genes.drop_duplicates("symbol").reset_index(drop=True)


# ───────────────────────────── gene-set sources ─────────────────────────────

def go_metadata(force: bool) -> dict:
    """GO id -> (name, namespace, depth) from go-basic.obo."""
    path = fetch(GO_OBO_URL, "go-basic.obo", force)
    terms, cur = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = {"is_a": []}
        elif line.startswith("[") and line.endswith("]"):
            cur = None
        elif cur is not None and ": " in line:
            k, v = line.split(": ", 1)
            if k == "id":
                cur["id"] = v
                terms[v] = cur
            elif k in ("name", "namespace"):
                cur[k] = v
            elif k == "is_a":
                cur["is_a"].append(v.split(" ! ")[0])
    depth_cache: dict[str, int] = {}

    def depth(t: str, seen=()) -> int:
        if t in depth_cache:
            return depth_cache[t]
        parents = [p for p in terms.get(t, {}).get("is_a", []) if p in terms and p not in seen]
        d = 0 if not parents else 1 + min(depth(p, seen + (t,)) for p in parents)
        depth_cache[t] = d
        return d

    return {t: (v.get("name", t), v.get("namespace", ""), depth(t)) for t, v in terms.items()}


def msigdb_sets(release: str, fname: str, force: bool) -> dict:
    url = f"{MSIGDB_BASE}{release}/{fname}"
    return json.load(open(fetch(url, fname, force)))


def build_sets(organism: str, genes: pd.DataFrame, orth: dict, force: bool) -> tuple[list, dict]:
    spec = ORGANISMS[organism]
    sp = spec["msigdb_species"]
    rel = latest_msigdb_release(sp)
    ver = rel.split(".")[0] + "." + rel.split(".")[1]
    v = f"v{rel}"
    pfx = "m" if sp == "Mm" else "c"
    records: list[dict] = []
    manifest: dict = {"msigdb_release": rel}
    entrez2sym = dict(zip(genes["entrez"].dropna(), genes.loc[genes["entrez"].notna(), "symbol"]))

    def add(db, set_id, term_key, name, genes_, url="", depth=None, source=""):
        g = sorted(set(x for x in genes_ if x))
        if len(g) >= MIN_SET_SIZE:
            records.append({"set_id": set_id, "term_key": term_key, "database": db, "name": name,
                            "name_norm": norm_name(name), "genes": g, "url": url,
                            "go_depth": depth, "source": source})

    # GO (MSigDB C5/M5 GO, propagated; ids from exactSource, names/depth from GO OBO)
    go = go_metadata(force)
    for sub, db in (("bp", "GO_BP"), ("mf", "GO_MF"), ("cc", "GO_CC")):
        d = msigdb_sets(rel, f"{pfx}5.go.{sub}.{v}.json", force)
        for nm, info in d.items():
            gid = info.get("exactSource", "")
            name, _, dep = go.get(gid, (pretty_msigdb(nm, nm.split("_")[0] + "_"), "", None))
            add(db, gid or nm, gid or nm, name, info["geneSymbols"],
                f"https://amigo.geneontology.org/amigo/term/{gid}", dep, f"MSigDB {rel}")
    manifest["GO"] = f"MSigDB {rel} ({pfx.upper()}5 GO); names: go-basic.obo {date.today()}"

    # Reactome (MSigDB CP:REACTOME; names from Reactome)
    rnames = pd.read_csv(fetch(REACTOME_NAMES_URL, "ReactomePathways.txt", force), sep="\t",
                         header=None, names=["stid", "name", "species"], dtype=str)
    rname = dict(zip(rnames["stid"], rnames["name"]))
    d = msigdb_sets(rel, f"{pfx}2.cp.reactome.{v}.json", force)
    for nm, info in d.items():
        st = info.get("exactSource", "")
        num = st.split("-")[-1] if st else nm
        add("REACTOME", f"REACT:{st or nm}", f"REACT:{num}", rname.get(st, pretty_msigdb(nm, "REACTOME_")),
            info["geneSymbols"], f"https://reactome.org/content/detail/{st}", source=f"MSigDB {rel}")
    manifest["REACTOME"] = f"MSigDB {rel} CP:REACTOME"

    # Hallmark
    d = msigdb_sets(rel, f"{'mh' if sp == 'Mm' else 'h'}.all.{v}.json", force)
    for nm, info in d.items():
        add("HALLMARK", f"HALLMARK:{nm}", f"HALLMARK:{norm_name(nm.replace('HALLMARK_', ''))}",
            pretty_msigdb(nm, "HALLMARK_"), info["geneSymbols"],
            f"https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName={nm}", source=f"MSigDB {rel}")
    manifest["HALLMARK"] = f"MSigDB {rel}"

    # KEGG (KEGG REST; academic use — see KEGG licence)
    k = spec["kegg"]
    names = pd.read_csv(fetch(f"{KEGG_REST}list/pathway/{k}", f"kegg_list_{k}.tsv", force), sep="\t",
                        header=None, names=["pid", "name"], dtype=str)
    links = pd.read_csv(fetch(f"{KEGG_REST}link/{k}/pathway", f"kegg_link_{k}.tsv", force), sep="\t",
                        header=None, names=["pid", "gene"], dtype=str)
    links["pid"] = links["pid"].str.replace("path:", "", regex=False)
    links["sym"] = links["gene"].str.split(":").str[1].map(entrez2sym)
    kname = {p: re.sub(r" - [A-Z][a-z]+ [a-z]+ \(.*\)$", "", n) for p, n in zip(names["pid"], names["name"])}
    for pid, grp in links.dropna().groupby("pid"):
        add("KEGG", f"KEGG:{pid}", f"KEGG:{pid[-5:]}", kname.get(pid, pid), grp["sym"].tolist(),
            f"https://www.kegg.jp/pathway/{pid}", source=f"KEGG REST {date.today()}")
    manifest["KEGG"] = f"KEGG REST API, downloaded {date.today()} (academic-use terms apply)"

    # WikiPathways (Entrez ids -> symbols)
    wp_file = latest_wikipathways(spec["reactome"])
    for line in open(fetch(WP_BASE + wp_file, wp_file, force), encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        parts = f[0].split("%")
        wid = parts[2] if len(parts) > 2 else parts[0]
        add("WIKIPATHWAYS", f"WP:{wid}", f"WP:{wid}", parts[0], [entrez2sym.get(x) for x in f[2:]],
            f[1], source=wp_file)
    manifest["WIKIPATHWAYS"] = wp_file

    # DisGeNET (Enrichr snapshot, human symbols; orthologs for mouse)
    lib = fetch(ENRICHR_LIB + "DisGeNET", "enrichr_DisGeNET.txt", force)
    h2m = orth.get("human_to_mouse", {})
    for line in open(lib, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        hs = [x.split(",")[0] for x in f[2:] if x]
        gs = hs if organism == "hsapiens" else [m for x in hs for m in h2m.get(x, "").split("|")]
        add("DISGENET", f"DISGENET:{f[0]}", f"DISGENET:{norm_name(f[0])}", f[0], gs,
            source="Enrichr DisGeNET snapshot")
    manifest["DISGENET"] = "Enrichr 'DisGeNET' library snapshot" + (
        " (human->mouse via MGI orthologs)" if organism != "hsapiens" else "")
    manifest["msigdb_version"] = ver
    return records, manifest


# ───────────────────────────── compile ──────────────────────────────────────

def compile_library(organism: str, force: bool = False) -> Path:
    out = library_path(organism)
    (out / "gmt").mkdir(parents=True, exist_ok=True)
    print(f"[library] building {organism} …")
    t0 = time.time()
    mouse_genes, orth = mouse_gene_table(force)
    genes = mouse_genes if organism == "mmusculus" else human_gene_table(force, orth)
    records, manifest = build_sets(organism, genes, orth, force)

    known = set(genes["symbol"])
    set_genes = sorted({g for r in records for g in r["genes"]})
    gene_index = sorted(known | set(set_genes))
    gpos = {g: i for i, g in enumerate(gene_index)}
    rows, cols = [], []
    for i, r in enumerate(records):
        for g in r["genes"]:
            rows.append(i)
            cols.append(gpos[g])
    mat = sparse.csr_matrix((np.ones(len(rows), dtype=np.bool_), (rows, cols)),
                            shape=(len(records), len(gene_index)))
    sparse.save_npz(out / "membership.npz", mat)
    (out / "gene_index.txt").write_text("\n".join(gene_index))

    sets = pd.DataFrame([{k: v for k, v in r.items() if k != "genes"} for r in records])
    sets["size"] = [len(r["genes"]) for r in records]
    sets.to_parquet(out / "sets.parquet", index=False)
    genes.to_parquet(out / "genes.parquet", index=False)

    by_db = defaultdict(list)
    for r in records:
        by_db[r["database"]].append(r)
    for db, recs in by_db.items():
        with open(out / "gmt" / f"{db}.gmt", "w") as fh:
            for r in recs:
                fh.write("\t".join([r["set_id"], r["name"].replace("\t", " ")] + r["genes"]) + "\n")

    manifest.update({
        "organism": organism, "built": str(date.today()),
        "n_sets": {db: len(v) for db, v in by_db.items()},
        "n_genes_in_sets": len(set_genes), "n_genes_table": int(len(genes)),
        "build_seconds": round(time.time() - t0, 1),
    })
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[library] {organism}: {len(records)} sets, {len(set_genes)} genes "
          f"({time.time() - t0:.0f}s) -> {out}")
    for db, n in manifest["n_sets"].items():
        print(f"    {db:<14} {n:>6} sets")
    return out
