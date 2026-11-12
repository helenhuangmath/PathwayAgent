"""Convert the MSigDB .gmt files in data_source/ into data/*.js files the web page can load offline.
Run once when updating databases:  python3 tools/make_data.py
(The app itself never needs this script.)"""
import json, pathlib
root = pathlib.Path(__file__).resolve().parent.parent
SRC = root / "data_source"; OUT = root / "data"; OUT.mkdir(exist_ok=True)
FILES = {
    ("Hs", "Hallmark"): "h.all.v2026.1.Hs.symbols.gmt",
    ("Hs", "GO_BP"): "c5.go.bp.v2026.1.Hs.symbols.gmt",
    ("Hs", "Reactome"): "c2.cp.reactome.v2026.1.Hs.symbols.gmt",
    ("Hs", "WikiPathways"): "c2.cp.wikipathways.v2026.1.Hs.symbols.gmt",
    ("Mm", "Hallmark"): "mh.all.v2026.1.Mm.symbols.gmt",
    ("Mm", "GO_BP"): "m5.go.bp.v2026.1.Mm.symbols.gmt",
    ("Mm", "Reactome"): "m2.cp.reactome.v2026.1.Mm.symbols.gmt",
    ("Mm", "WikiPathways"): "m2.cp.wikipathways.v2026.1.Mm.symbols.gmt",
}
for (species, db), fname in FILES.items():
    lines = []
    for line in (SRC / fname).read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            lines.append("\t".join([parts[0]] + parts[2:]))  # drop the URL column to save space
    info = {"species": species, "db": db, "source": fname, "version": "MSigDB " + fname.split(".v")[1].split(".symbols")[0],
            "downloaded": "2026-09-24", "gmt": "\n".join(lines)}
    (OUT / f"{species}_{db}.js").write_text("PQ.registerLibrary(" + json.dumps(info) + ");\n")
    print(species, db, len(lines), "pathways")
