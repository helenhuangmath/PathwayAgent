# Quick start

A first analysis in three steps. For details on each step see [Installation](INSTALL.md) and
[Inputs and command line](usage.md).

## 1. Install

```bash
git clone https://github.com/helenhuangmath/PathwayAgent.git
cd PathwayAgent
python -m venv .venv
.venv/bin/pip install -e ".[dev]" scikit-learn
```

Python 3.11 or newer is required. `scikit-learn` is optional; it is used to build the term-landscape
maps behind figure 08.

## 2. Build the gene-set library (one time)

```bash
.venv/bin/python -m pathwayagent db build
```

This downloads the reference gene sets (GO, KEGG, Reactome, WikiPathways, MSigDB Hallmark, DisGeNET)
for mouse and human and takes about five minutes. Add `--organism mmusculus` to build one organism
only. Check the result with:

```bash
.venv/bin/python -m pathwayagent db list
```

## 3. Run

```bash
.venv/bin/python -m pathwayagent run testing/gene_lists/genels_ex_module123.txt --out results/modules
```

The example file has three gene modules, one per column. Progress is printed per tool, and the run
ends with the path of the report:

```
Report: results/modules/report.html
Summary: results/modules/summary.txt
```

Open `report.html` in a browser. It has one tab per gene list, sortable tables, a network you can
hover, and a link to the PDF of every figure.

## Common variations

```bash
# local tools only, a few seconds, nothing sent to web services
.venv/bin/python -m pathwayagent run genes.txt --preset quick --local-only

# a DE table with a background of expressed genes, every tool
.venv/bin/python -m pathwayagent run de_table.csv --universe expressed_genes.txt --preset everything

# only some columns of a multi-column file
.venv/bin/python -m pathwayagent run lists.csv --lists TEX_G1 TEX_G2

# the same analysis in the browser, at http://127.0.0.1:8765
.venv/bin/python -m pathwayagent serve
```

!!! tip "Supply a background"
    For real experiments pass the genes that were measured or expressed with `--universe`. Without
    it every annotated gene is used as background, which inflates tissue-specific terms. The
    [tutorial](tutorial.md) shows what this looks like in practice.

## Without R

The R tools (`clusterprofiler`, `clusterprofiler_go`, `reactomepa`, `fgsea`) need R with
clusterProfiler, fgsea, ReactomePA and org.Mm.eg.db / org.Hs.eg.db. If R is missing those tools are
reported as skipped and the rest of the run continues. Point to a specific R with
`PATHWAYAGENT_RSCRIPT=/path/to/Rscript`.

## Next

- [Tutorial](tutorial.md): a complete run on real data and how to read the outputs.
- [Tools and presets](tools.md): what each preset runs.
- [Python API](PYTHON_API.md): run the analysis from a script or notebook.
