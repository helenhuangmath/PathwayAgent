# Inputs and command line

## Input formats

The format is detected automatically. Files can be plain text, CSV, TSV or Excel.

| Input | Example | Analysed as |
|---|---|---|
| One gene per line | `Stat1`<br>`Irf7`<br>`Isg15` | Gene list (over-representation) |
| Gene + score | `Stat1,8.0`<br>`Irf7,6.0` | Ranked list |
| DE table | columns `gene`, `log2FoldChange`, `pvalue` | Ranked list, by the statistic column when there is one, otherwise by signed −log10 p |
| Several columns under a header row | `TEX_G1,TEX_G2,TEX_G5G6` | One gene list per column |

A header row is recognised by a header word (gene, module, cluster, list, set, up, down) or by column
names built on a shared stem (`TEX_G1`, `TEX_G2`). Check the parsed lists when your headers are
unusual; `id_mapping.tsv` in the output shows what each input gene was mapped to.

### Ranked or unranked?

Whether the list carries scores decides which tools can run:

- **Unranked list**: over-representation tools only. The rank-based tools are skipped.
- **Ranked, genome-wide (at least 1,000 genes)**: over-representation, GSEA (`local_gsea`, `fgsea`)
  and rank tests.
- **Ranked, short** (for example a module with its Wald statistics): GSEA is skipped and the ordered
  test (`local_mhg`) is used instead.

Higher scores mean stronger positive regulation.

## Commands

```
pathwayagent run FILE [FILE ...]   analyse one or more gene-list files
pathwayagent db build | list       build or list the local gene-set library
pathwayagent serve                 start the local web app
pathwayagent tools                 list available tools and presets
```

`python -m pathwayagent ...` is equivalent to the `pathwayagent` command.

### `run`

```bash
pathwayagent run lists.csv --lists TEX_G1 TEX_G2 --preset everything \
    --universe expressed_genes.txt --out results/my_run --title "My modules"
```

| Option | Default | Meaning |
|---|---|---|
| `FILE ...` | required | One or more input files. Multi-column files become one list per column. |
| `--organism` | `auto` | `auto`, `mmusculus` or `hsapiens`. |
| `--preset` | `standard` | `quick`, `standard` or `everything`. See [Tools and presets](tools.md). |
| `--tools` | | Comma-separated tool names. Overrides `--preset`. |
| `--db` | all eight | Comma-separated databases: `GO_BP`, `GO_MF`, `GO_CC`, `KEGG`, `REACTOME`, `WIKIPATHWAYS`, `HALLMARK`, `DISGENET`. |
| `--lists` | every list | Analyse only these list names from a multi-column file. |
| `--universe` | all annotated genes | File with the background genes. Recommended. |
| `--profile` | | Project profile YAML (genes and pathways of interest, background). |
| `--out` | `results/<timestamp>` | Output directory. Use a fresh one per run. |
| `--min-size`, `--max-size` | `10`, `500` | Gene-set size bounds. |
| `--alpha` | `0.05` | Adjusted P threshold. |
| `--deadline` | `60` | Time budget in seconds for the tool stage. Tools still running are marked as timed out. |
| `--local-only` | off | Never send genes to external services. |
| `--metascape` | | `PATH` or `LIST=PATH`: Metascape result package(s) to include. See [Metascape](tools.md#metascape). |
| `--title` | `PathwayAgent run` | Report title. |
| `--print-summary` | off | Print `summary.txt` when the run ends. |

### `db`

```bash
pathwayagent db build                        # mouse and human
pathwayagent db build --organism mmusculus   # one organism
pathwayagent db build --force                # re-download the sources
pathwayagent db list                         # what is built, with set counts
```

The build also computes the term-landscape maps (needs scikit-learn) and, when R is available,
caches the org.db GO annotation used by `clusterprofiler_go`.

### `serve`

```bash
pathwayagent serve --port 8765 --host 127.0.0.1 --no-browser
```

See [Web app and REST API](webapp.md).

## Where files go

| Variable | Purpose |
|---|---|
| `PATHWAYAGENT_LIBRARY` | Folder holding the built gene-set library. |
| `PATHWAYAGENT_RESULTS` | Default parent folder of run outputs. |
| `PATHWAYAGENT_DATA` | Parent of `library/` and `results/`, when both live together. |
| `PATHWAYAGENT_CACHE` | SQLite cache of web-service responses. |
| `PATHWAYAGENT_RSCRIPT` | Path to `Rscript` (default `/usr/local/bin/Rscript`). |
| `MPLCONFIGDIR` | Matplotlib cache folder; set it on servers where the default is not writable. |

Without these, the library and results go to `library/` and `results/` inside the repository (both
ignored by Git), or to a `large_data/` folder beside the repository when one exists. See
[Installation](INSTALL.md) for details.
