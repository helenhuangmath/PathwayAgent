# Tools and presets

`pathwayagent tools` prints this list for your installation.

## Tools

Each tool belongs to a method family. The consensus counts families as well as tools, so eight
over-representation tools do not count as eight independent votes.

### Over-representation (Fisher / hypergeometric)

| Tool | Runs | Databases |
|---|---|---|
| `local_ora` | local | GO BP/MF/CC, KEGG, Reactome, WikiPathways, Hallmark, DisGeNET |
| `enrichr` | web | GO BP/MF/CC, KEGG, Reactome, WikiPathways, Hallmark, DisGeNET |
| `gprofiler` | web | GO BP/MF/CC, KEGG, Reactome, WikiPathways |
| `reactome` | web | Reactome |
| `string` | web | GO BP/MF/CC, KEGG, Reactome, WikiPathways |
| `clusterprofiler` | local, R | same gene sets as `local_ora` |
| `clusterprofiler_go` | local, R | GO BP/MF/CC from org.db |
| `reactomepa` | local, R | Reactome |
| `metascape` | imported result | the ontologies present in the package |

`clusterprofiler` in shared-library mode reproduces `local_ora` exactly (this is tested), so it is
reported as a cross-check and counted once.

### Pre-ranked GSEA

| Tool | Runs | Databases |
|---|---|---|
| `local_gsea` (blitzGSEA) | local | GO BP/MF/CC, KEGG, Reactome, WikiPathways, Hallmark |
| `fgsea` (fgseaMultilevel) | local, R | GO BP/MF/CC, KEGG, Reactome, WikiPathways, Hallmark |

These need a genome-wide ranking of at least 1,000 genes.

### Rank tests

| Tool | Runs | Databases |
|---|---|---|
| `local_mhg` (minimum hypergeometric, ordered) | local | all eight |
| `local_rank` (Wilcoxon rank-sum) | local | all eight |
| `gprofiler_ordered` (g:Profiler ordered query) | web | GO BP/MF/CC, KEGG, Reactome, WikiPathways |

These need scores. Unranked lists skip them.

## Presets

| Preset | Tools | Typical time |
|---|---|---|
| `quick` | `local_ora`, `local_mhg`, `local_gsea`, `local_rank` | a few seconds |
| `standard` (default) | `quick` + `enrichr`, `gprofiler`, `reactome`, `string`, `clusterprofiler`, `fgsea` | |
| `everything` | `standard` + `gprofiler_ordered`, `clusterprofiler_go`, `reactomepa`, and `metascape` when a result is supplied | about 150 s for three lists in the tutorial |

Choose individual tools with `--tools local_ora,enrichr,gprofiler`. With `--local-only`, tools that
send genes to a web service are skipped whatever the preset.

A tool that cannot run (no scores, no R, service unavailable, deadline reached) is reported with
status `skipped`, `error` or `timeout`, and the run continues with the others.

## Metascape

Metascape has no public API, so its result is imported:

1. Analyse the lists at [metascape.org](https://metascape.org).
2. On the report page choose **All in One Zip File**.
3. Pass the zip, or the unpacked folder, with `--metascape PATH`.

A multi-list package serves every list it contains, matched by list name. For one package per list
use `--metascape TEX_G1=g1.zip TEX_G2=g2.zip`. Only the ontologies present in the package count as
covered by Metascape. Without a package the tool is skipped, unless `PATHWAYAGENT_MSBIO` points at
an MSBio installation whose Docker container is running (this path has not been tested).

The [tutorial](tutorial.md#6-adding-metascape) shows how Metascape compares with the other tools.

## Not yet implemented

| Tool | Reason |
|---|---|
| WebGestalt | needs the WebGestaltR package |
| DAVID | web service requires registration |
| GREAT | needs genomic regions (BED) as input |
