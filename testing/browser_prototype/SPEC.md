# Pathway Quick-Run — SPEC (version 1)

## Purpose
For a wet-lab / computational biologist who has one or more gene lists (e.g. up- and
down-regulated genes from DESeq2) and wants a fast, Enrichr-style answer to:
**"Which known pathways are over-represented in each of my gene lists, and how do the lists compare?"**
Runs offline by double-clicking `index.html` (macOS and Windows).

## What it does
1. Takes one gene list, or several named groups, by pasting into a box or loading a `.txt` / `.csv` file.
2. Species choice: **Human** or **Mouse** (each uses that species' own MSigDB gene sets).
3. Checkboxes for four bundled databases: **Hallmark, GO Biological Process, Reactome, WikiPathways**.
4. Before running, lists problems per group: Excel-mangled names (`1-Sep`, `Mar-1`, …),
   unrecognised symbols, duplicates. It never silently "fixes" a name.
5. **Run** button: for each group × database, an over-representation test (one-sided
   hypergeometric = Fisher's exact test, as in Enrichr / clusterProfiler), then
   Benjamini–Hochberg FDR within each group × database.
6. Results table per group: Database, Pathway, Overlap (k/K), Genes in overlap, P-value, FDR.
   Sorted by P-value; rows with FDR < 0.05 highlighted. Download as CSV.
7. One **dot plot** comparing groups: columns = groups, rows = top 5 pathways per group
   (by P-value), dot size = overlap count, colour = −log10(FDR). Download as **SVG** and **PNG**.
8. Results appear within ~2 seconds for 3 lists × 4 databases on a normal laptop.

## What it doesn't do (v1)
- No pathway network, bar plots, cell-type/tissue fields (planned v2).
- No GSEA / ranked lists, no fold-changes, no KEGG (licence), no other software (IPA, clusterProfiler).
- No internet, no AI model, no saving of sessions. No automatic correction of gene names.

## Inputs
**A. Pasted text or `.txt` file.** Genes separated by new lines, commas, tabs or spaces.
A line starting with `>` starts a new named group. With no `>` line, it is one group called "List1".
Real example: `example_Tex_Tmem.txt` in this folder (mouse, 20 literature markers per group), start of file:
```
>Tex
Tox
Pdcd1
Havcr2
...
>Tmem
Tcf7
Lef1
Sell
...
```
**B. `.csv` file (Metascape style)**: one column per group, first row = group names.
```
Tex,Tmem
Tox,Tcf7
Pdcd1,Lef1
Havcr2,Sell
```
**C. Bundled gene-set files (GMT, from MSigDB, CC BY 4.0).** One pathway per line:
name, description/URL, then genes, tab-separated (example shortened):
```
HALLMARK_INTERFERON_GAMMA_RESPONSE	https://www.gsea-msigdb.org/...	ADAR	APOL6	B2M	CD274	CXCL10	...
```
Gene matching ignores upper/lower case (`tox` = `TOX`); the table reports the database's spelling.

## Method (used by the test)
For one group and one database: **N** = number of distinct genes in that database (the "background");
**n** = number of the group's genes that are in N (unrecognised genes are left out);
for each pathway, **K** = its size, **k** = how many of the group's genes are in it.
P = probability of seeing k or more by chance = Σ_{x≥k} C(K,x)·C(N−K,n−x) / C(N,n).
Only pathways with k ≥ 1 are tested and shown. FDR = Benjamini–Hochberg over those pathways.

## Tiny test case (toy database, human, checkable by hand)
Toy database "TOY", 3 pathways, background N = 10 genes:
- EFFECTOR: GZMB, PRF1, IFNG (K=3)
- EXHAUSTION: TOX, PDCD1, HAVCR2, CD8A (K=4)
- MEMORY: TCF7, SELL, IL7R (K=3)

Input:
```
>G1
GZMB
PRF1
IFNG
TOX
1-Sep
FAKE1
gzmb
>G2
TCF7, SELL, IL7R, TOX, PDCD1
```
All-pairs count C(10,4) = 210, C(10,5) = 252.

**Expected answers**
| # | Group | Check / pathway | k/K | P-value | FDR |
|---|---|---|---|---|---|
| 1 | G1 | Warnings: `1-Sep` Excel-mangled; `FAKE1` unrecognised; `gzmb` duplicate of GZMB → n = 4 | | | |
| 2 | G1 | EFFECTOR | 3/3 | 7/210 = **0.03333** | **0.06667** |
| 3 | G1 | EXHAUSTION | 1/4 | 195/210 = **0.9286** | **0.9286** |
| 4 | G1 | MEMORY | 0/3 | not shown | |
| 5 | G2 | No warnings, n = 5 | | | |
| 6 | G2 | MEMORY | 3/3 | 21/252 = **0.08333** | **0.1667** |
| 7 | G2 | EXHAUSTION | 2/4 | 186/252 = **0.7381** | **0.7381** |
| 8 | G2 | EFFECTOR | 0/3 | not shown | |

Numbers must match to 4 significant figures.

## Open decisions (suggested default in bold)
1. Background N: **all genes in the chosen database** (clusterProfiler style) vs. a fixed ~20,000 genes (Enrichr).
2. "Recognised gene" = **appears in any bundled database for that species** (no extra gene list to download).
3. Pathway size limits: **none** (Enrichr style); option to add a 5–500 filter later.
4. MSigDB version: **latest release available at build time**, recorded in the app footer.
5. Dot plot: **top 5 per group**, dot size = overlap count.
6. Mouse: **MSigDB mouse collections (MH, M5 GO:BP, M2 Reactome, M2 WikiPathways)**, not converted human sets.

## Done when
- [ ] Double-clicking `index.html` opens the app in Chrome/Safari/Firefox with Wi-Fi turned off.
- [ ] Double-clicking `test.html` shows every row of the tiny test case as PASS.
- [ ] Loading `example_Tex_Tmem.txt`, choosing Mouse and clicking Run shows results in about 2 seconds.
- [ ] Tex hits exhaustion / immune-checkpoint-type pathways; Tmem hits memory / IL-7 / T cell differentiation-type pathways.
- [ ] `1-Sep` and a made-up gene show up as warnings before results.
- [ ] CSV, SVG and PNG downloads open in Excel / Illustrator / Preview.
- [ ] The footer shows the MSigDB version and date of the bundled databases.
