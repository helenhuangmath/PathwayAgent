# How the consensus works

The consensus turns the separate results of every tool into one table per gene list, with a score
and an evidence tier for each term.

## 1. Matching terms across tools

Terms are matched by ID (GO, KEGG, Reactome, WikiPathways) or, when there is no shared ID, by
normalised name. When a tool reports the same term more than once, its best result is kept.

## 2. Eligibility

A tool is **eligible** for a term when it ran and covers the term's database. Reactome's own
service cannot say anything about a Hallmark set, so it is neither for nor against it. In figure 03
an ineligible tool is drawn as a dash.

If an eligible tool did not report the term, that counts as not significant (adjusted P = 1).

## 3. Method families

Tools are grouped into families: over-representation, pre-ranked GSEA, and rank tests. The score is
averaged within each family first and then across families, so many tools of one kind do not outvote
a single tool of another kind.

## 4. The score

For each term, over its eligible tools:

| Part | Weight | Definition |
|---|---|---|
| Significance | 40% | Median over tools of −log10 adjusted P, capped at 10 and scaled to 0–1; averaged across families. |
| Support | 30% | Fraction of eligible tools with adjusted P < α; averaged across families. |
| Effect size | 15% | Median, over the reporting tools, of the term's within-tool percentile of effect size (absolute NES for GSEA). |
| Direction | 15% | Agreement in sign among significant directional results; 1 when there are none. |

```
score = 100 × specificity × (0.40 × significance + 0.30 × support + 0.15 × effect + 0.15 × direction)
```

**Specificity** penalises very large and very general terms: gene sets above 200 genes are scaled
by `1 / (1 + log10(size / 200))`, and GO terms at depth 3 or less are multiplied by 0.7.

## 5. Tiers

| Tier | Requirements |
|---|---|
| **ROBUST** | Significant in at least 2 tools; support of at least 60%; median adjusted P of 0.01 or lower; direction consistency of at least 0.8; and significant in at least 2 method families when more than one family was eligible. |
| **SUPPORTED** | Significant in at least 2 tools, with support of at least 30%. |
| **SINGLE-TOOL** | Significant in one tool only. |

Terms are ordered by tier, then by score, then by best adjusted P.

!!! warning "One method family"
    With an unranked gene list only over-representation tools run. ROBUST then means agreement
    between implementations and databases of the same test, not between independent methods.
    Supplying scores (gene + statistic, or a DE table) lets the GSEA and rank tools take part.

`clusterprofiler` in shared-library mode reproduces `local_ora` exactly, so it is reported but
counted once.

## 6. Rank consistency and stability

Two further columns help judge a term:

- `rra_p`, `rra_q`: a Robust Rank Aggregation p-value (Kolde et al. 2012) over the term's
  normalised rank in each eligible tool.
- `loo_rank_range`, `depends_on`: the consensus is recomputed leaving out one tool at a time. The
  first column gives the range of ranks the term took; the second names any tool whose removal makes
  the term lose significance.

## 7. Themes

Significant terms from all databases are grouped by the genes they hit: Cohen's kappa between the
terms' hit-gene vectors, average-linkage hierarchical clustering, cut at kappa 0.3 (the approach
used by DAVID and Metascape). A theme is scored by its best term, with a 10% bonus for each
additional database that contributes a term. Themes are numbered by score, so T1 is the best.

Themes are derived again for each gene list: T1 of one list is unrelated to T1 of another.

## 8. Key genes

A gene's weight is the sum, over ROBUST and SUPPORTED themes, of the theme score multiplied by the
fraction of the theme's significant results that contain the gene. When scores were supplied, genes
with extreme scores are weighted up.

A gene is flagged as a **driver** when removing it from the list makes the theme's best locally
tested term lose significance.

## Columns of `consensus_terms.tsv`

| Column | Meaning |
|---|---|
| `rank`, `term_key`, `term_name`, `database`, `term_size` | Identity of the term. |
| `consensus_score`, `tier`, `tier_label` | Score (0–100) and tier. |
| `tool_support`, `n_tools_sig`, `n_tools_eligible` | Tools significant out of tools eligible. |
| `n_families_sig`, `n_families_eligible` | The same, for method families. |
| `S_sig`, `S_sup`, `S_eff`, `direction_consistency`, `specificity` | The parts of the score. |
| `direction` | Sign of the significant directional results (+1, −1, or 0). |
| `best_padj`, `median_padj`, `padj__<tool>` | Adjusted P: best, median over eligible tools, and per tool. |
| `rra_p`, `rra_q`, `loo_rank_range`, `depends_on` | Rank consistency and leave-one-tool-out stability. |
| `tools_significant`, `hit_genes`, `n_hit_genes` | Which tools, and which of your genes. |
