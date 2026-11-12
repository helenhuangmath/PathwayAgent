"""Write every output of a run (design §10.3): tables, figures, text, Excel, HTML, manifest, methods."""
from __future__ import annotations

import json
import platform
from datetime import datetime
from pathlib import Path

import pandas as pd

from .. import __version__
from ..config import DATABASES, FAMILIES
from ..viz import static as viz
from .html import render_html

CONS_COLS = ["rank", "term_name", "database", "term_key", "tier_label", "consensus_score", "tool_support",
             "n_families_sig", "best_padj", "median_padj", "rra_q", "direction", "term_size", "n_hit_genes",
             "hit_genes", "tools_significant", "loo_rank_range", "depends_on", "theme_id", "S_sig", "S_sup",
             "S_eff", "direction_consistency", "specificity", "url"]


def _emit(cb, **ev):
    if cb:
        try:
            cb(ev)
        except Exception:
            pass


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)


def write_query(qr, run, qdir: Path) -> dict:
    qdir.mkdir(parents=True, exist_ok=True)
    (qdir / "tools").mkdir(exist_ok=True)
    figs = {}
    q = qr.query
    status = pd.DataFrame([r.summary() for r in qr.tool_results])
    status["databases"] = status["databases"].apply(lambda x: ",".join(x) if isinstance(x, list) else x)
    status.to_csv(qdir / "tool_status.tsv", sep="\t", index=False)
    if q.mapping is not None:
        q.mapping.to_csv(qdir / "id_mapping.tsv", sep="\t", index=False)
    if len(qr.records):
        qr.records.to_parquet(qdir / "all_results.parquet", index=False)
        qr.records.to_csv(qdir / "all_results.tsv.gz", sep="\t", index=False)
    cons = qr.consensus
    if len(cons):
        cols = [c for c in CONS_COLS if c in cons.columns] + [c for c in cons.columns if c.startswith("padj__")]
        cons.loc[cons["tier"] > 0, cols].to_csv(qdir / "consensus_terms.tsv", sep="\t", index=False)
        cons.to_parquet(qdir / "consensus_all.parquet", index=False)
    if len(qr.themes):
        qr.themes.to_csv(qdir / "themes.tsv", sep="\t", index=False)
    if qr.key_genes is not None and len(qr.key_genes):
        qr.key_genes.to_csv(qdir / "key_genes.tsv", sep="\t", index=False)

    # per-tool tables + figures (goal 1)
    for r in qr.tool_results:
        if r.status != "ok" or not len(r.records):
            continue
        rec = r.records.sort_values("padj")
        rec.to_csv(qdir / "tools" / f"{_safe(r.tool)}.tsv", sep="\t", index=False)
        f1 = viz.tool_dotplot(rec, f"{q.name}: {viz.tool_label(r.tool)} top terms", qdir / "figures" / "tools" /
                              f"{_safe(r.tool)}_dotplot")
        f2 = viz.tool_network(rec, f"{q.name}: {viz.tool_label(r.tool)} enrichment map", qdir / "figures" / "tools" /
                              f"{_safe(r.tool)}_network")
        figs[f"tool:{r.tool}"] = [x for x in (f1, f2) if x]

    from ..pipeline import REPLICATES
    ran = [r.tool for r in qr.tool_results if r.status in ("ok", "empty")]
    # one tool set, one order, in every figure: replicates (e.g. clusterProfiler == local ORA) left out
    tools_ran = [t for t in ran if not (t in REPLICATES and REPLICATES[t] in ran)]
    tool_dbs = {r.tool: r.databases for r in qr.tool_results if r.status in ("ok", "empty")}
    F = qdir / "figures"
    th = qr.themes
    if len(cons):
        figs["consensus_bar"] = viz.consensus_bar(cons, th, f"{q.name}: consensus pathways by theme",
                                                  F / "01_consensus_bar")
        png, cy = viz.theme_network(cons, th, f"{q.name}: enrichment map of consensus themes",
                                    F / "02_theme_network")
        figs["theme_network"] = png
        if cy:
            (qdir / "theme_network.cyjs").write_text(json.dumps({"elements": cy}))
            figs["_cy"] = cy
        figs["evidence_bubble"] = viz.evidence_bubble(cons, th, tools_ran, tool_dbs,
                                                      f"{q.name}: evidence per tool", F / "03_evidence_bubble")
        figs["theme_tool"] = viz.theme_tool_squares(qr.records, th, tools_ran,
                                                    f"{q.name}: theme support across tools", F / "04_theme_tool")
        figs["key_gene_bubble"] = viz.key_gene_bubble(qr.key_genes, th, f"{q.name}: key genes",
                                                      F / "05_key_gene_bubble")
        figs["gene_theme"] = viz.gene_theme_heatmap(qr.key_genes, th, f"{q.name}: key genes × themes",
                                                    F / "06_gene_theme_heatmap")
        figs["tool_overlap"] = viz.tool_overlap(qr.records, tools_ran, f"{q.name}: agreement between tools",
                                                F / "07_tool_overlap")
        ms = next((r for r in qr.tool_results if r.tool == "metascape" and r.status == "ok"), None)
        if ms is not None:
            figs["metascape_clusters"] = viz.metascape_clusters(
                cons, th, ms.extra.get("clusters", []), f"{q.name}: Metascape clusters vs consensus themes",
                F / "09_metascape_clusters")
        from ..db.landscape import load_landscape
        sigdb = cons[cons["tier"] > 0]["database"].value_counts()
        for db in [d for d in ("GO_BP", "REACTOME", "KEGG", "WIKIPATHWAYS") if sigdb.get(d, 0) >= 3][:2]:
            L = load_landscape(run.request.organism, db)
            figs[f"landscape_{db}"] = viz.term_landscape(L, cons, db, f"{q.name}: {viz.DB_LABEL.get(db, db)} term landscape",
                                                         F / f"08_landscape_{db}")

    # Excel workbook
    with pd.ExcelWriter(qdir / "results.xlsx", engine="openpyxl") as xw:
        pd.DataFrame({"PathwayAgent": [f"query: {q.name}", f"organism: {run.request.organism}",
                                       f"universe: {run.universe_size} genes ({run.universe_source})",
                                       "sheets: consensus, themes, key_genes, tool_status, one per tool",
                                       "tiers: ROBUST > SUPPORTED > SINGLE-TOOL (see methods.txt)"]}
                     ).to_excel(xw, sheet_name="README", index=False)
        if len(cons):
            cons.loc[cons["tier"] > 0, [c for c in CONS_COLS if c in cons.columns]].to_excel(
                xw, sheet_name="consensus", index=False)
        if len(qr.themes):
            qr.themes.to_excel(xw, sheet_name="themes", index=False)
        if qr.key_genes is not None and len(qr.key_genes):
            qr.key_genes.to_excel(xw, sheet_name="key_genes", index=False)
        status.to_excel(xw, sheet_name="tool_status", index=False)
        for r in qr.tool_results:
            if r.status == "ok" and len(r.records):
                r.records.sort_values("padj").head(2000).to_excel(xw, sheet_name=_safe(r.tool)[:31], index=False)
    qr.figures = figs
    return figs


def text_summary(run) -> str:
    L = []
    rq = run.request
    L.append("PathwayAgent — pathway analysis report")
    L.append("=" * 78)
    L.append(f"Run: {run.run_id}   organism: {rq.organism}   universe: {run.universe_size} genes "
             f"({run.universe_source})")
    L.append(f"Databases: {', '.join(rq.databases)}")
    L.append(f"Tools requested: {', '.join(rq.tools)}")
    L.append(f"Min/max set size: {rq.min_size}/{rq.max_size}   alpha: {rq.alpha}   timings: {run.timings}")
    for w in run.warnings:
        L.append(f"WARNING: {w}")
    for qr in run.queries:
        q = qr.query
        m = q.mapping
        n_in = len(m) if m is not None else len(q.genes)
        L.append("")
        L.append("#" * 78)
        L.append(f"# {q.name}: {n_in} input genes, {len(q.genes)} mapped"
                 + (" (ranked list)" if q.is_ranked else ""))
        L.append("#" * 78)
        L.append("Tools:")
        for r in qr.tool_results:
            s = r.summary()
            L.append(f"  {s['tool']:<20} {s['status']:<8} {s['runtime_s']:>6.1f}s  "
                     f"{s['n_significant']:>5} significant / {s['n_tested']:>6} tested  {s['message'][:90]}")
        for w in qr.warnings:
            L.append(f"  NOTE: {w}")
        c = qr.consensus
        if len(c):
            counts = c["tier_label"].value_counts()
            L.append(f"\nConsensus: {counts.get('ROBUST', 0)} ROBUST, {counts.get('SUPPORTED', 0)} SUPPORTED, "
                     f"{counts.get('SINGLE-TOOL', 0)} SINGLE-TOOL terms")
            L.append("\nTop consensus terms:")
            for r in c[c["tier"] > 0].head(20).itertuples():
                genes = ",".join(r.hit_genes.split(";")[:10])
                L.append(f"  {r.rank:>3}. [{r.tier_label:<11}] {r.term_name[:60]} ({r.database})  score={r.consensus_score:.1f}"
                         f"  tools={r.tool_support}  best padj={r.best_padj:.2g}  genes: {genes}")
        if len(qr.themes):
            L.append("\nThemes (terms clustered by shared genes across databases):")
            for t in qr.themes.head(10).itertuples():
                L.append(f"  T{t.theme_rank:<2} [{t.tier_label:<11}] {t.theme_name[:60]}  score={t.theme_score:.1f}"
                         f"  {t.n_terms} terms from {t.databases}  core genes: {','.join(t.core_genes.split(';')[:8])}")
        kg = qr.key_genes
        if kg is not None and len(kg):
            L.append("\nKey genes for follow-up:")
            for g in kg.head(15).itertuples():
                L.append(f"  {g.rank:>3}. {g.gene:<12} score={g.gene_score:>7.1f}  {g.why[:140]}")
        L.append("\nTop 5 terms per tool:")
        for r in qr.tool_results:
            if r.status == "ok" and len(r.records):
                top = r.records.sort_values("padj").head(5)
                L.append(f"  {r.tool}:")
                for x in top.itertuples():
                    L.append(f"      {x.term_name[:70]:<70} {x.database:<12} padj={x.padj:.2g}  k={x.overlap}")
    return "\n".join(L) + "\n"


def methods_text(run) -> str:
    rq = run.request
    ran = sorted({r.tool for qr in run.queries for r in qr.tool_results if r.status in ("ok", "empty")})
    vers = {}
    for qr in run.queries:
        for r in qr.tool_results:
            if r.version:
                vers[r.tool] = r.version
    lm = run.library_manifest
    tool_versions = ", ".join(f"{t} ({vers.get(t, 'n/a')})" for t in ran)
    return (
        f"Pathway enrichment was performed with PathwayAgent v{__version__} ({datetime.now():%Y-%m-%d}). "
        f"Gene symbols were mapped to official {rq.organism} symbols using MGI/HGNC tables (synonyms resolved "
        f"when unambiguous). Gene sets ({', '.join(rq.databases)}) were compiled locally from MSigDB "
        f"{lm.get('msigdb_release', '')} (GO, Reactome, Hallmark), the KEGG REST API, WikiPathways "
        f"({lm.get('WIKIPATHWAYS', '')}) and the Enrichr DisGeNET library. The background universe was "
        f"{run.universe_source} ({run.universe_size} genes); for each database only annotated genes were used "
        f"as background. Sets with {rq.min_size}-{rq.max_size} genes were tested. "
        f"The following tools were run: {tool_versions}. "
        "Local over-representation analysis used the one-sided hypergeometric test with Benjamini-Hochberg "
        "correction per database; ranked lists were additionally analysed with a minimum-hypergeometric "
        "(ordered) test. External services were queried with all results requested where supported. "
        "Terms were matched across tools by identifier (GO, KEGG, Reactome, WikiPathways) or normalised name. "
        "For each term a consensus score (0-100) combined family-balanced median significance (40%), "
        "fraction of eligible tools reporting padj < "
        f"{rq.alpha} (30%), effect-size percentile (15%) and direction consistency (15%), multiplied by a "
        "specificity factor penalising very large or very general terms; tools only count as evidence for "
        "terms in databases they cover. Tiers: ROBUST (significant in >=2 tools, >=60% support, median "
        "padj <= 0.01, consistent direction, >=2 method families when available), SUPPORTED (>=2 tools, "
        ">=30% support), SINGLE-TOOL. Rank consistency was summarised by Robust Rank Aggregation (Kolde et "
        "al., Bioinformatics 2012). Significant terms were clustered into themes by Cohen's kappa of shared "
        "hit genes (average linkage, kappa 0.3), and key genes were ranked by their weighted membership in "
        "ROBUST/SUPPORTED themes.\n\nCite the underlying tools: Enrichr (Kuleshov et al. 2016), g:Profiler "
        "(Kolberg et al. 2023), Reactome (Milacic et al. 2024), STRING (Szklarczyk et al. 2023), "
        "clusterProfiler (Wu et al. 2021), ReactomePA (Yu & He 2016), fgsea (Korotkevich et al. 2021), "
        "blitzGSEA (Lachmann et al. 2022), Metascape (Zhou et al. 2019), MSigDB (Liberzon et al. 2015).\n")


def write_outputs(run, ctx, on_event=None):
    out = run.outdir
    for qr in run.queries:
        _emit(on_event, type="stage", stage="figures", query=qr.query.name)
        write_query(qr, run, out / _safe(qr.query.name))
    (out / "summary.txt").write_text(text_summary(run))
    (out / "methods.txt").write_text(methods_text(run))
    (out / "report.html").write_text(render_html(run))


def write_manifest(run):
    rq = run.request
    man = {
        "pathwayagent_version": __version__, "run_id": run.run_id, "created": datetime.now().isoformat(),
        "python": platform.python_version(), "organism": rq.organism,
        "parameters": {"databases": rq.databases, "tools": rq.tools, "min_size": rq.min_size,
                       "max_size": rq.max_size, "alpha": rq.alpha, "deadline_s": rq.deadline_s,
                       "local_only": rq.local_only, "min_ranked_for_gsea": rq.min_ranked_for_gsea},
        "universe": {"size": run.universe_size, "source": run.universe_source},
        "library": run.library_manifest, "timings": run.timings, "warnings": run.warnings,
        "queries": [{"name": qr.query.name, "n_genes_mapped": len(qr.query.genes), "ranked": qr.query.is_ranked,
                     "tools": [r.summary() for r in qr.tool_results]} for qr in run.queries],
        "families": FAMILIES, "database_labels": DATABASES,
    }
    (run.outdir / "manifest.json").write_text(json.dumps(man, indent=2, default=str))
