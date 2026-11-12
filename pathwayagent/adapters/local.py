"""Local engine adapters (design §7.2): the Enrichr-speed backbone."""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..core import stats
from ..core.schema import QueryInput
from .base import Adapter, RunContext


def _per_db(ctx: RunContext, dbs, fn) -> pd.DataFrame:
    """Run fn(rows) per database so BH is applied within each database (its own test family)."""
    lib = ctx.library
    frames = []
    for db in dbs:
        rows = lib.db_rows([db])
        df = fn(rows)
        if df is None or df.empty:
            continue
        df["database"] = db
        df["padj"] = stats.bh(df["pval"].values)
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _decorate(df: pd.DataFrame, ctx: RunContext, base: dict) -> pd.DataFrame:
    if df.empty:
        return df
    s = ctx.library.sets
    df["term_id"] = s["set_id"].values[df["row"].values]
    df["term_key"] = s["term_key"].values[df["row"].values]
    df["term_name"] = s["name"].values[df["row"].values]
    for k, v in base.items():
        df[k] = v
    return df.drop(columns=["row"])


class LocalORA(Adapter):
    name, family, method, mode = "local_ora", "F1_ORA", "hypergeometric", "shared_library"
    databases = ("GO_BP", "GO_MF", "GO_CC", "KEGG", "REACTOME", "WIKIPATHWAYS", "HALLMARK", "DISGENET")
    version = "pathwayagent"

    def _run(self, q: QueryInput, ctx: RunContext):
        lib, r = ctx.library, ctx.request
        qidx = lib.gene_indices(q.genes)
        df = _per_db(ctx, self.covered(ctx), lambda rows: stats.ora(
            lib.matrix, rows, ctx.universe_idx, qidx, lib.genes, r.min_size, r.max_size))
        df = _decorate(df, ctx, {"query": q.name, "tool": self.name, "method": self.method,
                                 "family": self.family, "mode": self.mode, "padj_method": "BH",
                                 "effect_type": "odds_ratio", "direction": 0})
        return df, {}


class LocalMHG(Adapter):
    """Ordered ORA over prefixes of the ranked list; works for partial rankings (e.g. a module)."""
    name, family, method, mode = "local_mhg", "F3_RANK", "min-hypergeometric (ordered)", "shared_library"
    databases = LocalORA.databases
    needs_scores = True
    version = "pathwayagent"

    def _run(self, q: QueryInput, ctx: RunContext):
        lib, r = ctx.library, ctx.request
        ranked = lib.gene_indices(q.ranked_genes())
        df = _per_db(ctx, self.covered(ctx), lambda rows: stats.mhg(
            lib.matrix, rows, ctx.universe_idx, ranked, lib.genes, r.min_size, r.max_size))
        df = _decorate(df, ctx, {"query": q.name, "tool": self.name, "method": self.method,
                                 "family": self.family, "mode": self.mode,
                                 "padj_method": "BH (mHG Bonferroni over cut-offs)",
                                 "effect_type": "fold_enrichment_at_cutoff", "direction": 0})
        return df, {}


class LocalRank(Adapter):
    name, family, method, mode = "local_rank", "F3_RANK", "Wilcoxon rank-sum", "shared_library"
    databases = LocalORA.databases
    needs_ranked = True
    version = "pathwayagent"

    def _run(self, q: QueryInput, ctx: RunContext):
        lib, r = ctx.library, ctx.request
        genes = [g for g in q.genes if g in lib.gene_pos and g in ctx.universe_symbols]
        idx = lib.gene_indices(genes)
        sc = np.array([q.scores[g] for g in genes], dtype=float)
        df = _per_db(ctx, self.covered(ctx), lambda rows: stats.rank_test(
            lib.matrix, rows, idx, sc, lib.genes, r.min_size, r.max_size))
        df = _decorate(df, ctx, {"query": q.name, "tool": self.name, "method": self.method,
                                 "family": self.family, "mode": self.mode, "padj_method": "BH",
                                 "effect_type": "AUC"})
        return df, {}


class LocalGSEA(Adapter):
    name, family, method, mode = "local_gsea", "F2_GSEA", "pre-ranked GSEA (blitzGSEA)", "shared_library"
    databases = ("KEGG", "REACTOME", "WIKIPATHWAYS", "HALLMARK", "GO_BP", "GO_MF", "GO_CC")
    needs_ranked = True
    version = "blitzgsea"

    def _run(self, q: QueryInput, ctx: RunContext):
        lib, r = ctx.library, ctx.request
        sig = pd.Series({g: s for g, s in q.scores.items() if g in ctx.universe_symbols})
        frames = []
        for db in self.covered(ctx):
            rows = lib.db_rows([db])
            gsets = {lib.sets["set_id"].iat[i]: lib.set_genes(i) for i in rows}
            res = stats.gsea_blitz(sig, gsets, r.min_size, r.max_size)
            if res.empty:
                continue
            res["database"] = db
            frames.append(res)
        if not frames:
            return pd.DataFrame(), {}
        res = pd.concat(frames, ignore_index=True)
        s = lib.sets.set_index("set_id")
        out = pd.DataFrame({
            "query": q.name, "tool": self.name, "method": self.method, "family": self.family,
            "mode": self.mode, "database": res["database"], "term_id": res["set_id"],
            "term_key": s.loc[res["set_id"], "term_key"].values,
            "term_name": s.loc[res["set_id"], "name"].values,
            "term_size": res.get("geneset_size"), "query_size": len(sig), "universe_size": len(sig),
            "overlap": res.get("geneset_size"), "pval": res["pval"], "padj": res["fdr"],
            "padj_method": "blitzGSEA FDR", "effect": res["nes"], "effect_type": "NES",
            "direction": np.sign(res["nes"]).astype(int),
            "hit_genes": res["leading_edge"].apply(lambda x: str(x).replace(",", ";")),
        })
        return out, {}
