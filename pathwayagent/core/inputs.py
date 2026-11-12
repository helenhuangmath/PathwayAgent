"""
Parse user input into QueryInput objects (design §3.1).

Recognised shapes (auto-detected):
  * one column of genes                         -> gene list (ORA)
  * several columns with a header row           -> one gene list per column (e.g. modules)
  * gene + numeric column                       -> ranked list (gene, score)
  * DE table (gene, log2FoldChange, pvalue...)  -> ranked list by signed -log10(p) or stat;
                                                   plus up / down lists
"""
from __future__ import annotations

import io
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .schema import QueryInput

HEADER_WORDS = {"gene", "genes", "symbol", "gene_symbol", "genesymbol", "gene_name", "id",
                "gene_id", "name", "feature", "features"}
FC_COLS = ["log2foldchange", "log2fc", "logfc", "avg_log2fc", "avg_logfc"]
P_COLS = ["pvalue", "p_val", "pval", "p.value", "p"]
STAT_COLS = ["stat", "wald", "t", "maxwaldstat", "score", "statistic"]


def _clean(tokens) -> list[str]:
    out, seen = [], set()
    for t in tokens:
        t = str(t).strip().strip('"').strip("'").strip()
        if not t or t.lower() in ("nan", "na", "none"):
            continue
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out


def _read_table(text: str) -> pd.DataFrame:
    text = text.replace("\r\n", "\n").replace("\r", "\n").strip("\n")
    first = text.split("\n", 1)[0]
    sep = "\t" if "\t" in first else ("," if "," in first else None)
    if sep is None:
        # whitespace / single column
        return pd.DataFrame({0: re.split(r"[\s;]+", text.strip())})
    return pd.read_csv(io.StringIO(text), sep=sep, header=None, dtype=str, keep_default_na=False)


def _is_numeric(col: pd.Series, frac: float = 0.9) -> bool:
    vals = pd.to_numeric(col.replace("", np.nan), errors="coerce")
    nonempty = (col.astype(str).str.strip() != "").sum()
    return nonempty > 0 and vals.notna().sum() >= frac * nonempty


def parse_text(text: str, name: str = "query") -> list[QueryInput]:
    df = _read_table(text)
    if df.shape[1] == 1:
        genes = _clean(df.iloc[:, 0])
        if genes and genes[0].lower() in HEADER_WORDS:
            genes = genes[1:]
        return [QueryInput(name=name, genes=genes, raw_genes=list(genes))]

    header = [str(x).strip() for x in df.iloc[0].tolist()]
    body = df.iloc[1:].reset_index(drop=True)
    lower = [h.lower() for h in header]
    header_like = not _is_numeric(df.iloc[:1, 1]) if df.shape[1] > 1 else False

    # DE table / ranked list: first column genes, a numeric column present
    numeric_cols = [i for i in range(1, df.shape[1]) if _is_numeric(body.iloc[:, i])]
    if numeric_cols:
        if not header_like:  # no header row
            body, lower = df, [f"c{i}" for i in range(df.shape[1])]
        genes = body.iloc[:, 0].astype(str).str.strip()
        fc_i = next((lower.index(c) for c in FC_COLS if c in lower), None)
        p_i = next((lower.index(c) for c in P_COLS if c in lower), None)
        st_i = next((lower.index(c) for c in STAT_COLS if c in lower), None)
        if st_i is not None:
            score = pd.to_numeric(body.iloc[:, st_i], errors="coerce")
        elif fc_i is not None and p_i is not None:
            fc = pd.to_numeric(body.iloc[:, fc_i], errors="coerce")
            p = pd.to_numeric(body.iloc[:, p_i], errors="coerce").clip(lower=1e-300)
            score = np.sign(fc) * -np.log10(p)
        else:
            score = pd.to_numeric(body.iloc[:, numeric_cols[0]], errors="coerce")
        s = pd.Series(score.values, index=genes.values).dropna()
        s = s[~s.index.duplicated(keep="first")]
        s = s[[g != "" for g in s.index]]
        ordered = s.sort_values(ascending=False)
        q = QueryInput(name=name, genes=list(ordered.index), scores=ordered.to_dict(),
                       raw_genes=list(ordered.index))
        return [q]

    # several gene columns -> one list per column
    out = []
    has_header = any(h.lower() in HEADER_WORDS or re.search(r"module|cluster|list|set|up|down", h, re.I)
                     for h in header)
    if not has_header and len(body):
        # names built on a shared stem (TEX_G1, TEX_G2, ...) that the first data row does not share
        stem = os.path.commonprefix([h for h in header if h])
        row1 = [str(x).strip() for x in body.iloc[0].tolist()]
        has_header = len(stem) >= 3 and not any(x.startswith(stem) for x in row1 if x)
    cols = body if has_header else df
    for i in range(df.shape[1]):
        nm = header[i] if has_header and header[i] else f"{name}_{i + 1}"
        genes = _clean(cols.iloc[:, i])
        if genes:
            out.append(QueryInput(name=nm, genes=genes, raw_genes=list(genes)))
    return out


def parse_file(path: str | Path) -> list[QueryInput]:
    path = Path(path)
    if path.suffix.lower() in (".xlsx", ".xls"):
        text = pd.read_excel(path, header=None, dtype=str).to_csv(sep="\t", index=False, header=False)
    else:
        text = path.read_text(encoding="utf-8-sig", errors="replace")
    return parse_text(text, name=re.sub(r"[^\w.-]+", "_", path.stem))
