"""Self-contained HTML report (design §10.2): one tab per gene list, sortable tables,
interactive theme network (Cytoscape.js) with hover details, per-tool figures."""
from __future__ import annotations

import html
import json
from pathlib import Path

import pandas as pd

CYTOSCAPE = "https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.33.1/cytoscape.min.js"

CSS = """
/* Theme: flat Bootstrap-3 "Flatly" palette — dark slate navbar, teal links, white panels (matches the web app) */
@import url("https://fonts.googleapis.com/css?family=Lato:400,700,400italic&display=swap");
:root{--navy:#2c3e50;--navy-dark:#1a242f;--teal:#18bc9c;--teal-dark:#128f76;--blue:#3498db;--orange:#f39c12;
--red:#e74c3c;--grey:#95a5a6;--grey-2:#b4bcc2;--line:#ecf0f1;--line-2:#dddddd;--bg:#ffffff;--ink:#2c3e50;
--t1:#2c3e50;--t2:#18bc9c;--t3:#95a5a6;--ok:#18bc9c;--warn:#f39c12;--bad:#e74c3c;--surface:#ffffff;--ink-2:#7b8a8b}
*{box-sizing:border-box}html{scroll-padding-top:120px}
body{margin:0;padding-top:60px;background:var(--bg);color:var(--ink);font:15px/1.42857 "Lato","Helvetica Neue",Helvetica,Arial,sans-serif}
a{color:var(--teal);text-decoration:none}a:hover{color:var(--teal-dark);text-decoration:underline}
code{padding:2px 4px;font-size:90%;color:#c7254e;background:#f9f2f4;border-radius:4px}
.navbar{position:fixed;top:0;left:0;right:0;height:60px;background:var(--navy);z-index:1030}
.navbar .container{display:flex;align-items:center;height:60px}
.navbar-brand{color:#fff;font-size:19px;padding:0 15px;line-height:60px;white-space:nowrap}
.navbar-brand:hover{color:var(--teal);text-decoration:none}
.navbar-nav{display:flex;margin:0 0 0 auto;padding:0;list-style:none}
.navbar-nav a{display:block;color:#fff;padding:0 15px;line-height:60px;font-size:14px}
.navbar-nav a:hover{background:var(--navy-dark);color:var(--teal);text-decoration:none}
.container{width:100%;max-width:1400px;margin:0 auto;padding:0 15px}
header{padding:24px 0 6px}h1{margin:0 0 6px;font-size:32px;font-weight:400}
.meta{color:var(--ink-2);font-size:13px}
main{padding:0 0 30px}
.tabs{display:flex;gap:2px;flex-wrap:wrap;padding:8px 0;border-bottom:1px solid var(--line);position:sticky;top:60px;
background:var(--bg);z-index:5}
.tabs button{border:0;background:none;color:var(--teal);padding:10px 15px;border-radius:4px;cursor:pointer;font:inherit}
.tabs button:hover{background:var(--line)}.tabs button.active{background:var(--navy);color:#fff}
.q{display:none}.q.active{display:block}
h2{font-size:19px;font-weight:400;background:var(--line);color:var(--navy);padding:10px 15px;border-radius:4px;margin:28px 0 10px}
h3{font-size:17px;font-weight:400;margin:20px 0 6px;color:var(--navy)}
.cards{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}.card{background:#fff;border:1px solid var(--line);border-radius:4px;
padding:10px 16px;min-width:150px;box-shadow:0 1px 1px rgba(0,0,0,.05)}.card b{display:block;font-size:26px;font-weight:400;color:var(--navy)}
.card span{color:var(--ink-2);font-size:12px}
table{border-collapse:collapse;width:100%;font-size:13px;margin:0}
th,td{padding:8px;line-height:1.42857;border-top:1px solid var(--line);text-align:left;vertical-align:top}
thead th{cursor:pointer;color:var(--navy);border-top:0;border-bottom:2px solid var(--line);position:sticky;top:0;background:#fff}
tbody tr:nth-child(odd){background:#f9f9f9}tbody tr:hover{background:var(--line)}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.scroll{max-height:560px;overflow:auto;border:1px solid var(--line);border-radius:4px}
.badge{display:inline-block;padding:2px 7px;border-radius:10px;font-size:11px;font-weight:700;color:#fff;white-space:nowrap}
.b1{background:var(--t1)}.b2{background:var(--t2)}.b3{background:var(--t3)}
.st-ok{color:var(--ok);font-weight:700}.st-error,.st-timeout{color:var(--bad);font-weight:700}.st-skipped,.st-empty{color:var(--grey-2)}
.warn{padding:10px 15px;border-radius:4px;background:var(--orange);color:#fff;margin:8px 0;font-size:14px}
img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:4px;background:#fff}
figure{margin:0}figcaption{font-size:12px;text-align:right}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(460px,1fr));gap:15px}
details{border:1px solid var(--line);border-radius:4px;padding:8px 15px;margin:10px 0;background:#fff}
summary{cursor:pointer;color:var(--teal)}.genes{color:var(--ink-2);font-size:12px}
.cy{height:640px;border:1px solid var(--line);border-radius:4px;position:relative}
.tip{position:fixed;pointer-events:none;background:#fff;border:1px solid var(--line-2);border-radius:4px;
padding:8px 10px;font-size:12px;max-width:360px;box-shadow:0 4px 14px rgba(0,0,0,.15);display:none;z-index:20}
.footer-border{margin-top:30px;border-top:1px dashed var(--line-2)}
.footer-copyright{margin:15px 0 30px;line-height:24px;text-align:right;font-size:14px}
@media (max-width:600px){.grid2{grid-template-columns:1fr}.navbar-nav a{padding:0 8px;font-size:12px}}
"""

JS = """
function showTab(i){document.querySelectorAll('.q').forEach((e,j)=>e.classList.toggle('active',i===j));
document.querySelectorAll('.tabs button').forEach((e,j)=>e.classList.toggle('active',i===j));
}
document.addEventListener("toggle",e=>{const d=e.target;if(d.tagName==="DETAILS"&&d.open){const c=d.querySelector(".cy");if(c)drawNet(+c.id.slice(2));}},true);
document.addEventListener("click",e=>{const th=e.target.closest('th');if(!th)return;const tb=th.closest('table');
const i=[...th.parentNode.children].indexOf(th);const rows=[...tb.tBodies[0].rows];const asc=th.dataset.asc!=='1';
th.dataset.asc=asc?'1':'0';rows.sort((a,b)=>{let x=a.cells[i].dataset.v??a.cells[i].innerText,
y=b.cells[i].dataset.v??b.cells[i].innerText;const nx=parseFloat(x),ny=parseFloat(y);
if(!isNaN(nx)&&!isNaN(ny)){x=nx;y=ny}return (x>y?1:x<y?-1:0)*(asc?1:-1)});rows.forEach(r=>tb.tBodies[0].appendChild(r));});
const drawn={};
function drawNet(i){const el=document.getElementById('cy'+i);if(!el||drawn[i]||typeof cytoscape==='undefined')return;
drawn[i]=true;const data=NETS[i];if(!data)return;const cs=getComputedStyle(document.documentElement);
const col={1:cs.getPropertyValue('--t1'),2:cs.getPropertyValue('--t2'),3:cs.getPropertyValue('--t3'),0:'#b8b7b2'};
const cy=cytoscape({container:el,elements:data,layout:{name:'preset'},wheelSensitivity:0.2,style:[
{selector:'node',style:{'background-color':n=>col[n.data('tier')]||'#b8b7b2','width':n=>8+n.data('score')/4,
'height':n=>8+n.data('score')/4,'border-width':2,'border-color':cs.getPropertyValue('--surface')}},
{selector:'edge',style:{'width':e=>1+3*e.data('weight'),'line-color':cs.getPropertyValue('--line')}},
{selector:'node:selected',style:{'label':'data(label)','font-size':10,'color':cs.getPropertyValue('--ink')}}]});
const tip=document.getElementById('tip');
cy.on('mouseover','node',ev=>{const d=ev.target.data();tip.innerHTML='<b>'+d.label+'</b><br>'+d.db+
' · tier '+d.tier+' · score '+d.score.toFixed(1)+' · tools '+d.support+'<br><span class=genes>'+
(d.genes||'').split(';').slice(0,25).join(', ')+'</span>';tip.style.display='block'});
cy.on('mousemove',ev=>{tip.style.left=(ev.originalEvent.clientX+14)+'px';tip.style.top=(ev.originalEvent.clientY+14)+'px'});
cy.on('mouseout','node',()=>tip.style.display='none');}
window.addEventListener('load',()=>showTab(0));
"""


def _e(x) -> str:
    return html.escape("" if x is None or (isinstance(x, float) and pd.isna(x)) else str(x))


def _fmt(v, col=""):
    if isinstance(v, float):
        if pd.isna(v):
            return '<td class="num"></td>'
        if "padj" in col or col.endswith("_p") or col.endswith("_q") or col in ("pval",):
            return f'<td class="num" data-v="{v}">{v:.2g}</td>'
        return f'<td class="num" data-v="{v}">{v:.2f}</td>'
    if isinstance(v, (int,)) and not isinstance(v, bool):
        return f'<td class="num" data-v="{v}">{v}</td>'
    return f"<td>{_e(v)}</td>"


def table(df: pd.DataFrame, cols: list[str], n: int = 100, badge_col: str | None = None,
          gene_cols=("hit_genes", "core_genes")) -> str:
    if df is None or df.empty:
        return "<p class=meta>No rows.</p>"
    cols = [c for c in cols if c in df.columns]
    out = ['<div class="scroll"><table><thead><tr>' + "".join(f"<th>{_e(c)}</th>" for c in cols) + "</tr></thead><tbody>"]
    for _, r in df.head(n).iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if c == badge_col:
                t = {"ROBUST": 1, "SUPPORTED": 2, "SINGLE-TOOL": 3}.get(str(v), 0)
                cells.append(f'<td data-v="{t}"><span class="badge b{t}">{_e(v)}</span></td>')
            elif c in gene_cols and isinstance(v, str):
                g = v.split(";")
                cells.append(f'<td class="genes">{_e(", ".join(g[:15]))}{" …" if len(g) > 15 else ""}</td>')
            elif c == "term_name" and "url" in df.columns and isinstance(r.get("url"), str) and r.get("url"):
                cells.append(f'<td><a href="{_e(r["url"])}" target="_blank" rel="noopener">{_e(v)}</a></td>')
            else:
                cells.append(_fmt(v, c))
        out.append("<tr>" + "".join(cells) + "</tr>")
    out.append("</tbody></table></div>")
    if len(df) > n:
        out.append(f"<p class=meta>Showing {n} of {len(df)} rows — full table in the TSV/Excel files.</p>")
    return "".join(out)


def _img(path, base: Path) -> str:
    """PNG preview + link to the publication PDF written next to it."""
    if not path:
        return ""
    rel = Path(path).relative_to(base)
    pdf = rel.with_suffix(".pdf").as_posix()
    return (f'<figure><img loading="lazy" src="{_e(rel.as_posix())}" alt="{_e(rel.stem)}">'
            f'<figcaption><a href="{_e(pdf)}" target="_blank">PDF</a></figcaption></figure>')


def render_html(run) -> str:
    base = run.outdir
    rq = run.request
    nets = []
    tabs = []
    bodies = []
    for i, qr in enumerate(run.queries):
        q = qr.query
        c = qr.consensus
        tabs.append(f'<button onclick="showTab({i})">{_e(q.name)}</button>')
        nets.append(qr.figures.get("_cy") or None)
        n_in = len(q.mapping) if q.mapping is not None else len(q.genes)
        cnt = c["tier_label"].value_counts() if len(c) else {}
        b = [f'<section class="q" id="q{i}"><h2>{_e(q.name)}</h2>']
        b.append('<div class="cards">' + "".join(
            f'<div class="card"><b>{v}</b><span>{k}</span></div>' for k, v in [
                ("input genes", n_in), ("mapped", len(q.genes)),
                ("ROBUST terms", cnt.get("ROBUST", 0) if len(c) else 0),
                ("SUPPORTED terms", cnt.get("SUPPORTED", 0) if len(c) else 0),
                ("themes", len(qr.themes)),
                ("tools with results", sum(r.status == "ok" for r in qr.tool_results))]) + "</div>")
        for w in qr.warnings:
            b.append(f'<div class="warn">{_e(w)}</div>')
        # tool status
        st = pd.DataFrame([r.summary() for r in qr.tool_results])
        rows = "".join(
            f'<tr><td>{_e(s.tool)}</td><td>{_e(s.family)}</td><td class="st-{_e(s.status)}">{_e(s.status)}</td>'
            f'<td class="num">{s.runtime_s:.1f}s</td><td class="num">{s.n_significant}</td>'
            f'<td class="num">{s.n_tested}</td><td class="genes">{_e(s.message)}</td></tr>'
            for s in st.itertuples())
        b.append("<h3>Tools</h3><div class=scroll><table><thead><tr><th>tool</th><th>family</th><th>status</th>"
                 "<th>time</th><th>significant</th><th>tested</th><th>note</th></tr></thead><tbody>"
                 + rows + "</tbody></table></div>")
        # consensus
        b.append("<h2>Consensus pathways</h2><p class=meta>Sorted by tier, then consensus score. "
                 "Click a column header to sort. Tool support = tools with padj &lt; "
                 f"{rq.alpha} / tools that covered the term's database.</p>")
        b.append('<div class="grid2">' + _img(qr.figures.get("consensus_bar"), base)
                 + _img(qr.figures.get("evidence_bubble"), base) + "</div>")
        if len(c):
            b.append(table(c[c["tier"] > 0], ["rank", "tier_label", "term_name", "database", "consensus_score",
                                              "tool_support", "best_padj", "median_padj", "rra_q",
                                              "loo_rank_range", "depends_on", "n_hit_genes", "hit_genes"],
                           n=150, badge_col="tier_label"))
        # themes
        b.append("<h2>Themes</h2><p class=meta>Significant terms clustered by shared hit genes across "
                 "databases (kappa ≥ 0.3). Hover a node for details; click to show its label.</p>")
        b.append(_img(qr.figures.get("theme_network"), base))
        if nets[-1]:
            b.append(f'<details><summary>Interactive network (hover for genes)</summary>'
                     f'<div class="cy" id="cy{i}"></div></details>')
        b.append('<div class="grid2">' + _img(qr.figures.get("theme_tool"), base) + "".join(
            _img(v, base) for k, v in qr.figures.items() if k.startswith("landscape_")) + "</div>")
        b.append(table(qr.themes, ["theme_rank", "tier_label", "theme_name", "theme_score", "n_terms",
                                   "databases", "best_padj", "core_genes", "member_terms"], n=40,
                       badge_col="tier_label"))
        # key genes
        b.append("<h2>Key genes for follow-up</h2>")
        b.append('<div class="grid2">' + _img(qr.figures.get("key_gene_bubble"), base)
                 + _img(qr.figures.get("gene_theme"), base) + "</div>")
        b.append(table(qr.key_genes, ["rank", "gene", "gene_score", "n_themes", "n_tools", "input_score",
                                      "driver_of", "why"], n=60))
        b.append("<h2>Per-tool results</h2>" + _img(qr.figures.get("tool_overlap"), base)
                 + _img(qr.figures.get("metascape_clusters"), base))
        for r in qr.tool_results:
            if r.status != "ok":
                continue
            figs = qr.figures.get(f"tool:{r.tool}", [])
            top = r.records.sort_values("padj").head(20)
            b.append(f"<details><summary>{_e(r.tool)} — {_e(r.method)} ({r.summary()['n_significant']} "
                     f"significant)</summary>" + '<div class="grid2">' + "".join(_img(f, base) for f in figs)
                     + "</div>" + table(top, ["database", "term_name", "term_id", "overlap", "term_size", "pval",
                                               "padj", "effect", "hit_genes"], n=20) + "</details>")
        if q.mapping is not None:
            un = q.mapping[q.mapping["symbol"].isna()]
            b.append(f"<details><summary>ID mapping ({len(un)} unmapped)</summary>"
                     + table(q.mapping[q.mapping["status"] != "exact"], ["input", "symbol", "status"], n=300)
                     + "</details>")
        b.append("</section>")
        bodies.append("".join(b))

    warn = "".join(f'<div class="warn">{_e(w)}</div>' for w in run.warnings)
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>PathwayAgent report</title>
<style>{CSS}</style><script src="{CYTOSCAPE}"></script></head><body>
<nav class="navbar"><div class="container"><a class="navbar-brand" href="#">PathwayAgent</a>
<ul class="navbar-nav"><li><a href="summary.txt" target="_blank">Summary</a></li><li><a href="methods.txt" target="_blank">Methods</a></li>
<li><a href="manifest.json" target="_blank">Manifest</a></li><li><a href="#about">About</a></li></ul></div></nav>
<div class="container">
<header><h1>PathwayAgent report</h1><div class="meta">Run {_e(run.run_id)} · organism {_e(rq.organism)} ·
universe {run.universe_size} genes ({_e(run.universe_source)}) · databases {_e(', '.join(rq.databases))} ·
total {run.timings.get('total_s', run.timings.get('tools_s', ''))}s · files: summary.txt, methods.txt, per-list folders
</div>{warn}</header>
<nav class="tabs">{''.join(tabs)}</nav><main>{''.join(bodies)}</main>
<details id="about"><summary>About PathwayAgent</summary><p>PathwayAgent runs several pathway-enrichment tools at once and combines
them into a cross-tool consensus: robust themes, key genes and publication-ready figures. Developed by Hua Huang, Berger and Wherry labs,
University of Pennsylvania. MIT License. Cite: Huang H. PathwayAgent: fast multi-tool pathway analysis with cross-tool consensus (2026).
Contact: <a href="mailto:helenhuang.math@gmail.com">helenhuang.math@gmail.com</a>. The methods paragraph with tool and database versions is in
<a href="methods.txt" target="_blank">methods.txt</a>.</p></details>
<div class="footer-border"></div><div class="footer-copyright">Copyright © 2026 Hua Huang. PathwayAgent is released under the MIT License.</div>
</div><div class="tip" id="tip"></div>
<script>const NETS={json.dumps(nets)};{JS}</script></body></html>"""
