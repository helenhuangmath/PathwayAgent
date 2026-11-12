"""
Command line (design §13.3).

  python -m pathwayagent db build [--organism mmusculus hsapiens] [--force]
  python -m pathwayagent db list
  python -m pathwayagent run FILE [FILE ...] [--preset standard] [--tools ...] [--db ...]
                                   [--organism auto] [--universe FILE] [--out DIR] [--local-only]
  python -m pathwayagent serve [--port 8765] [--no-browser]
  python -m pathwayagent tools
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .config import DEFAULT_DATABASES, DEFAULTS, ORGANISMS, PRESETS, library_path


def _progress(ev: dict):
    t = ev.get("type")
    if t == "tool_done":
        mark = {"ok": "✓", "empty": "–", "skipped": "·", "error": "✗", "timeout": "⏱"}.get(ev["status"], "?")
        msg = f"  {mark} {ev['query'][:22]:<22} {ev['tool']:<20} {ev['status']:<8} {ev['runtime_s']:>6.1f}s " \
              f"{ev['n_significant']:>5} sig"
        if ev["status"] in ("error", "timeout", "skipped") and ev.get("message"):
            msg += f"   ({ev['message'][:100]})"
        print(msg, flush=True)
    elif t == "stage":
        if ev["stage"] == "prepared":
            print(f"[PathwayAgent] organism={ev['organism']} universe={ev['universe']} genes, "
                  f"lists={ev['queries']}", flush=True)
        elif ev["stage"] == "figures":
            print(f"[PathwayAgent] writing tables & figures for {ev['query']} …", flush=True)
    elif t == "done":
        print(f"[PathwayAgent] done in {ev['timings']['total_s']}s -> {ev['outdir']}", flush=True)


def cmd_run(a):
    from .core.inputs import parse_file
    from .core.schema import AnalysisRequest
    from .knowledge.profile import load_profile
    from .pipeline import run_analysis

    queries = []
    for f in a.files:
        qs = parse_file(f)
        if a.lists:
            qs = [q for q in qs if q.name in a.lists]
        queries.extend(qs)
    if not queries:
        sys.exit("No gene lists found in the input.")
    universe = None
    if a.universe:
        universe = [x.strip() for x in Path(a.universe).read_text().replace(",", "\n").split() if x.strip()]
    tools = a.tools.split(",") if a.tools else PRESETS[a.preset]
    profile = load_profile(a.profile)
    metascape = {}
    for item in a.metascape or []:                  # LIST=PATH, or PATH alone for every list
        name, _, path = item.rpartition("=")
        metascape[name or "*"] = path
    req = AnalysisRequest(queries=queries, organism=a.organism, universe=universe,
                          databases=a.db.split(",") if a.db else list(DEFAULT_DATABASES), tools=tools,
                          min_size=a.min_size, max_size=a.max_size, alpha=a.alpha, deadline_s=a.deadline,
                          local_only=a.local_only, profile=profile, metascape=metascape or None,
                          title=a.title or "PathwayAgent run")
    run = run_analysis(req, outdir=a.out, on_event=_progress)
    print(f"\nReport: {run.outdir / 'report.html'}\nSummary: {run.outdir / 'summary.txt'}")
    if a.print_summary:
        print((run.outdir / "summary.txt").read_text())
    return run


def cmd_db(a):
    if a.action == "build":
        from .db.build_library import compile_library
        for org in a.organism:
            compile_library(org, force=a.force)
            try:
                from .db.landscape import build_all
                build_all(org)                      # term-landscape maps (needs scikit-learn)
            except ImportError as exc:
                print(f"  (skipping term landscapes: {exc})")
            from .adapters.rtools import RWorker, r_available
            if r_available():                       # one-time org.db GO table for clusterprofiler_go
                print(f"  caching org.db GO annotation for R ({org}; first time only, ~2 min) …")
                w = RWorker()
                try:
                    w.call("enrichgo", {"genes": ["x"], "universe": [], "organism": org, "onts": ["BP", "MF", "CC"],
                                        "cache_dir": str(library_path(org) / "rcache"), "min_size": 10,
                                        "max_size": 500}, timeout=900)
                except Exception as exc:
                    print(f"  (R GO cache not built: {exc})")
                finally:
                    w.stop()
    else:
        for org in ORGANISMS:
            p = library_path(org) / "manifest.json"
            if p.exists():
                m = json.loads(p.read_text())
                print(f"{org}: built {m['built']}, MSigDB {m.get('msigdb_release')}")
                for db, n in m["n_sets"].items():
                    print(f"   {db:<14} {n:>6} sets")
            else:
                print(f"{org}: not built (python -m pathwayagent db build --organism {org})")


def cmd_tools(_a):
    from .adapters import PLANNED, REGISTRY
    for name, cls in REGISTRY.items():
        print(f"  {name:<20} {cls.family:<10} {'web' if cls.web else 'local':<6} {cls.method}  "
              f"[{','.join(cls.databases)}]")
    print("\nPlanned / not yet available:")
    for k, v in PLANNED.items():
        print(f"  {k:<20} {v}")
    print("\nPresets:")
    for k, v in PRESETS.items():
        print(f"  {k:<10} {','.join(v)}")


def cmd_serve(a):
    from .server.app import serve
    serve(port=a.port, open_browser=not a.no_browser, host=a.host)


def main(argv=None):
    p = argparse.ArgumentParser(prog="pathwayagent", description="Fast multi-tool pathway analysis with consensus")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="analyse one or more gene-list files")
    r.add_argument("files", nargs="+")
    r.add_argument("--organism", default="auto", choices=["auto", *ORGANISMS])
    r.add_argument("--preset", default="standard", choices=list(PRESETS))
    r.add_argument("--tools", help="comma-separated tool names (overrides --preset); see `tools`")
    r.add_argument("--db", help=f"comma-separated databases (default: {','.join(DEFAULT_DATABASES)})")
    r.add_argument("--lists", nargs="*", help="only analyse these list names (multi-column files)")
    r.add_argument("--universe", help="background gene file (recommended)")
    r.add_argument("--profile", help="project profile YAML (genes/pathways of interest, background)")
    r.add_argument("--out", help="output directory (default results/<timestamp>)")
    r.add_argument("--min-size", type=int, default=DEFAULTS["min_size"])
    r.add_argument("--max-size", type=int, default=DEFAULTS["max_size"])
    r.add_argument("--alpha", type=float, default=DEFAULTS["alpha"])
    r.add_argument("--deadline", type=float, default=DEFAULTS["deadline_s"])
    r.add_argument("--local-only", action="store_true", help="never send genes to external services")
    r.add_argument("--metascape", nargs="*", metavar="[LIST=]PATH",
                   help="Metascape result package(s) (zip or folder) to include as the `metascape` tool")
    r.add_argument("--title")
    r.add_argument("--print-summary", action="store_true")
    r.set_defaults(fn=cmd_run)

    d = sub.add_parser("db", help="build or list the local gene-set library")
    d.add_argument("action", choices=["build", "list"])
    d.add_argument("--organism", nargs="+", default=["mmusculus", "hsapiens"], choices=list(ORGANISMS))
    d.add_argument("--force", action="store_true", help="re-download sources")
    d.set_defaults(fn=cmd_db)

    s = sub.add_parser("serve", help="start the local web app")
    s.add_argument("--port", type=int, default=8765)
    s.add_argument("--host", default="127.0.0.1")
    s.add_argument("--no-browser", action="store_true")
    s.set_defaults(fn=cmd_serve)

    t = sub.add_parser("tools", help="list available tools and presets")
    t.set_defaults(fn=cmd_tools)

    a = p.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
