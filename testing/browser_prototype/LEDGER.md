# LEDGER

## 2026-09-24 — Pathway Quick-Run v1 (build of SPEC.md)

### Asked
Build exactly what SPEC.md describes: an offline, double-click web page (`index.html`) that takes one or
more gene lists, runs Enrichr-style over-representation against bundled MSigDB 2026.1 Hallmark, GO BP,
Reactome and WikiPathways (human or mouse), warns about bad gene names, and gives a results table (CSV)
and a group-comparison dot plot (SVG/PNG). Save the tiny test case and its answers as files, with a
`test.html` to double-click.

### Agent did
- Downloaded (with permission) 8 MSigDB 2026.1 `.gmt` files (~12 MB) into `data_source/`, and converted
  them with `tools/make_data.py` into `data/*.js` so the page can read them offline.
- Wrote `engine.js` (input reading, gene-name checks, hypergeometric test, Benjamini–Hochberg FDR),
  `plot.js` (dot plot, drawn identically as SVG and PNG), `app.js` + `app.css` + `index.html` (the page).
- Saved the tiny test case as `test_case/toy_database.js`, `test_case/input.js`,
  `test_case/expected_answers.js`, and wrote `test.html`, which runs it through the same `engine.js` the app uses.
- Fixed two display bugs found by screenshot: long pathway names squashing the table, and
  overlapping "0" / "FDR 0.05" legend labels.

### Checked how
- **Tiny test:** `test.html` opened from disk in real Chrome (headless): **All 8 checks PASS**. Also re-ran
  the same files through Node: same numbers. The code matched the hand-worked answers first time;
  expected answers were not touched.
- **Independent maths check:** recomputed every result for the Tex/Tmem example (mouse, all 4 databases)
  in R 4.1.3 using `phyper` + `p.adjust("BH")`, straight from the original downloaded `.gmt` files:
  2,394/2,394 rows match, k/K/n/N identical, P and FDR differ by at most 2×10⁻¹⁰ (relative);
  773 pathways with FDR < 0.05 in both.
- **Speed (real clock, Node):** 3 lists × 4 mouse databases = 0.07 s, plus 0.2 s the first time to read the databases.
- **App in Chrome (headless, screenshots):** Tex/Tmem example loads, 40/40 genes recognised, dot plot and
  both tables render. With `1-Sep`, `FAKE1` and a duplicate `tox` added, all three are flagged and left out.
  SVG export parses as valid SVG; PNG export produces an image.
- **CSV input:** the SPEC's CSV example parses to groups Tex (3 genes) and Tmem (2 genes).
- **Offline:** no web addresses in the app code except the standard SVG namespace label.

### Confidently wrong
- The first version of the results table *looked* finished in code but was unreadable on screen
  (long GO names without spaces pushed columns off the edge). Only caught by taking a screenshot.
- The built-in preview pane showed `test.html` stuck on "Running…". It displays local files as static
  snapshots, so this was not an app bug. Real Chrome was used instead.
- Headless Chrome reports "Done in 0.00 s" because it fakes the clock; real timing was measured in Node.
- Biology note (not a bug): with all four databases ticked, the dot plot's top 5 are all GO BP terms,
  and the Tex top hits are broad "T cell activation" terms rather than an "exhaustion" pathway.
  None of the four databases has a dedicated T cell exhaustion gene set.

### Keep
- Hand-worked tiny test + an independent recomputation in R caught nothing wrong in the maths.
  Keep both whenever the engine changes (`test.html`, and the R check described above).
- Always screenshot the real page: the code-level checks all passed while the table was unreadable.
- To update databases: put new `.gmt` files in `data_source/`, edit the file names in
  `tools/make_data.py`, run `python3 tools/make_data.py`, then re-open `test.html`.
