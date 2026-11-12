// Pathway Quick-Run — calculation engine (see SPEC.md "Method").
// Shared by index.html (the app) and test.html (the tiny test), so both run the same code.
(function (root) {
  "use strict";

  // ---------- Gene-set libraries ----------
  const registered = [];   // raw data from data/*.js
  const parsedCache = {};

  function registerLibrary(info) { registered.push(info); }

  function listLibraries(species) {
    return registered.filter(r => r.species === species);
  }

  // GMT text: one pathway per line, "NAME<TAB>GENE<TAB>GENE..." (URL column already removed,
  // but a full MSigDB line with a URL in column 2 is also accepted).
  function parseGMT(text) {
    const terms = [];
    const universe = new Map();            // UPPERCASE symbol -> database spelling
    for (const line of text.split(/\r?\n/)) {
      if (!line.trim()) continue;
      let parts = line.split("\t").map(s => s.trim()).filter(Boolean);
      if (parts.length > 1 && /^https?:\/\//.test(parts[1])) parts.splice(1, 1);
      const name = parts[0];
      const genes = [];
      const seen = new Set();
      for (const g of parts.slice(1)) {
        const up = g.toUpperCase();
        if (seen.has(up)) continue;
        seen.add(up);
        genes.push(up);
        if (!universe.has(up)) universe.set(up, g);
      }
      terms.push({ name, genes });
    }
    return { terms, universe };
  }

  function getLibrary(species, db) {
    const key = species + "|" + db;
    if (!parsedCache[key]) {
      const info = registered.find(r => r.species === species && r.db === db);
      if (!info) throw new Error("Library not found: " + key);
      parsedCache[key] = Object.assign({ db, version: info.version, source: info.source }, parseGMT(info.gmt));
    }
    return parsedCache[key];
  }

  // "Recognised gene" = appears in any bundled database for that species (SPEC open decision 2).
  function recognisedFor(libs) {
    const map = new Map();
    for (const lib of libs) for (const [up, g] of lib.universe) if (!map.has(up)) map.set(up, g);
    return map;
  }

  // ---------- Reading the user's input ----------
  function cleanToken(t) { return t.trim().replace(/^["']+|["']+$/g, ""); }

  // Pasted text / .txt: genes separated by new lines, commas, tabs or spaces; ">Name" starts a group.
  function parseGroupsText(text) {
    const groups = [];
    let cur = null;
    for (const raw of text.split(/\r?\n/)) {
      const line = raw.trim();
      if (!line) continue;
      if (line.startsWith(">")) {
        cur = { name: line.slice(1).trim() || "List" + (groups.length + 1), genes: [] };
        groups.push(cur);
        continue;
      }
      if (!cur) { cur = { name: "List1", genes: [] }; groups.push(cur); }
      for (const t of line.split(/[\s,;]+/)) {
        const g = cleanToken(t);
        if (g) cur.genes.push(g);
      }
    }
    return groups.filter(g => g.genes.length > 0);
  }

  function splitCSVLine(line) {
    const out = [];
    let cell = "", inQ = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (inQ) {
        if (c === '"' && line[i + 1] === '"') { cell += '"'; i++; }
        else if (c === '"') inQ = false;
        else cell += c;
      } else if (c === '"') inQ = true;
      else if (c === ",") { out.push(cell); cell = ""; }
      else cell += c;
    }
    out.push(cell);
    return out;
  }

  // .csv (Metascape style): one column per group, first row = group names.
  function parseGroupsCSV(text) {
    const rows = text.split(/\r?\n/).filter(r => r.trim()).map(splitCSVLine);
    if (!rows.length) return [];
    const groups = rows[0].map((h, i) => ({ name: cleanToken(h) || "List" + (i + 1), genes: [] }));
    for (const row of rows.slice(1)) {
      row.forEach((cell, i) => {
        const g = cleanToken(cell);
        if (g && groups[i]) groups[i].genes.push(g);
      });
    }
    return groups.filter(g => g.genes.length > 0);
  }

  function groupsToText(groups) {
    return groups.map(g => ">" + g.name + "\n" + g.genes.join("\n")).join("\n");
  }

  // ---------- Checking gene names ----------
  const MONTHS = "jan|feb|mar|apr|may|jun|jul|aug|sep|sept|oct|nov|dec";
  const EXCEL_RE = new RegExp("^(\\d{1,2}-(" + MONTHS + ")|(" + MONTHS + ")-\\d{1,2})$", "i");

  function isExcelMangled(g) { return EXCEL_RE.test(g); }

  // Returns the de-duplicated, recognised genes (database spelling) and the warnings.
  // Never "fixes" a name: problem names are only reported and left out.
  function checkGenes(genes, recognised) {
    const seen = new Set();
    const clean = [];
    const warnings = { excel: [], unrecognised: [], duplicates: [] };
    for (const g of genes) {
      const up = g.toUpperCase();
      if (seen.has(up)) { warnings.duplicates.push(g); continue; }
      seen.add(up);
      if (isExcelMangled(g)) { warnings.excel.push(g); continue; }
      if (!recognised.has(up)) { warnings.unrecognised.push(g); continue; }
      clean.push(recognised.get(up));
    }
    return { clean, warnings };
  }

  // ---------- Statistics ----------
  const LF = [0]; // LF[i] = log(i!)
  function logFact(n) {
    for (let i = LF.length; i <= n; i++) LF[i] = LF[i - 1] + Math.log(i);
    return LF[n];
  }
  function logChoose(a, b) { return logFact(a) - logFact(b) - logFact(a - b); }

  // P(X >= k) for X ~ Hypergeometric(N background genes, K in pathway, n drawn).
  function hyperUpper(N, K, n, k) {
    const hi = Math.min(K, n);
    const denom = logChoose(N, n);
    let s = 0;
    for (let x = k; x <= hi; x++) {
      if (n - x > N - K) continue;
      s += Math.exp(logChoose(K, x) + logChoose(N - K, n - x) - denom);
    }
    return Math.min(1, s);
  }

  // Benjamini–Hochberg adjusted P-values (same order as input).
  function bh(ps) {
    const m = ps.length;
    const order = ps.map((p, i) => i).sort((a, b) => ps[a] - ps[b]);
    const q = new Array(m);
    let prev = 1;
    for (let r = m - 1; r >= 0; r--) {
      const i = order[r];
      prev = Math.min(prev, ps[i] * m / (r + 1));
      q[i] = prev;
    }
    return q;
  }

  // One group x one database. Only pathways with k >= 1 are tested and returned.
  function enrich(cleanGenes, lib) {
    const inN = cleanGenes.filter(g => lib.universe.has(g.toUpperCase()));
    const q = new Set(inN.map(g => g.toUpperCase()));
    const N = lib.universe.size, n = q.size;
    const rows = [];
    if (n === 0) return { N, n, rows };
    for (const t of lib.terms) {
      const hits = [];
      for (const g of t.genes) if (q.has(g)) hits.push(lib.universe.get(g));
      if (hits.length === 0) continue;
      rows.push({ db: lib.db, pathway: t.name, k: hits.length, K: t.genes.length, n, N,
                  genes: hits, p: hyperUpper(N, t.genes.length, n, hits.length) });
    }
    const fdr = bh(rows.map(r => r.p));
    rows.forEach((r, i) => { r.fdr = fdr[i]; });
    return { N, n, rows };
  }

  // Everything for all groups: warnings, n per database, and results sorted by P-value.
  function runAll(groups, libs, recognised) {
    return groups.map(g => {
      const { clean, warnings } = checkGenes(g.genes, recognised);
      const perDb = [];
      let rows = [];
      for (const lib of libs) {
        const r = enrich(clean, lib);
        perDb.push({ db: lib.db, n: r.n, N: r.N });
        rows = rows.concat(r.rows);
      }
      rows.sort((a, b) => a.p - b.p || a.fdr - b.fdr);
      return { name: g.name, input: g.genes.length, clean, warnings, perDb, rows };
    });
  }

  const PQ = { registerLibrary, listLibraries, getLibrary, parseGMT, recognisedFor,
               parseGroupsText, parseGroupsCSV, groupsToText, isExcelMangled, checkGenes,
               hyperUpper, bh, enrich, runAll };
  root.PQ = PQ;
  if (typeof module !== "undefined" && module.exports) module.exports = PQ;
})(typeof window !== "undefined" ? window : globalThis);
