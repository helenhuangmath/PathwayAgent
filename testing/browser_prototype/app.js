// Pathway Quick-Run — page behaviour (buttons, file loading, tables, downloads).
(function () {
  "use strict";
  const $ = id => document.getElementById(id);
  const DB_LABEL = { Hallmark: "Hallmark", GO_BP: "GO BP", Reactome: "Reactome", WikiPathways: "WikiPathways" };
  let lastResults = null, lastPlot = null;

  const species = () => document.querySelector('input[name="species"]:checked').value;
  const chosenDbs = () => Array.from(document.querySelectorAll("#dbs input:checked")).map(i => i.value);
  const allLibs = sp => PQ.listLibraries(sp).map(r => PQ.getLibrary(sp, r.db));

  // ---------- warnings (shown before running) ----------
  function renderWarnings() {
    const groups = PQ.parseGroupsText($("genes").value);
    $("groupInfo").textContent = groups.length
      ? groups.length + " group(s): " + groups.map(g => g.name + " (" + g.genes.length + ")").join(", ")
      : "";
    const box = $("warnings");
    box.innerHTML = "";
    if (!groups.length) return;
    const recognised = PQ.recognisedFor(allLibs(species()));
    let any = false;
    for (const g of groups) {
      const { clean, warnings: w } = PQ.checkGenes(g.genes, recognised);
      const parts = [];
      if (w.excel.length) parts.push("Excel-mangled (likely a date, not a gene): " + w.excel.join(", "));
      if (w.unrecognised.length) parts.push("Not recognised for " + (species() === "Mm" ? "mouse" : "human") + ": " + w.unrecognised.join(", "));
      if (w.duplicates.length) parts.push("Duplicates (counted once): " + w.duplicates.join(", "));
      if (parts.length) {
        any = true;
        const d = document.createElement("div");
        d.className = "warn";
        d.innerHTML = "<b>" + esc(g.name) + "</b> — " + clean.length + " of " + g.genes.length + " names used. " +
          parts.map(esc).join(" · ") + ". These are left out, not corrected.";
        box.appendChild(d);
      }
    }
    if (!any) {
      const d = document.createElement("div");
      d.className = "ok";
      d.textContent = "All gene names recognised; no duplicates or Excel-mangled names.";
      box.appendChild(d);
    }
  }

  // ---------- run ----------
  function run() {
    const groups = PQ.parseGroupsText($("genes").value);
    if (!groups.length) { $("status").textContent = "Paste or load at least one gene list first."; return; }
    const dbs = chosenDbs();
    if (!dbs.length) { $("status").textContent = "Tick at least one database."; return; }
    $("run").disabled = true;
    $("status").textContent = "Running…";
    setTimeout(() => {
      const t0 = performance.now();
      const sp = species();
      const recognised = PQ.recognisedFor(allLibs(sp));
      const libs = dbs.map(db => PQ.getLibrary(sp, db));
      lastResults = PQ.runAll(groups, libs, recognised);
      renderPlot();
      renderTables();
      const secs = (performance.now() - t0) / 1000;
      $("status").textContent = "Done in " + secs.toFixed(2) + " s — " + groups.length + " group(s) × " + dbs.length + " database(s).";
      $("results").hidden = false;
      $("run").disabled = false;
    }, 20);
  }

  function renderPlot() {
    lastPlot = PQPlot.buildDotPlot(lastResults);
    const box = $("plot");
    if (!lastPlot) { box.textContent = "No pathways had any hits, so there is nothing to plot."; return; }
    box.innerHTML = PQPlot.toSVG(lastPlot);
  }

  function fmtP(v) { return v >= 0.001 ? v.toPrecision(4) : v.toExponential(3); }
  function esc(s) { return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;"); }

  function renderTables() {
    const box = $("tables");
    box.innerHTML = "";
    for (const g of lastResults) {
      const div = document.createElement("div");
      div.className = "group";
      const nInfo = g.perDb.map(d => DB_LABEL[d.db] + " n = " + d.n + " (background N = " + d.N + ")").join(" · ");
      const nSig = g.rows.filter(r => r.fdr < 0.05).length;
      let html = "<h3>" + esc(g.name) + "</h3><p class='hint'>" + g.clean.length + " genes used of " + g.input +
        " entered. " + nInfo + ". " + g.rows.length + " pathways with ≥1 hit; " + nSig + " with FDR &lt; 0.05.</p>";
      if (!g.rows.length) { html += "<p class='hint'>No pathway contained any of these genes.</p>"; div.innerHTML = html; box.appendChild(div); continue; }
      html += "<div class='tablewrap'><table><thead><tr><th>Database</th><th>Pathway</th><th>Overlap (k/K)</th><th>Genes in overlap</th><th>P-value</th><th>FDR</th></tr></thead><tbody>";
      for (const r of g.rows) {
        html += "<tr" + (r.fdr < 0.05 ? " class='sig'" : "") + "><td>" + DB_LABEL[r.db] + "</td><td class='path'>" + esc(r.pathway).replace(/_/g, "_<wbr>") +
          "</td><td class='num'>" + r.k + "/" + r.K + "</td><td class='genes'>" + esc(r.genes.join(", ")) +
          "</td><td class='num'>" + fmtP(r.p) + "</td><td class='num'>" + fmtP(r.fdr) + "</td></tr>";
      }
      html += "</tbody></table></div>";
      div.innerHTML = html;
      box.appendChild(div);
    }
  }

  // ---------- downloads ----------
  function download(blob, name) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.appendChild(a);
    a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
  }
  function csvCell(v) { const s = String(v); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; }
  function downloadCsv() {
    if (!lastResults) return;
    const lines = [["Group", "Database", "Pathway", "Overlap_k", "Pathway_size_K", "Genes_used_n", "Background_N", "Genes_in_overlap", "P_value", "FDR"].join(",")];
    for (const g of lastResults) for (const r of g.rows) {
      lines.push([g.name, r.db, r.pathway, r.k, r.K, r.n, r.N, r.genes.join(";"), r.p, r.fdr].map(csvCell).join(","));
    }
    download(new Blob([lines.join("\r\n") + "\r\n"], { type: "text/csv" }), "pathway_results.csv");
  }
  function downloadSvg() {
    if (!lastPlot) return;
    download(new Blob([PQPlot.toSVG(lastPlot)], { type: "image/svg+xml" }), "pathway_dotplot.svg");
  }
  function downloadPng() {
    if (!lastPlot) return;
    PQPlot.toCanvas(lastPlot, 3).toBlob(b => download(b, "pathway_dotplot.png"), "image/png");  // 3x for print quality
  }

  // ---------- file loading ----------
  function loadFile(ev) {
    const f = ev.target.files[0];
    if (!f) return;
    const reader = new FileReader();
    reader.onload = () => {
      const text = String(reader.result);
      $("genes").value = /\.csv$/i.test(f.name) ? PQ.groupsToText(PQ.parseGroupsCSV(text)) : text;
      renderWarnings();
    };
    reader.readAsText(f);
    ev.target.value = "";
  }

  // ---------- footer ----------
  function renderFooter() {
    const hs = PQ.listLibraries("Hs")[0], mm = PQ.listLibraries("Mm")[0];
    const ver = [hs, mm].filter(Boolean).map(x => x.version).join(" / ");
    const date = (hs || mm || {}).downloaded || "";
    $("footer").textContent = "Gene sets: " + ver + " (Hallmark, GO:BP, Reactome, WikiPathways; CC BY 4.0), downloaded " + date +
      ". Test: one-sided hypergeometric, Benjamini–Hochberg FDR; background = all genes in each database.";
  }

  let timer = null;
  $("genes").addEventListener("input", () => { clearTimeout(timer); timer = setTimeout(renderWarnings, 250); });
  document.querySelectorAll('input[name="species"]').forEach(r => r.addEventListener("change", renderWarnings));
  $("file").addEventListener("change", loadFile);
  $("clear").addEventListener("click", () => { $("genes").value = ""; renderWarnings(); });
  $("run").addEventListener("click", run);
  $("dlCsv").addEventListener("click", downloadCsv);
  $("dlSvg").addEventListener("click", downloadSvg);
  $("dlPng").addEventListener("click", downloadPng);
  renderFooter();
})();
