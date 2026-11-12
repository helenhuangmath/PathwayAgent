// Pathway Quick-Run — dot plot (SPEC.md item 7).
// Builds one list of shapes, then draws it as SVG (for Illustrator) or onto a canvas (for PNG),
// so both downloads show exactly the same picture.
(function (root) {
  "use strict";
  const FONT = "Arial, Helvetica, sans-serif";
  const TOP_N = 5;

  // Sequential colour scale for -log10(FDR): pale yellow -> orange -> dark red.
  const STOPS = [[255, 237, 160], [253, 141, 60], [189, 0, 38]];
  function colour(t) {
    t = Math.max(0, Math.min(1, t));
    const seg = t < 0.5 ? 0 : 1, f = t < 0.5 ? t / 0.5 : (t - 0.5) / 0.5;
    const a = STOPS[seg], b = STOPS[seg + 1];
    const c = a.map((v, i) => Math.round(v + (b[i] - v) * f));
    return "rgb(" + c.join(",") + ")";
  }

  let measureCtx = null;
  function textWidth(s, size, bold) {
    if (!measureCtx) measureCtx = document.createElement("canvas").getContext("2d");
    measureCtx.font = (bold ? "bold " : "") + size + "px " + FONT;
    return measureCtx.measureText(s).width;
  }

  // results = PQ.runAll(...) output. Returns {width, height, shapes} or null if nothing to plot.
  function buildDotPlot(results) {
    // Rows: union of the top 5 pathways (by P-value) of each group, in the order groups first need them.
    const rowKeys = [];
    const seen = new Set();
    for (const g of results) {
      for (const r of g.rows.slice(0, TOP_N)) {
        const key = r.db + "|" + r.pathway;
        if (!seen.has(key)) { seen.add(key); rowKeys.push({ key, db: r.db, pathway: r.pathway }); }
      }
    }
    if (!rowKeys.length) return null;

    const cell = {};
    let maxK = 1, maxL = 0;
    for (const g of results) for (const r of g.rows) {
      const key = r.db + "|" + r.pathway;
      if (!seen.has(key)) continue;
      const l = -Math.log10(Math.max(r.fdr, 1e-300));
      cell[g.name + "\u0000" + key] = { k: r.k, l };
      maxK = Math.max(maxK, r.k); maxL = Math.max(maxL, l);
    }
    const colMax = Math.max(maxL, 1.30103);   // scale reaches at least FDR = 0.05

    const FS = 12, ROW_H = 24, R_MAX = 10;
    const labels = rowKeys.map(r => r.pathway + "  [" + r.db + "]");
    const labelW = Math.max(...labels.map(s => textWidth(s, FS))) + 12;
    const colW = Math.max(56, ...results.map(g => textWidth(g.name, FS, true) + 16));
    const top = 34, left = 12 + labelW;
    const plotW = colW * results.length, plotH = ROW_H * rowKeys.length;
    const legendX = left + plotW + 28;
    const width = legendX + 150;
    const height = Math.max(top + plotH + 16, top + 250);
    const radius = k => 3 + (R_MAX - 3) * Math.sqrt(k / maxK);

    const S = [];
    S.push({ t: "rect", x: 0, y: 0, w: width, h: height, fill: "#ffffff" });
    // grid + frame
    for (let i = 0; i < rowKeys.length; i++) {
      const y = top + i * ROW_H + ROW_H / 2;
      S.push({ t: "line", x1: left, y1: y, x2: left + plotW, y2: y, stroke: "#e6e6e6", sw: 1 });
      S.push({ t: "text", x: left - 8, y: y + 4, s: labels[i], size: FS, anchor: "end", fill: "#222" });
    }
    results.forEach((g, j) => {
      const x = left + j * colW + colW / 2;
      S.push({ t: "line", x1: x, y1: top, x2: x, y2: top + plotH, stroke: "#e6e6e6", sw: 1 });
      S.push({ t: "text", x, y: top - 12, s: g.name, size: FS, anchor: "middle", fill: "#222", bold: true });
    });
    S.push({ t: "rectStroke", x: left, y: top, w: plotW, h: plotH, stroke: "#999", sw: 1 });
    // dots
    rowKeys.forEach((rk, i) => {
      results.forEach((g, j) => {
        const c = cell[g.name + "\u0000" + rk.key];
        if (!c) return;
        S.push({ t: "circle", cx: left + j * colW + colW / 2, cy: top + i * ROW_H + ROW_H / 2,
                 r: radius(c.k), fill: colour(c.l / colMax), stroke: "#555", sw: 0.6 });
      });
    });
    // legend: size
    let y = top + 4;
    S.push({ t: "text", x: legendX, y, s: "Genes in overlap", size: FS, anchor: "start", fill: "#222", bold: true });
    const sizes = Array.from(new Set([1, Math.round(maxK / 2), maxK])).filter(v => v >= 1);
    y += 18;
    for (const k of sizes) {
      const r = radius(k);
      S.push({ t: "circle", cx: legendX + R_MAX, cy: y + R_MAX - 4, r, fill: "#bbbbbb", stroke: "#555", sw: 0.6 });
      S.push({ t: "text", x: legendX + 2 * R_MAX + 8, y: y + R_MAX, s: String(k), size: FS, anchor: "start", fill: "#222" });
      y += 2 * R_MAX + 6;
    }
    // legend: colour
    y += 12;
    S.push({ t: "text", x: legendX, y, s: "−log10(FDR)", size: FS, anchor: "start", fill: "#222", bold: true });
    y += 8;
    const barH = 100, steps = 25;
    for (let s = 0; s < steps; s++) {
      const t = 1 - s / (steps - 1);
      S.push({ t: "rect", x: legendX, y: y + s * barH / steps, w: 14, h: barH / steps + 0.5, fill: colour(t) });
    }
    S.push({ t: "rectStroke", x: legendX, y, w: 14, h: barH, stroke: "#999", sw: 0.6 });
    const fmt = v => (v >= 10 ? v.toFixed(0) : v.toFixed(1));
    S.push({ t: "text", x: legendX + 20, y: y + 9, s: fmt(colMax), size: FS - 1, anchor: "start", fill: "#222" });
    const y05 = y + barH * (1 - 1.30103 / colMax);
    if (y + barH - y05 > 14) S.push({ t: "text", x: legendX + 20, y: y + barH, s: "0", size: FS - 1, anchor: "start", fill: "#222" });
    S.push({ t: "line", x1: legendX - 3, y1: y05, x2: legendX + 17, y2: y05, stroke: "#222", sw: 1 });
    S.push({ t: "text", x: legendX + 20, y: y05 + 4, s: "FDR 0.05", size: FS - 1, anchor: "start", fill: "#222" });

    return { width: Math.ceil(width), height: Math.ceil(height), shapes: S };
  }

  function esc(s) { return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;"); }

  function toSVG(plot) {
    const out = ['<svg xmlns="http://www.w3.org/2000/svg" width="' + plot.width + '" height="' + plot.height +
                 '" viewBox="0 0 ' + plot.width + " " + plot.height + '" font-family="' + FONT + '">'];
    for (const s of plot.shapes) {
      if (s.t === "rect") out.push('<rect x="' + s.x + '" y="' + s.y + '" width="' + s.w + '" height="' + s.h + '" fill="' + s.fill + '"/>');
      else if (s.t === "rectStroke") out.push('<rect x="' + s.x + '" y="' + s.y + '" width="' + s.w + '" height="' + s.h + '" fill="none" stroke="' + s.stroke + '" stroke-width="' + s.sw + '"/>');
      else if (s.t === "line") out.push('<line x1="' + s.x1 + '" y1="' + s.y1 + '" x2="' + s.x2 + '" y2="' + s.y2 + '" stroke="' + s.stroke + '" stroke-width="' + s.sw + '"/>');
      else if (s.t === "circle") out.push('<circle cx="' + s.cx + '" cy="' + s.cy + '" r="' + s.r.toFixed(2) + '" fill="' + s.fill + '" stroke="' + s.stroke + '" stroke-width="' + s.sw + '"/>');
      else if (s.t === "text") out.push('<text x="' + s.x + '" y="' + s.y + '" font-size="' + s.size + '" text-anchor="' + s.anchor + '" fill="' + s.fill + '"' + (s.bold ? ' font-weight="bold"' : "") + ">" + esc(s.s) + "</text>");
    }
    out.push("</svg>");
    return out.join("\n");
  }

  function toCanvas(plot, scale) {
    const c = document.createElement("canvas");
    c.width = Math.ceil(plot.width * scale); c.height = Math.ceil(plot.height * scale);
    const ctx = c.getContext("2d");
    ctx.scale(scale, scale);
    for (const s of plot.shapes) {
      if (s.t === "rect") { ctx.fillStyle = s.fill; ctx.fillRect(s.x, s.y, s.w, s.h); }
      else if (s.t === "rectStroke") { ctx.strokeStyle = s.stroke; ctx.lineWidth = s.sw; ctx.strokeRect(s.x, s.y, s.w, s.h); }
      else if (s.t === "line") { ctx.strokeStyle = s.stroke; ctx.lineWidth = s.sw; ctx.beginPath(); ctx.moveTo(s.x1, s.y1); ctx.lineTo(s.x2, s.y2); ctx.stroke(); }
      else if (s.t === "circle") { ctx.beginPath(); ctx.arc(s.cx, s.cy, s.r, 0, 2 * Math.PI); ctx.fillStyle = s.fill; ctx.fill(); ctx.strokeStyle = s.stroke; ctx.lineWidth = s.sw; ctx.stroke(); }
      else if (s.t === "text") {
        ctx.font = (s.bold ? "bold " : "") + s.size + "px " + FONT;
        ctx.fillStyle = s.fill;
        ctx.textAlign = s.anchor === "end" ? "right" : s.anchor === "middle" ? "center" : "left";
        ctx.textBaseline = "alphabetic";
        ctx.fillText(s.s, s.x, s.y);
      }
    }
    return c;
  }

  root.PQPlot = { buildDotPlot, toSVG, toCanvas, TOP_N };
})(window);
