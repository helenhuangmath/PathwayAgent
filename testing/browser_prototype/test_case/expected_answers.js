// Tiny test case — the right answers, worked out by hand and confirmed (SPEC.md table, rows 1–8).
// DO NOT edit these to make a test pass. If the app disagrees, the app is wrong.
// Numbers are compared to 4 significant figures.
var EXPECTED_ANSWERS = [
  { row: 1, group: "G1", check: "warnings", excel: ["1-Sep"], unrecognised: ["FAKE1"], duplicates: ["gzmb"], n: 4 },
  { row: 2, group: "G1", pathway: "EFFECTOR",   k: 3, K: 3, p: 0.03333, fdr: 0.06667 },
  { row: 3, group: "G1", pathway: "EXHAUSTION", k: 1, K: 4, p: 0.9286,  fdr: 0.9286 },
  { row: 4, group: "G1", pathway: "MEMORY",     notShown: true },
  { row: 5, group: "G2", check: "warnings", excel: [], unrecognised: [], duplicates: [], n: 5 },
  { row: 6, group: "G2", pathway: "MEMORY",     k: 3, K: 3, p: 0.08333, fdr: 0.1667 },
  { row: 7, group: "G2", pathway: "EXHAUSTION", k: 2, K: 4, p: 0.7381,  fdr: 0.7381 },
  { row: 8, group: "G2", pathway: "EFFECTOR",   notShown: true }
];
