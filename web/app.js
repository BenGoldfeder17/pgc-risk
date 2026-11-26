// Page logic for the pgc-risk browser tool.  The files a visitor chooses are handed to a
// background thread of this same page (worker.js); nothing is sent over the network.
const $ = (id) => document.getElementById(id);
const el = (tag, cls, text) => {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (text !== undefined) n.textContent = text;
  return n;
};
let last = null;
let busy = 0;

const worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });

// ---------------------------------------------------------------- formatting
const pct = (x) => {
  const v = 100 * x;
  if (v < 0.095) return v.toFixed(2) + "%";
  if (v < 9.95) return v.toFixed(1) + "%";
  return Math.floor(v + 0.5) + "%";
};
const ordinal = (p) => {
  const v = Math.min(99, Math.max(1, Math.round(100 * p)));
  const suf = (v % 100 >= 10 && v % 100 <= 20) ? "th" : ({ 1: "st", 2: "nd", 3: "rd" }[v % 10] || "th");
  return v + suf;
};
const num = (n) => n.toLocaleString("en-US");
const FORMATS = { "23andme": "23andMe-style export", ancestry: "AncestryDNA export", csv: "MyHeritage / FamilyTreeDNA export" };
const BASIS = { published: "published", scaled: "published, scaled", "hold-out": "hold-out test", expected: "estimated from study size" };

// ---------------------------------------------------------------- choosing files
const drop = $("pr-drop"), input = $("pr-file");
drop.addEventListener("click", () => input.click());
drop.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); input.click(); } });
input.addEventListener("change", () => { take(input.files); input.value = ""; });
["dragenter", "dragover"].forEach((t) => document.addEventListener(t, (e) => { e.preventDefault(); drop.classList.add("over"); }));
["dragleave", "drop"].forEach((t) => document.addEventListener(t, (e) => { e.preventDefault(); drop.classList.remove("over"); }));
document.addEventListener("drop", (e) => { if (e.dataTransfer) take(e.dataTransfer.files); });

function take(files) {
  if (!files || !files.length) return;
  $("pr-error").hidden = true;
  for (const file of files) {
    busy++;
    worker.postMessage({ cmd: "file", file });
  }
  progress("Reading " + (files.length === 1 ? files[0].name : files.length + " files"));
}
function progress(text) {
  $("pr-progress").hidden = false;
  $("pr-stage").textContent = text;
}
function settle() {
  busy = Math.max(0, busy - 1);
  if (!busy) $("pr-progress").hidden = true;
}
function setData(info, kept) {
  const status = $("pr-data-status");
  status.replaceChildren();
  const badge = el("span", "pr-status " + (info ? "good" : ""));
  badge.append(el("i", "", info ? "✓" : "–"), document.createTextNode(info ? "Loaded" : "Not loaded yet"));
  status.append(badge);
  if (info) {
    status.append(document.createTextNode(`  ·  ${info.description}; ${num(info.snps)} SNPs, ${info.disorders} disorders, built ${info.created}` +
      (kept ? ". Kept in this browser." : ". Not kept: this browser would not store it, so load it again next time.")));
  }
  $("pr-forget").hidden = !info;
  $("pr-drop-title").textContent = info ? "Drop your genome file here, or click to choose it"
                                        : "Drop the reference data file and your genome file here, or click to choose them";
}

worker.onmessage = (e) => {
  const m = e.data;
  document.documentElement.dataset.prReady = "1";       // the scorer thread is up and has answered
  if (m.type === "bundle") { setData(m.info, m.kept); settle(); }
  else if (m.type === "no-bundle") { setData(null); }
  else if (m.type === "need-bundle") {
    settle();
    fail("Load the reference data first (step 1). " + m.name + " will be scored as soon as it is loaded.", "One more file needed");
  }
  else if (m.type === "progress") { progress(m.text); }
  else if (m.type === "result") { settle(); render(m.result); }
  else if (m.type === "error") { settle(); fail(m.message); }
};
worker.onerror = () => fail("This browser could not start the scorer. A current version of Chrome, Edge, Firefox or Safari is needed.");
worker.postMessage({ cmd: "restore" });
$("pr-forget").addEventListener("click", () => worker.postMessage({ cmd: "forget" }));

function fail(message, title) {
  $("pr-error-title").textContent = title || "Could not score this file";
  $("pr-error-text").textContent = message;
  $("pr-error").hidden = false;
}

// ---------------------------------------------------------------- results
function render(res) {
  last = res;
  const g = res.genome, fit = res.ancestry_fit;
  $("pr-error").hidden = true;

  const warn = $("pr-warnings");
  warn.replaceChildren();
  for (const w of res.warnings) {
    const box = el("div", "pr-card pr-msg " + (res.valid ? "warning" : "critical"));
    box.append(el("span", "pr-ico", "!"));
    const body = el("div");
    body.append(el("b", "", res.valid ? "Treat these numbers as rough" : "This genome cannot be scored"), el("span", "", w));
    box.append(body);
    warn.append(box);
  }

  const facts = $("pr-facts");
  facts.replaceChildren();
  const fact = (k, v) => { facts.append(el("dt", "", k)); const dd = el("dd"); if (typeof v === "string") dd.textContent = v; else dd.append(v); facts.append(dd); };
  fact("File", g.file);
  fact("Format", (FORMATS[g.format] || g.format) + ", GRCh" + g.build);
  fact("Genotypes", num(g.called) + " autosomal SNPs with a call");
  const level = !(fit.value <= fit.invalid_above) ? "critical" : (fit.value <= fit.good_below ? "good" : "warning");
  const badge = el("span", "pr-status " + level);
  badge.append(el("i", "", level === "good" ? "✓" : "!"),
    document.createTextNode({ good: "Matches the European reference", warning: "Borderline", critical: "Does not match the European reference" }[level]));
  const fitLine = el("span");
  fitLine.append(badge, document.createTextNode(Number.isFinite(fit.value)
    ? `  ·  ${fit.value.toFixed(3)} over ${num(fit.snps)} SNPs (1.00 is typical; above ${fit.invalid_above} is not scorable)` : ""));
  fact("Ancestry check", fitLine);
  fact("Reference data", `${res.reference.name}, built ${res.reference.created}`);

  drawChart(res);
  drawTable(res);

  const src = $("pr-sources");
  src.replaceChildren();
  for (const d of res.disorders) {
    const s = res.sources[d.key];
    if (!s) continue;
    src.append(el("dt", "", d.name));
    src.append(el("dd", "", "Study: " + s.citation + (s.ref ? " (" + s.ref + ")" : "") + "." + (s.notes ? " " + s.notes : "")));
    src.append(el("dd", "", "Baseline: " + s.prevalence_source + "."));
    src.append(el("dd", "", "Accuracy: " + s.accuracy_source));
  }
  $("pr-low-note").hidden = !res.disorders.some((d) => d.low_power);
  $("pr-results").hidden = false;
  $("pr-results").scrollIntoView({ block: "start" });
}

const TICKS = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1];

function drawChart(res) {
  const rows = $("pr-rows"), ticksEl = $("pr-ticks");
  rows.replaceChildren();
  ticksEl.replaceChildren();
  $("pr-legend-estimate").hidden = !res.valid;
  const vals = [];
  for (const d of res.disorders) {
    vals.push(d.prevalence);
    if (res.valid) vals.push(d.probability);
  }
  if (!vals.length) return;
  // the axis runs between the tick marks that enclose the data
  const lo = Math.max(...TICKS.filter((t) => t <= Math.min(...vals) * 1.0001), TICKS[0]);
  const hi = Math.min(...TICKS.filter((t) => t >= Math.max(...vals) * 0.9999), 1);
  const x = (p) => 100 * (Math.log10(Math.min(hi, Math.max(lo, p))) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo));
  const ticks = TICKS.filter((t) => t >= lo && t <= hi);
  const every = ticks.length > 8 ? 2 : 1;
  ticks.forEach((t, i) => {
    if (i % every && i !== ticks.length - 1) return;
    const s = el("span", "", (100 * t >= 1 ? (100 * t).toFixed(0) : String(+(100 * t).toPrecision(1))) + "%");
    s.style.left = x(t) + "%";
    ticksEl.append(s);
  });
  for (const d of res.disorders) {
    const shown = res.valid;
    const row = el("div", "pr-row");
    row.tabIndex = 0;
    const name = el("div", "pr-name", d.name);
    if (d.low_power) name.append(el("small", "", "little information"));
    const track = el("div", "pr-track");
    for (const t of ticks) { const gl = el("div", "pr-g"); gl.style.left = x(t) + "%"; track.append(gl); }
    if (shown) {
      const a = x(d.prevalence), b = x(d.probability);
      const link = el("div", "pr-link");
      link.style.left = Math.min(a, b) + "%";
      link.style.width = Math.abs(a - b) + "%";
      track.append(link);
    }
    const pb = el("div", "pr-pt b"); pb.style.left = x(d.prevalence) + "%"; track.append(pb);
    if (shown) { const pe = el("div", "pr-pt e"); pe.style.left = x(d.probability) + "%"; track.append(pe); }
    const fig = el("div", "pr-fig", shown ? pct(d.probability) : "–");
    if (shown) fig.append(document.createTextNode(" "), el("small", "", d.relative_risk.toFixed(2) + "×"));
    row.append(name, track, fig);
    row.setAttribute("aria-label", d.name + ": baseline " + pct(d.prevalence) +
      (shown ? ", your estimate " + pct(d.probability) + ", " + ordinal(d.percentile) + " percentile" : ", no estimate"));
    const tipFor = (ev) => showTip(row, d, shown, ev);
    row.addEventListener("pointermove", tipFor);
    row.addEventListener("focus", tipFor);
    row.addEventListener("pointerleave", hideTip);
    row.addEventListener("blur", hideTip);
    rows.append(row);
  }
}

function showTip(row, d, shown, ev) {
  const tip = $("pr-tip"), chart = $("pr-chart");
  tip.replaceChildren(el("b", "", d.name));
  const line = (value, label, key) => {
    const r = el("div");
    const lab = el("span");
    if (key) lab.append(el("i", key));
    lab.append(document.createTextNode(label));
    r.append(el("span", "", value), lab);
    tip.append(r);
  };
  if (shown) line(pct(d.probability), "Your estimate", "e");
  line(pct(d.prevalence), "Population baseline", "b");
  if (shown) {
    line(d.relative_risk.toFixed(2) + "×", "Relative to average");
    line(ordinal(d.percentile), "Score percentile");
  }
  line(num(d.snps_used), "SNPs used");
  tip.hidden = false;
  const box = chart.getBoundingClientRect(), r = row.getBoundingClientRect();
  const px = ev && ev.clientX ? ev.clientX - box.left : r.left - box.left + r.width / 2;
  tip.style.left = Math.max(0, Math.min(box.width - tip.offsetWidth, px + 14)) + "px";
  const below = r.bottom - box.top + 6;
  tip.style.top = (below + tip.offsetHeight > box.height + 40 ? r.top - box.top - tip.offsetHeight - 6 : below) + "px";
}
function hideTip() { $("pr-tip").hidden = true; }

function drawTable(res) {
  const body = $("pr-tbody");
  body.replaceChildren();
  for (const d of res.disorders) {
    const shown = res.valid;
    const tr = el("tr");
    const first = el("td", "", d.name);
    if (d.low_power) first.append(el("small", "pr-tag", "little information"));
    for (const w of d.warnings) first.append(el("small", "pr-tag", w));
    tr.append(first,
      el("td", "", pct(d.prevalence)),
      el("td", "", shown ? pct(d.probability) : "–"),
      el("td", "", shown ? d.relative_risk.toFixed(2) + "×" : "–"),
      el("td", "", shown ? ordinal(d.percentile) : "–"),
      el("td", "", num(d.snps_used)));
    const acc = el("td", "", (100 * d.r2_liability).toFixed(d.r2_liability < 0.01 ? 2 : 1) + "% ");
    if (d.r2_basis) acc.append(el("small", "", BASIS[d.r2_basis] || d.r2_basis));
    tr.append(acc);
    body.append(tr);
  }
}

// ---------------------------------------------------------------- actions
function asText(res) {
  const pad = (s, n) => String(s).padStart(n);
  const lines = [`pgc-risk - ${res.genome.file}`,
    `Ancestry fit ${Number.isFinite(res.ancestry_fit.value) ? res.ancestry_fit.value.toFixed(3) : "n/a"} over ${num(res.ancestry_fit.snps)} SNPs`,
    ...res.warnings.map((w) => "WARNING: " + w), "",
    "Disorder".padEnd(34) + pad("Baseline", 10) + pad("Estimate", 10) + pad("Relative", 10) + pad("Percentile", 12) + pad("SNPs", 8)];
  for (const d of res.disorders) {
    lines.push((d.name + (d.low_power ? " *" : "")).padEnd(34) + pad(pct(d.prevalence), 10) +
      pad(res.valid ? pct(d.probability) : "-", 10) + pad(res.valid ? d.relative_risk.toFixed(2) + "x" : "-", 10) +
      pad(res.valid ? ordinal(d.percentile) : "-", 12) + pad(num(d.snps_used), 8));
  }
  lines.push("", "* little information: the study behind this score is small.",
             "A statistical expectation for people with a similar score, not a diagnosis.");
  return lines.join("\n");
}
$("pr-again").addEventListener("click", () => { last = null; $("pr-results").hidden = true; $("pr-drop").scrollIntoView({ block: "center" }); });
$("pr-print").addEventListener("click", () => window.print());
$("pr-copy").addEventListener("click", async (e) => {
  try { await navigator.clipboard.writeText(asText(last)); e.target.textContent = "Copied"; }
  catch (err) { e.target.textContent = "Copy failed"; }
  setTimeout(() => { e.target.textContent = "Copy as text"; }, 1500);
});
$("pr-save").addEventListener("click", () => {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([JSON.stringify(last, null, 2)], { type: "application/json" }));
  a.download = "pgc-risk-results.json";
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
});
