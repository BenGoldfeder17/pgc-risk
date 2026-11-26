// The browser tool in a real (headless) browser, including the promise that nothing is uploaded.
//
// Serves a folder over http on this machine, opens the page in Microsoft Edge, then:
//   1. takes the browser OFFLINE before any file is chosen,
//   2. hands the page a genome file first (it must ask for the reference data) and then the
//      reference bundle (it must then score the waiting genome by itself),
//   3. compares every figure in the page's table with the scorer run directly under Node,
//   4. goes back online, reloads, and checks the reference data is still there (kept in the
//      browser's own storage) and a genome alone is enough; then removes it again.
// It fails if the page logs an error, if the browser reports a blocked resource, or if the page
// made any request at all after its own files had loaded, or any request to another host.
//
// usage: node tools/web_browser_check.mjs <folder to serve, or https://site> <page, e.g. index.html> <bundle> <genome> [screenshot.png]
import { spawn } from "node:child_process";
import { existsSync, mkdtempSync, readFileSync, rmSync, statSync, writeFileSync } from "node:fs";
import { createServer } from "node:http";
import { tmpdir } from "node:os";
import { basename, extname, join, resolve } from "node:path";
import { readBundle, scoreGenomeBytes, unwrap } from "../web/pgcrisk-core.js";

const [root, pageArg, bundlePath, genomePath, shot] = process.argv.slice(2);
const pagePath = "/" + pageArg.replace(/^\/+/, "");      // given without the leading slash (shells rewrite "/x")
const remote = /^https?:/.test(root);                      // a deployed site instead of a local folder
const port = 8700 + Math.floor(Math.random() * 200);
const origin = remote ? root.replace(/[/]+$/, "") : `http://127.0.0.1:${port}`;
// a plain static file server, standing in for the web host
const MIME = { ".html": "text/html; charset=utf-8", ".js": "text/javascript", ".mjs": "text/javascript", ".css": "text/css",
               ".svg": "image/svg+xml", ".woff2": "font/woff2", ".woff": "font/woff", ".webp": "image/webp", ".json": "application/json" };
const server = createServer((req, res) => {
  let file = join(resolve(root), decodeURIComponent(new URL(req.url, origin).pathname).toLowerCase());
  if (existsSync(file) && statSync(file).isDirectory()) file = join(file, "index.html");
  if (!file.startsWith(resolve(root)) || !existsSync(file)) { res.writeHead(404); res.end("not found"); return; }
  res.writeHead(200, { "Content-Type": MIME[extname(file)] || "application/octet-stream" });
  res.end(readFileSync(file));
});
if (!remote) await new Promise((r) => server.listen(port, "127.0.0.1", r));
const EDGE = ["C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe", "C:/Program Files/Microsoft/Edge/Application/msedge.exe"];
const debugPort = 9300 + Math.floor(Math.random() * 500);
const profile = mkdtempSync(join(tmpdir(), "pgc-risk-edge-"));
const edge = spawn(EDGE.find((p) => existsSync(p)), ["--headless=new", "--disable-gpu", `--remote-debugging-port=${debugPort}`,
  `--user-data-dir=${profile}`, "--window-size=1200,1500", "about:blank"], { stdio: "ignore" });

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const pct = (x) => { const v = 100 * x; return v < 0.095 ? v.toFixed(2) + "%" : v < 9.95 ? v.toFixed(1) + "%" : Math.floor(v + 0.5) + "%"; };
const ordinal = (p) => { const v = Math.min(99, Math.max(1, Math.round(100 * p)));
  return v + ((v % 100 >= 10 && v % 100 <= 20) ? "th" : ({ 1: "st", 2: "nd", 3: "rd" }[v % 10] || "th")); };

let target;
for (let i = 0; i < 100 && !target; i++) {
  await sleep(100);
  try { target = (await (await fetch(`http://127.0.0.1:${debugPort}/json`)).json()).find((t) => t.type === "page"); } catch { /* starting */ }
}
const ws = new WebSocket(target.webSocketDebuggerUrl);
await new Promise((r) => ws.addEventListener("open", r));
let seq = 0;
const pending = new Map();
const problems = [], requests = [];
ws.addEventListener("message", (m) => {
  const msg = JSON.parse(m.data);
  if (msg.id && pending.has(msg.id)) { pending.get(msg.id)(msg); pending.delete(msg.id); return; }
  if (msg.method === "Runtime.exceptionThrown") problems.push("exception: " + msg.params.exceptionDetails.text);
  if (msg.method === "Log.entryAdded" && ["error", "warning"].includes(msg.params.entry.level)) problems.push(`${msg.params.entry.level}: ${msg.params.entry.text}`);
  if (msg.method === "Runtime.consoleAPICalled" && msg.params.type === "error") problems.push("console.error: " + msg.params.args.map((a) => a.value).join(" "));
  if (msg.method === "Network.requestWillBeSent") requests.push({ url: msg.params.request.url, method: msg.params.request.method, body: !!msg.params.request.hasPostData });
});
const send = (method, params = {}) => new Promise((res) => { const id = ++seq; pending.set(id, res); ws.send(JSON.stringify({ id, method, params })); });
const evaluate = async (expression) => (await send("Runtime.evaluate", { expression, returnByValue: true, awaitPromise: true })).result.result.value;
const visible = (id) => evaluate(`!document.getElementById("${id}").hidden`);
const give = async (files) => {
  const doc = await send("DOM.getDocument");
  const input = await send("DOM.querySelector", { nodeId: doc.result.root.nodeId, selector: "#pr-file" });
  await send("DOM.setFileInputFiles", { nodeId: input.result.nodeId, files: files.map((f) => resolve(f)) });
  await evaluate(`document.getElementById("pr-file").dispatchEvent(new Event("change"))`);
};
const until = async (what, fn, ms = 120000) => {
  for (let t = 0; t < ms; t += 150) { if (await fn()) return true; await sleep(150); }
  throw new Error("timed out waiting for " + what);
};
const table = () => evaluate(`[...document.querySelectorAll("#pr-tbody tr")].map((tr) => [...tr.children].map((td) => td.textContent.trim()))`);

let ok = true;
const check = (label, pass, detail = "") => { console.log(`${pass ? "ok  " : "FAIL"} ${label}${detail ? " - " + detail : ""}`); ok = ok && pass; };
try {
  for (const d of ["Runtime", "Log", "Network", "Page", "DOM"]) await send(d + ".enable");
  await send("Page.navigate", { url: origin + pagePath });
  await until("the page", () => evaluate(`document.documentElement.dataset.prReady === "1" && document.getElementById("pr-data-status").textContent.includes("Not loaded")`));
  await sleep(700);                                          // let late page assets (fonts, icon) finish
  const note = await evaluate(`document.querySelector(".pr-privacy").textContent.replace(/\\s+/g, " ").trim()`);
  check("page states that the data is not uploaded", /not uploaded/i.test(note) && /own\s+device/i.test(note), note.slice(0, 96) + "...");
  const loaded = requests.length;
  check("page loaded only its own files", requests.every((r) => r.url.startsWith(origin) && r.method === "GET"),
        `${loaded} requests: ${[...new Set(requests.map((r) => basename(new URL(r.url).pathname) || "/"))].join(", ")}`);

  await send("Network.emulateNetworkConditions", { offline: true, latency: 0, downloadThroughput: -1, uploadThroughput: -1 });
  await give([genomePath]);
  await until("the request for reference data", () => visible("pr-error"));
  check("genome before reference data: page asks for the reference data", (await evaluate(`document.getElementById("pr-error-text").textContent`)).includes("reference data"));
  await give([bundlePath]);
  await until("results", () => visible("pr-results"));
  check("scored offline once the reference data arrived", true, await evaluate(`document.getElementById("pr-data-status").textContent.trim().slice(0, 110)`));

  const bundle = readBundle(await unwrap(new Uint8Array(readFileSync(bundlePath))));
  const want = await scoreGenomeBytes(bundle, new Uint8Array(readFileSync(genomePath)), basename(genomePath));
  const rows = await table();
  let same = rows.length === want.disorders.length;
  want.disorders.forEach((d, i) => {
    const r = rows[i] || [];
    const exp = want.valid ? [pct(d.prevalence), pct(d.probability), d.relative_risk.toFixed(2) + "×", ordinal(d.percentile), d.snps_used.toLocaleString("en-US")]
                           : [pct(d.prevalence), "–", "–", "–", d.snps_used.toLocaleString("en-US")];
    if (!r[0] || !r[0].startsWith(d.name) || exp.some((v, k) => r[k + 1] !== v)) { same = false; console.log("   mismatch:", r.join(" | "), "expected", exp.join(" | ")); }
  });
  check(`table equals the scorer run directly (${rows.length} disorders)`, same, rows.slice(0, 3).map((r) => r.slice(0, 5).join(" ")).join(" ; "));
  check("chart drawn", (await evaluate(`document.querySelectorAll("#pr-rows .pr-row").length`)) === rows.length &&
        (await evaluate(`document.querySelectorAll("#pr-rows .pr-pt.e").length`)) === (want.valid ? rows.length : 0));
  await evaluate(`document.querySelector("#pr-rows .pr-row").dispatchEvent(new PointerEvent("pointermove", {clientX: 600, clientY: 400, bubbles: true}))`);
  check("tooltip shows on hover", !(await evaluate(`document.getElementById("pr-tip").hidden`)));
  await evaluate(`document.getElementById("pr-save").click()`);
  await sleep(300);
  if (shot) {
    const metrics = await send("Page.getLayoutMetrics");
    await send("Emulation.setDeviceMetricsOverride", { width: 1200, height: Math.ceil(metrics.result.cssContentSize.height), deviceScaleFactor: 1, mobile: false });
    writeFileSync(shot, Buffer.from((await send("Page.captureScreenshot", { format: "png" })).result.data, "base64"));
    await send("Emulation.clearDeviceMetricsOverride");
  }
  const after = requests.slice(loaded).filter((r) => !r.url.startsWith("blob:") && !r.url.startsWith("data:"));
  check("no network request of any kind after the files were chosen", after.length === 0, after.map((r) => r.method + " " + r.url).join(" "));

  await send("Network.emulateNetworkConditions", { offline: false, latency: 0, downloadThroughput: -1, uploadThroughput: -1 });
  await send("Page.reload");
  await until("the remembered reference data", () => evaluate(`document.documentElement.dataset.prReady === "1" && document.getElementById("pr-data-status").textContent.includes("Loaded")`));
  check("after a reload the reference data is still in the browser", true);
  await give([genomePath]);
  await until("results after reload", () => visible("pr-results"));
  check("genome alone is enough the second time", (await table()).length === rows.length);
  await evaluate(`document.getElementById("pr-forget").click()`);
  await until("removal", () => evaluate(`document.getElementById("pr-data-status").textContent.includes("Not loaded")`));
  check("reference data can be removed from the browser", true);
  check("no request to any other host, and none that carries data",
        requests.every((r) => (r.url.startsWith(origin) || r.url.startsWith("blob:") || r.url.startsWith("data:")) && r.method === "GET" && !r.body));
  check("no errors or blocked resources reported by the browser", problems.length === 0, problems.join(" ; "));
} catch (e) {
  check("run completed", false, String(e && e.message ? e.message : e));
  try {
    console.log("   page:", await evaluate(`location.href + " | " + document.title + " | " + (document.body ? document.body.innerText.slice(0, 160).replace(/\s+/g, " ") : "(no body)")`));
    console.log("   requests:", requests.map((r) => r.url).join(" "));
    console.log("   problems:", problems.join(" ; "));
  } catch { /* browser gone */ }
} finally {
  ws.close();
  edge.kill();
  if (!remote) server.close();
  await sleep(500);
  try { rmSync(profile, { recursive: true, force: true }); } catch { /* browser still letting go */ }
}
console.log(ok ? "PASS" : "FAIL");
process.exit(ok ? 0 : 1);
