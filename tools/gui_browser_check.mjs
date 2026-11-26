// Drive the real page in a real (headless) browser: choose a file, wait for the result, report.
//
// Starts Microsoft Edge headless with remote debugging, opens the running GUI, puts a genome
// file into the page's file input (as a user's file dialog would), and waits for the results.
// Fails if the page logs an error, if the browser blocks anything under the page's content
// security policy, or if the page asks the network for anything but this server.
//
// usage: node tools/gui_browser_check.mjs <genome file> [port or URL] [screenshot.png] [access code]
//        (pgc-risk gui --no-browser must be running; needs Node 22+ for the built-in WebSocket)
// With a URL such as https://192.168.1.20:8765/ and an access code it plays "another device":
// it accepts the self-signed certificate and types the code into the access-code page.
import { spawn } from "node:child_process";
import { mkdtempSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const file = resolve(process.argv[2]);
const where = process.argv[3] || "8765";
const shot = process.argv[4] && process.argv[4] !== "-" ? process.argv[4] : null;
const code = process.argv[5];
const base = where.includes("://") ? where : `http://127.0.0.1:${where}/`;
const EDGE = ["C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe",
              "C:/Program Files/Microsoft/Edge/Application/msedge.exe"];
const debugPort = 9300 + Math.floor(Math.random() * 500);
const profile = mkdtempSync(join(tmpdir(), "pgc-risk-edge-"));
const { existsSync } = await import("node:fs");
const edge = spawn(EDGE.find((p) => existsSync(p)), ["--headless=new", "--disable-gpu",
  `--remote-debugging-port=${debugPort}`, `--user-data-dir=${profile}`, "--window-size=1100,1400",
  ...(base.startsWith("https") ? ["--ignore-certificate-errors"] : []), "about:blank"],
  { stdio: "ignore" });

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
let target;
for (let i = 0; i < 100 && !target; i++) {
  await sleep(100);
  try {
    const list = await (await fetch(`http://127.0.0.1:${debugPort}/json`)).json();
    target = list.find((t) => t.type === "page");
  } catch { /* browser still starting */ }
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
  if (msg.method === "Log.entryAdded" && ["error", "warning"].includes(msg.params.entry.level))
    problems.push(`${msg.params.entry.level}: ${msg.params.entry.text}`);
  if (msg.method === "Runtime.consoleAPICalled" && msg.params.type === "error")
    problems.push("console.error: " + msg.params.args.map((a) => a.value).join(" "));
  if (msg.method === "Network.requestWillBeSent") requests.push(msg.params.request.url);
});
const send = (method, params = {}) => new Promise((res) => {
  const id = ++seq;
  pending.set(id, res);
  ws.send(JSON.stringify({ id, method, params }));
});
const evaluate = async (expression) =>
  (await send("Runtime.evaluate", { expression, returnByValue: true, awaitPromise: true })).result.result.value;

let ok = false;
try {
  for (const d of ["Runtime", "Log", "Network", "Page", "DOM"]) await send(d + ".enable");
  await send("Page.navigate", { url: base });
  await sleep(800);
  if (await evaluate(`!!document.getElementById("code")`)) {
    console.log("access-code page shown:", await evaluate(`document.querySelector("p").textContent`));
    await evaluate(`(() => { document.getElementById("code").value = "000000";
                             document.getElementById("form").requestSubmit(); })()`);
    await sleep(600);
    console.log("wrong code ->", await evaluate(`document.getElementById("msg").textContent`));
    await evaluate(`(() => { document.getElementById("code").value = ${JSON.stringify(code || "")};
                             document.getElementById("form").requestSubmit(); })()`);
    await sleep(1500);
    problems.length = 0;      // the wrong code above is answered with 403 on purpose
    console.log("after the right code, upload page shown:", await evaluate(`!!document.getElementById("drop")`),
                "| note:", await evaluate(`document.querySelector(".private span") ? document.querySelector(".private span").textContent.slice(0, 150) : ""`));
  }
  const doc = await send("DOM.getDocument");
  const input = await send("DOM.querySelector", { nodeId: doc.result.root.nodeId, selector: "#file" });
  await send("DOM.setFileInputFiles", { nodeId: input.result.nodeId, files: [file] });
  await evaluate(`document.getElementById("file").dispatchEvent(new Event("change"))`);
  const stages = new Set();
  let state = "";
  for (let i = 0; i < 1500; i++) {
    state = await evaluate(`(() => {
      const vis = (id) => !document.getElementById(id).hidden;
      return vis("results") ? "results" : vis("error") ? "error" : "stage:" + document.getElementById("stage").textContent;
    })()`);
    if (state === "results" || state === "error") break;
    stages.add(state);
    await sleep(100);
  }
  console.log("progress stages shown:", [...stages].join(" | "));
  if (state === "error") {
    console.log("page shows error:", await evaluate(`document.getElementById("error-text").textContent`));
  } else if (state === "results") {
    const rows = await evaluate(`[...document.querySelectorAll("#tbody tr")].map((tr) =>
      [...tr.children].map((td) => td.textContent.trim()).join(" | ")).join("\\n")`);
    console.log(rows);
    console.log("chart rows:", await evaluate(`document.querySelectorAll("#rows .row").length`),
                "| estimate markers:", await evaluate(`document.querySelectorAll("#rows .pt.e").length`),
                "| facts:", await evaluate(`[...document.querySelectorAll("#facts dd")].map((d) => d.textContent).join(" / ")`));
    // hover a row: the tooltip must appear with that row's values
    await evaluate(`document.querySelector("#rows .row").dispatchEvent(new PointerEvent("pointermove", {clientX: 500, clientY: 300, bubbles: true}))`);
    console.log("tooltip:", await evaluate(`document.getElementById("tip").hidden ? "(hidden)" : document.getElementById("tip").textContent`));
    ok = true;
  }
  if (shot) {
    const metrics = await send("Page.getLayoutMetrics");
    const h = Math.ceil(metrics.result.cssContentSize.height);
    await send("Emulation.setDeviceMetricsOverride", { width: 1100, height: h, deviceScaleFactor: 1, mobile: false });
    const png = await send("Page.captureScreenshot", { format: "png" });
    writeFileSync(shot, Buffer.from(png.result.data, "base64"));
    console.log("screenshot:", shot);
  }
  const outside = requests.filter((u) => !u.startsWith(base) && !u.startsWith("data:") && !u.startsWith("blob:") && u !== "about:blank");
  console.log("requests made by the page:", requests.length, "| to anywhere but this server:", outside.length, outside.join(" "));
  console.log("errors or policy violations reported by the browser:", problems.length, problems.join(" ; "));
  ok = ok && outside.length === 0 && problems.length === 0;
} finally {
  ws.close();
  edge.kill();
  await sleep(500);
  try { rmSync(profile, { recursive: true, force: true }); } catch { /* browser still letting go */ }
}
console.log(ok ? "PASS" : "FAIL");
process.exit(ok ? 0 : 1);
