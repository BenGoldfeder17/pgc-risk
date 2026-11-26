// Background thread of the pgc-risk page: reads the files it is handed and scores, on this device.
// It keeps the reference data in this browser's own storage (IndexedDB) so it only has to be
// loaded once; a genome file is never stored.
import { PgcRiskError, isBundle, readBundle, readGenome, scoreCounts, unwrap } from "./pgcrisk-core.js";

let bundle = null;
let waiting = null;             // a genome file that arrived before the reference data

// ---------------------------------------------------------------- browser storage (reference data only)
function store() {
  return new Promise((resolve, reject) => {
    const req = indexedDB.open("pgc-risk", 1);
    req.onupgradeneeded = () => req.result.createObjectStore("kv");
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
}
async function stored(mode, fn) {
  const db = await store();
  try {
    return await new Promise((resolve, reject) => {
      const tx = db.transaction("kv", mode);
      const req = fn(tx.objectStore("kv"));
      tx.oncomplete = () => resolve(req.result);
      tx.onerror = tx.onabort = () => reject(tx.error);
    });
  } finally {
    db.close();
  }
}

function describe(kept) {
  const h = bundle.header;
  return { type: "bundle", kept, info: { name: h.name, created: h.created, snps: h.snps, disorders: h.disorders.length,
                                         description: h.description } };
}

async function useBundle(bytes, blob) {
  bundle = readBundle(bytes);
  let kept = false;
  if (blob) {
    try { await stored("readwrite", (s) => s.put(blob, "bundle")); kept = true; } catch (e) { /* private window */ }
  }
  postMessage(describe(kept || !blob));
  if (waiting) { const file = waiting; waiting = null; await score(file); }
}

async function score(file) {
  postMessage({ type: "progress", text: "Reading " + file.name });
  const plain = await unwrap(new Uint8Array(await file.arrayBuffer()));
  if (isBundle(plain)) return useBundle(plain, file);
  if (!bundle) { waiting = file; postMessage({ type: "need-bundle", name: file.name }); return; }
  const { counts, genome } = readGenome(new TextDecoder("utf-8").decode(plain), bundle, file.name);
  postMessage({ type: "progress", text: "Scoring" });
  const result = scoreCounts(bundle, counts);
  result.genome = genome;
  result.reference = { name: bundle.header.name, created: bundle.header.created, tool: bundle.header.tool };
  result.sources = {};
  for (const d of bundle.header.disorders) {
    result.sources[d.key] = { citation: d.citation, ref: d.ref, notes: d.notes,
                              prevalence_source: d.prevalence_source, accuracy_source: d.r2_source };
  }
  postMessage({ type: "result", result });
}

onmessage = async (e) => {
  const msg = e.data;
  try {
    if (msg.cmd === "restore") {
      let blob = null;
      try { blob = await stored("readonly", (s) => s.get("bundle")); } catch (err) { /* storage unavailable */ }
      if (!blob) return postMessage({ type: "no-bundle" });
      await useBundle(await unwrap(new Uint8Array(await blob.arrayBuffer())), null);
    } else if (msg.cmd === "forget") {
      bundle = null;
      waiting = null;
      try { await stored("readwrite", (s) => s.delete("bundle")); } catch (err) { /* nothing stored */ }
      postMessage({ type: "no-bundle" });
    } else if (msg.cmd === "file") {
      await score(msg.file);
    }
  } catch (err) {
    postMessage({ type: "error", message: err instanceof PgcRiskError ? err.message
      : "Something went wrong reading that file (" + (err && err.message ? err.message : err) + ")." });
  }
};
