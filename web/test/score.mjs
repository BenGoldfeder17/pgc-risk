// Run the browser scorer under Node: node web/test/score.mjs <bundle.pgcr> <genome file> [...]
// Prints one JSON line per genome file (or {"error": ...}), then a line of normalCdf values
// for tools/check_web.py to compare with Python.
import { readFileSync } from "node:fs";
import { basename } from "node:path";
import { PgcRiskError, normalCdf, readBundle, scoreGenomeBytes, unwrap } from "../pgcrisk-core.js";

const [bundlePath, ...genomes] = process.argv.slice(2);
const bundle = readBundle(await unwrap(new Uint8Array(readFileSync(bundlePath))));
for (const g of genomes) {
  try {
    const t0 = performance.now();
    const result = await scoreGenomeBytes(bundle, new Uint8Array(readFileSync(g)), basename(g));
    result.ms = Math.round(performance.now() - t0);
    console.log(JSON.stringify(result));
  } catch (e) {
    if (!(e instanceof PgcRiskError)) throw e;
    console.log(JSON.stringify({ error: e.message }));
  }
}
const grid = [];
for (let z = -9; z <= 9.0001; z += 0.125) grid.push(normalCdf(z));
console.log(JSON.stringify({ normalCdf: grid }));
