// pgc-risk in the browser: read a reference bundle and a genome file, score, all on this device.
//
// Nothing in this file touches the network.  It mirrors pgc_risk/bundle.py and
// pgc_risk/genome.py; web/test/compare.mjs checks that both give the same numbers.

const MAGIC = "PGCRISK1";
const FORMAT = 1;
const BASES = { A: 0, C: 1, G: 2, T: 3, a: 0, c: 1, g: 2, t: 3 };

export class PgcRiskError extends Error {}

// ---------------------------------------------------------------------------- bytes
async function inflate(bytes, kind) {
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream(kind));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

/** The largest file inside a zip archive (23andMe and AncestryDNA downloads are zips). */
async function unzipLargest(bytes) {
  const v = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let eocd = -1;
  for (let i = bytes.length - 22; i >= Math.max(0, bytes.length - 65557); i--) {
    if (v.getUint32(i, true) === 0x06054b50) { eocd = i; break; }
  }
  if (eocd < 0) throw new PgcRiskError("The zip archive could not be read.");
  const entries = v.getUint16(eocd + 10, true);
  let p = v.getUint32(eocd + 16, true);
  let best = null;
  for (let n = 0; n < entries && p + 46 <= bytes.length && v.getUint32(p, true) === 0x02014b50; n++) {
    const entry = { method: v.getUint16(p + 10, true), packed: v.getUint32(p + 20, true),
                    size: v.getUint32(p + 24, true), local: v.getUint32(p + 42, true) };
    const nameLen = v.getUint16(p + 28, true);
    const isFolder = bytes[p + 46 + nameLen - 1] === 47;
    if (!isFolder && (!best || entry.size > best.size)) best = entry;
    p += 46 + nameLen + v.getUint16(p + 30, true) + v.getUint16(p + 32, true);
  }
  if (!best) throw new PgcRiskError("The zip archive is empty.");
  const start = best.local + 30 + v.getUint16(best.local + 26, true) + v.getUint16(best.local + 28, true);
  const data = bytes.subarray(start, start + best.packed);
  if (best.method === 0) return data;
  if (best.method === 8) return inflate(data, "deflate-raw");
  throw new PgcRiskError("The zip archive uses a compression method this page cannot read.");
}

/** Undo gzip / zip wrapping, whatever the file is called. */
export async function unwrap(bytes) {
  if (bytes.length > 2 && bytes[0] === 0x1f && bytes[1] === 0x8b) return inflate(bytes, "gzip");
  if (bytes.length > 4 && bytes[0] === 0x50 && bytes[1] === 0x4b && bytes[2] === 3 && bytes[3] === 4) return unzipLargest(bytes);
  return bytes;
}

export function isBundle(bytes) {
  if (bytes.length < 12) return false;
  for (let i = 0; i < 8; i++) if (bytes[i] !== MAGIC.charCodeAt(i)) return false;
  return true;
}

// ---------------------------------------------------------------------------- bundle
const TYPES = { u1: Uint8Array, u2: Uint16Array, u4: Uint32Array, i2: Int16Array, i4: Int32Array };

/** Parse an unwrapped bundle (see pgc_risk/bundle.py for the layout). */
export function readBundle(bytes) {
  if (!isBundle(bytes)) throw new PgcRiskError("This is not a pgc-risk reference data file.");
  // typed arrays need an aligned buffer of their own
  const buf = bytes.byteOffset === 0 && bytes.byteLength === bytes.buffer.byteLength
    ? bytes.buffer : bytes.slice().buffer;
  const headLen = new DataView(buf).getUint32(8, true);
  const header = JSON.parse(new TextDecoder().decode(new Uint8Array(buf, 12, headLen)));
  if (header.format !== FORMAT) {
    throw new PgcRiskError(`This reference data file is format ${header.format}; this page reads format ${FORMAT}. ` +
                           "Download the file that goes with this version of the page.");
  }
  const base = 12 + headLen;
  const arr = (name) => {
    const s = header.sections[name];
    return new TYPES[s.type](buf, base + s.offset, s.length);
  };
  const rsid = arr("snp.rsid");
  const bundle = { header, rsid, bAllele: arr("snp.b"), freqRaw: arr("snp.freq"), disorders: [] };
  for (const meta of header.disorders) {
    const delta = arr(`${meta.key}.idx_delta`);
    const idx = new Uint32Array(delta.length);
    let at = 0;
    for (let i = 0; i < delta.length; i++) { at += delta[i]; idx[i] = at; }
    bundle.disorders.push({ meta, idx, beta: arr(`${meta.key}.beta`), ld: arr(`${meta.key}.ld`),
                            clumpLen: arr(`${meta.key}.clump_len`) });
  }
  // rs number -> SNP index (+1), open addressing
  let bits = 1;
  while ((1 << bits) < rsid.length * 2) bits++;
  const table = new Int32Array(1 << bits);
  const mask = table.length - 1;
  for (let i = 0; i < rsid.length; i++) {
    let h = (Math.imul(rsid[i], 0x9e3779b1) >>> (32 - bits)) & mask;
    while (table[h] !== 0) h = (h + 1) & mask;
    table[h] = i + 1;
  }
  bundle.lookup = (rs) => {
    let h = (Math.imul(rs, 0x9e3779b1) >>> (32 - bits)) & mask;
    for (let t = table[h]; t !== 0; t = table[h = (h + 1) & mask]) if (rsid[t - 1] === rs) return t - 1;
    return -1;
  };
  return bundle;
}

// ---------------------------------------------------------------------------- genome
const UNSEEN = -1, AMBIGUOUS = -2, NO_CALL = -3;
const AUTOSOME = new Set();
for (let c = 1; c <= 22; c++) { AUTOSOME.add(String(c)); AUTOSOME.add("chr" + c); }

function buildFromText(text) {
  const t = text.toLowerCase();
  for (const [build, marks] of [[38, ["grch38", "hg38", "build 38", "b38", "hs38", "_v38"]],
                                [37, ["grch37", "hg19", "build 37", "b37", "hs37", "_v37"]],
                                [36, ["ncbi36", "hg18", "build 36", "b36"]]]) {
    if (marks.some((m) => t.includes(m))) return build;
  }
  return null;
}

/**
 * Read an array export (23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA, LivingDNA) and count,
 * for every SNP of the bundle, the copies of its A allele.  Returns the counts (Int8Array, a
 * negative value where the file has no usable call) and what was learnt about the file.
 */
export function readGenome(text, bundle, name = "") {
  const first = text.slice(0, 200);
  if (first.startsWith("##fileformat=VCF") || /^#CHROM\t/m.test(text.slice(0, 200000))) {
    throw new PgcRiskError(
      "This is a VCF. The browser version reads genotyping-chip exports only (23andMe, AncestryDNA, " +
      "MyHeritage, FamilyTreeDNA, LivingDNA). Sequencing files need the reference genome to fill in the " +
      "sites they leave out, which the desktop version of pgc-risk does.");
  }
  const counts = new Int8Array(bundle.rsid.length).fill(UNSEEN);
  const { lookup, bAllele } = bundle;
  let fmt = null, delim = "\t", nCols = 0, called = 0, records = 0;
  const comments = [];
  let pos = 0;
  const n = text.length;
  while (pos < n) {
    let end = text.indexOf("\n", pos);
    if (end < 0) end = n;
    let line = text.slice(pos, end);
    pos = end + 1;
    if (line.length === 0) continue;
    if (line.charCodeAt(line.length - 1) === 13) line = line.slice(0, -1);
    if (line.length === 0) continue;
    if (line.charCodeAt(0) === 35) { if (comments.length < 60) comments.push(line); continue; }
    if (fmt === null) {
      const stripped = line.trim();
      if (!stripped) continue;
      delim = stripped.includes("\t") ? "\t" : (stripped.includes(",") ? "," : null);
      const cols = (delim ? stripped.split(delim) : stripped.split(/\s+/)).map((c) => c.trim().replace(/^"|"$/g, ""));
      nCols = cols.length;
      if (nCols !== 4 && nCols !== 5) {
        throw new PgcRiskError(`${name || "This file"}: not a recognised genome file (expected a 23andMe, ` +
                               "AncestryDNA, MyHeritage, FamilyTreeDNA or LivingDNA raw-data export).");
      }
      fmt = nCols === 5 ? "ancestry" : (delim === "," ? "csv" : "23andme");
      if (["rsid", "rs", "snp", "name", "# rsid"].includes(cols[0].toLowerCase())) continue;   // header row
    }
    let parts = delim ? line.split(delim) : line.trim().split(/\s+/);
    if (parts.length !== nCols) continue;
    if (delim === ",") parts = parts.map((c) => c.trim().replace(/^"|"$/g, ""));
    if (!AUTOSOME.has(parts[1])) continue;
    if (!/^-?\d+$/.test(parts[2].trim())) continue;
    let u1 = -1, u2 = -1;
    if (nCols === 5) {
      const x = BASES[parts[3].trim()], y = BASES[parts[4].trim()];
      if (x !== undefined && y !== undefined) { u1 = x; u2 = y; }
    } else {
      const g = parts[3].trim();
      if (g.length === 2) {
        const x = BASES[g[0]], y = BASES[g[1]];
        if (x !== undefined && y !== undefined) { u1 = x; u2 = y; }
      }
    }
    records++;
    if (u1 >= 0) called++;
    const id = parts[0];
    if (id.length < 3 || (id[0] !== "r" && id[0] !== "R") || (id[1] !== "s" && id[1] !== "S")) continue;
    const m = /^\d+/.exec(id.slice(2));
    if (!m) continue;
    const rs = Number(m[0]);
    if (rs >= 4294967296) continue;
    const i = lookup(rs);
    if (i < 0) continue;
    if (counts[i] !== UNSEEN) { counts[i] = AMBIGUOUS; continue; }    // an rs number listed twice is not used
    if (u1 < 0) { counts[i] = NO_CALL; continue; }
    const b = bAllele[i];
    if ((u1 === 0 || u1 === b) && (u2 === 0 || u2 === b)) {
      counts[i] = (u1 === 0) + (u2 === 0);
    } else {
      const f1 = 3 - u1, f2 = 3 - u2;                                 // the same call read on the other strand
      counts[i] = ((f1 === 0 || f1 === b) && (f2 === 0 || f2 === b)) ? (f1 === 0) + (f2 === 0) : NO_CALL;
    }
  }
  if (fmt === null || records === 0) {
    throw new PgcRiskError(`${name || "This file"}: no autosomal SNP genotypes found (expected a 23andMe, ` +
                           "AncestryDNA, MyHeritage, FamilyTreeDNA or LivingDNA raw-data export).");
  }
  return { counts, genome: { file: name, format: fmt, build: buildFromText(comments.join(" ")) || 37, called } };
}

// ---------------------------------------------------------------------------- statistics
/** erfc(x) for x >= 0: a positive-term series below 3, a continued fraction above. */
function erfcPositive(x) {
  if (x < 3) {
    let term = x, sum = x;
    for (let k = 1; k < 200; k++) {
      term *= 2 * x * x / (2 * k + 1);
      sum += term;
      if (term < sum * 1e-17) break;
    }
    return 1 - 2 / Math.sqrt(Math.PI) * Math.exp(-x * x) * sum;
  }
  let f = 0;
  for (let k = 80; k >= 1; k--) f = (k / 2) / (x + f);
  return Math.exp(-x * x) / (Math.sqrt(Math.PI) * (x + f));
}

/** Standard normal distribution function. */
export function normalCdf(z) {
  const half = 0.5 * erfcPositive(Math.abs(z) / Math.SQRT2);
  return z > 0 ? 1 - half : half;
}

/** Lifetime risk for a standardised polygenic score (liability-threshold model). */
export function liabilityRisk(z, threshold, prevalence, r2) {
  if (!(r2 > 0)) return prevalence;
  return 1 - normalCdf((threshold - Math.sqrt(r2) * z) / Math.sqrt(1 - r2));
}

// ---------------------------------------------------------------------------- scoring
/** Score every disorder from A-allele counts; the same computation as bundle.score_counts. */
export function scoreCounts(bundle, counts) {
  const h = bundle.header;
  const { freqRaw } = bundle;
  const scale = h.freq_scale;
  let fitSum = 0, fitSnps = 0;
  for (let i = 0; i < counts.length; i++) {
    const d = counts[i];
    if (d < 0) continue;
    const f = freqRaw[i] / scale;
    if (f < 0.05 || f > 0.95) continue;
    const dev = d - 2 * f;
    fitSum += dev * dev / (2 * f * (1 - f));
    fitSnps++;
  }
  const fit = fitSnps ? fitSum / fitSnps : NaN;
  let valid = true;
  const warnings = [];
  if (fitSnps < h.min_fit_snps) {
    valid = false;
    warnings.push(`Only ${fitSnps.toLocaleString("en-US")} of this file's SNPs are in the reference data - too few to score.`);
  } else if (fit > h.fit_poor) {
    valid = false;
    warnings.push(`Ancestry fit ${fit.toFixed(3)}: this genome does not match the European-ancestry reference the ` +
      "scores are calibrated on (non-European or strongly admixed ancestry). In such genomes the scores sit " +
      "several standard deviations off for reasons unrelated to risk, so no percentiles or probabilities are given.");
  } else if (fit > h.fit_good) {
    warnings.push(`Ancestry fit ${fit.toFixed(3)}: this genome sits some distance from the European-ancestry ` +
      "reference (about 1 European in 20 does; so do partly non-European genomes). Treat the numbers as rough.");
  }
  const disorders = [];
  for (const dis of bundle.disorders) {
    const m = dis.meta;
    const { idx, beta, ld, clumpLen } = dis;
    let raw = 0, mean = 0, indep = 0, cross = 0, used = 0, row = 0;
    for (let c = 0; c < clumpLen.length; c++) {
      const stop = row + clumpLen[c];
      for (let r = row; r < stop; r++) {                       // the clump's best SNP that the file has
        const snp = idx[r];
        const d = counts[snp];
        if (d < 0) continue;
        const b = beta[r] * m.beta_scale, f = freqRaw[snp] / scale;
        const v = 2 * f * (1 - f) * b * b;
        raw += b * d;
        mean += 2 * f * b;
        indep += v;
        cross += v * (ld[r] * h.ld_step - 1);
        used++;
        break;
      }
      row = stop;
    }
    const coverage = m.full_indep > 0 ? Math.min(1, indep / m.full_indep) : 0;
    const variance = indep + Math.pow(coverage, h.ld_share_power) * cross;
    const z = variance > 0 ? (raw - mean) / Math.sqrt(variance) : 0;
    const probability = liabilityRisk(z, m.liability_threshold, m.prevalence, m.r2);
    const res = { key: m.key, name: m.name, snps_used: used, clumps: clumpLen.length, coverage, z,
                  percentile: normalCdf(z), prevalence: m.prevalence, r2_liability: m.r2, r2_basis: m.r2_basis,
                  probability, relative_risk: probability / m.prevalence, low_power: m.r2 < h.low_power_r2,
                  warnings: [] };
    if (coverage < 0.2) {
      res.warnings.push("this file has few of the SNPs the score is built from; the estimate is less accurate " +
                        "than for a typical genotyping chip");
    }
    disorders.push(res);
  }
  return { valid, warnings, ancestry_fit: { value: fit, snps: fitSnps, good_below: h.fit_good, invalid_above: h.fit_poor },
           disorders };
}

/** Bytes of a genome file (possibly zipped or gzipped) -> full result. */
export async function scoreGenomeBytes(bundle, bytes, name) {
  const plain = await unwrap(bytes);
  if (isBundle(plain)) throw new PgcRiskError("That is the reference data file, not a genome file.");
  const text = new TextDecoder("utf-8").decode(plain);
  const { counts, genome } = readGenome(text, bundle, name);
  const result = scoreCounts(bundle, counts);
  result.genome = genome;
  return result;
}
