// check every docs page's math spans through the real katex engine.
// exits nonzero if any span fails to parse, so CI can gate on it.
// katex ships bundled inside the @quartz-community/latex plugin dist; we
// slice that module out of the bundle and run it in a vm.
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const dist = path.join(__dirname, "node_modules/@quartz-community/latex/dist/index.js");
const bundle = fs.readFileSync(dist, "utf8");

const start = bundle.indexOf("// node_modules/katex/dist/katex.mjs");
if (start < 0) throw new Error("katex module marker not found in latex plugin bundle");
const markers = [];
const re = /\/\/ node_modules\//g;
let mch;
while ((mch = re.exec(bundle)) !== null) markers.push(mch.index);
const idx = markers.indexOf(start);
const end = idx >= 0 && idx + 1 < markers.length ? markers[idx + 1] : bundle.length;
const mod = { exports: {} };
const ctx = { module: mod, exports: mod.exports, console, process };
vm.createContext(ctx);
vm.runInContext(bundle.slice(start, end) + "\nmodule.exports = katex;", ctx);
const katex = mod.exports;
if (!katex || typeof katex.renderToString !== "function") {
  throw new Error("katex extraction failed");
}

const docsDir = path.join(__dirname, "..", "docs");
const pages = [];
const walk = (dir) => {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, entry.name);
    if (entry.isDirectory()) walk(p);
    else if (entry.name.endsWith(".md")) pages.push(p);
  }
};
walk(docsDir);

let totalOk = 0, totalFail = 0;
for (const page of pages.sort()) {
  const src = fs.readFileSync(page, "utf8");
  const spans = [];
  for (const m of src.matchAll(/\$\$([^]+?)\$\$/g)) spans.push(m[1].trim());
  for (const m of src.matchAll(/(?<!\$)\$([^$\n]+?)\$(?!\$)/g)) spans.push(m[1].trim());
  if (!spans.length) continue;
  let ok = 0, fail = 0;
  for (const tex of spans) {
    try { katex.renderToString(tex, { throwOnError: true }); ok++; }
    catch (e) { fail++; console.log(`FAIL ${path.relative(docsDir, page)}: ${JSON.stringify(tex).slice(0, 60)} -> ${e.message.slice(0, 80)}`); }
  }
  totalOk += ok; totalFail += fail;
}
console.log(`docs math: ${totalOk} spans parsed OK, ${totalFail} failed across ${pages.length} pages`);
process.exit(totalFail > 0 ? 1 : 0);
