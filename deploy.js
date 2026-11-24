"use strict";
const fs = require("node:fs");
const path = require("node:path");

const ROOT = fs.realpathSync(__dirname);
const OUTPUT = path.join(ROOT, "_site");
const MARKER = ".astroinformatics-html";
const OWNER = "Astroinformatics HTML build\n";
const tools = JSON.parse(fs.readFileSync(path.join(ROOT, "html-tools.json"), "utf8"));
const escape = text => text.replace(/[&<>"']/g, character =>
  ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[character]);
const decodeText = text => text.replace(/&(#(?:x[\da-f]+|\d+)|amp|lt|gt|quot|apos);/gi, (_, entity) => {
  if (entity[0] === "#") return String.fromCodePoint(entity[1].toLowerCase() === "x" ?
    parseInt(entity.slice(2), 16) : Number(entity.slice(1)));
  return ({ amp: "&", lt: "<", gt: ">", quot: '"', apos: "'" })[entity.toLowerCase()];
});

function planAssets(catalogue = tools) {
  const assets = new Map(), paths = new Set();
  if (!Array.isArray(catalogue) || !catalogue.length) throw new Error("The HTML catalogue is empty");
  for (const tool of catalogue) {
    for (const field of ["path", "title", "description", "lecture", "section", "group"]) {
      if (typeof tool[field] !== "string" || !tool[field]) throw new Error(`Missing tool ${field}`);
    }
    if (!/^[A-Za-z][A-Za-z0-9_-]*$/.test(tool.path) || paths.has(tool.path)) {
      throw new Error(`Invalid or duplicate tool path: ${tool.path}`);
    }
    paths.add(tool.path);
    if (!Array.isArray(tool.static_files) || !tool.static_files.includes("index.html")) {
      throw new Error(`${tool.path}: missing static_files manifest or index.html`);
    }
    for (const asset of tool.static_files) {
      if (typeof asset !== "string" || /[:\\?#]/.test(asset) || asset.startsWith("/")) {
        throw new Error(`${tool.path}: invalid asset path ${asset}`);
      }
      const name = path.posix.normalize(`${tool.path}/${asset}`);
      if (name.startsWith("../") || !/\.(html|css|js)$/.test(name) || /shinylive|webr/i.test(name)) {
        throw new Error(`${tool.path}: not a local HTML runtime asset: ${asset}`);
      }
      const source = fs.realpathSync(path.join(ROOT, name));
      if (!source.startsWith(ROOT + path.sep)) throw new Error(`Asset escapes the repository: ${asset}`);
      if (!assets.has(name)) assets.set(name, fs.readFileSync(source));
    }
    const html = assets.get(`${tool.path}/index.html`).toString("utf8");
    const title = html.match(/<title>([\s\S]*?)<\/title>/i)?.[1];
    if (title === undefined || decodeText(title).trim() !== tool.title || !html.includes(tool.lecture)) {
      throw new Error(`${tool.path}: page title or lecture link differs from the catalogue`);
    }
    for (const match of html.matchAll(/<(?:script|link)\b[^>]*\b(?:src|href)=["']([^"']+)["'][^>]*>/gi)) {
      if (!tool.static_files.includes(match[1])) {
        throw new Error(`${tool.path}: ${match[1]} is not in the static_files manifest`);
      }
    }
  }
  return assets;
}

function gallery() {
  const groups = [...new Set(tools.map(tool => tool.group))];
  const sections = groups.map(group => {
    const cards = tools.filter(tool => tool.group === group).map((tool, i) => `
      <article class="card">
        <span class="eyebrow">"AstroInformatics Demonstration / ${String(i + 1).padStart(2, "0")}"</span>
        <h3><a href="${tool.path}/">${escape(tool.title)}</a></h3>
        <p>${escape(tool.description)}</p>
        <div class="card-links">
          <a class="explore" href="${tool.path}/">Explore demonstration &rarr;</a>
          <a href="${escape(tool.lecture)}">Lecture: ${escape(tool.section)}</a>
        </div>
      </article>`).join("\n");
    return `<section><h2>${escape(group)}</h2><div class="grid">${cards}</div></section>`;
  }).join("\n");
  return `<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Astroinformatics | See the statistics</title>
  <meta name="description" content="Interactive HTML demonstrations of astronomical selection, Bayesian inference, sampling, simulation and statistical foundations.">
  <style>
    *{box-sizing:border-box}body{margin:0;background:#f5f7fb;color:#25344b;font:16px/1.6 system-ui,-apple-system,sans-serif}
    header{background:linear-gradient(125deg,#182a46,#255e73);color:white;padding:64px max(24px,calc((100vw - 1180px)/2))}
    header p{max-width:760px;color:#dce8f1;font-size:1.15rem}h1{font-size:clamp(2rem,5vw,3.8rem);line-height:1.1;margin:12px 0}
    header a{color:#c6eee9}.eyebrow{font-size:.73rem;letter-spacing:.1em;text-transform:uppercase;font-weight:700;color:#087f8c}
    header .eyebrow{color:#b1e6df}main{max-width:1228px;margin:auto;padding:36px 24px 60px}
    section{margin:36px 0}h2{font-size:1.45rem;margin:0 0 18px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,310px),1fr));gap:20px}
    .card{background:white;border:1px solid #e0e5ef;border-radius:16px;padding:24px;display:flex;flex-direction:column;box-shadow:0 5px 16px #25344b05}
    .card h3{line-height:1.25;font-size:1.25rem;margin:12px 0}.card h3 a{color:#25344b;text-decoration:none}.card p{color:#536176;margin:0 0 24px}
    a{color:#365b99;text-underline-offset:3px}a:hover{text-decoration:underline}.card-links{display:grid;gap:12px;margin-top:auto;font-size:.84rem}
    .explore{font-size:.95rem;font-weight:700}footer{border-top:1px solid #dfe5ef;padding:24px;color:#536176;text-align:center}
    a:focus-visible{outline:3px solid #d56534;outline-offset:4px}.note{color:#536176;font-size:.9rem}
  </style>
</head>
<body>
  <header>
    <span class="eyebrow">Astroinformatics / Interactive lecture companions</span>
    <h1>See the statistics.</h1>
    <p>Explore principles behind astronomical inference, one visual experiment at a time. Change a few assumptions, reveal the hidden truth, and follow the evidence.</p>
    <a href="https://anguswright.github.io/LectureNotes/IntroductionToStatistics_materials.html">Open the lecture series</a>
  </header>
  <main>
    <p>All ${tools.length} demonstrations run entirely in your browser, with no installation or runtime downloads. Explore biases in inference, methods of sampling, and more. Every page links to relevant lecture notes, and explains any underlying assumptions.</p>
    ${sections}
    <p class="note">These are controlled teaching simulations, not analyses of real survey data.</a>.</p>
  </main>
  <footer>Astroinformatics &middot; Statistical principles for astronomers and physicists</footer>
</body>
</html>
`;
}

function buildSite(output = OUTPUT) {
  output = path.resolve(output);
  if (output === ROOT || ROOT.startsWith(output + path.sep) || output === path.parse(output).root) {
    throw new Error("The HTML output must not contain the source repository");
  }
  const assets = planAssets();
  assets.set("index.html", Buffer.from(gallery()));
  assets.set(".nojekyll", Buffer.alloc(0));
  assets.set(MARKER, Buffer.from(OWNER));
  if (fs.existsSync(output)) {
    if (!fs.lstatSync(output).isDirectory() ||
        (fs.readdirSync(output).length && (!fs.existsSync(path.join(output, MARKER)) ||
         fs.readFileSync(path.join(output, MARKER), "utf8") !== OWNER))) {
      throw new Error(`Refusing to replace an output directory not owned by the HTML build: ${output}`);
    }
  }
  const staging = fs.mkdtempSync(path.join(path.dirname(output), ".astro-html-"));
  try {
    for (const [name, content] of assets) {
      const target = path.join(staging, name);
      fs.mkdirSync(path.dirname(target), { recursive: true });
      fs.writeFileSync(target, content);
    }
    fs.rmSync(output, { recursive: true, force: true });
    fs.renameSync(staging, output);
  } finally {
    fs.rmSync(staging, { recursive: true, force: true });
  }
  return output;
}

if (require.main === module) {
  console.log(`Built ${tools.length} HTML tools in ${buildSite()}`);
}
module.exports = { tools, planAssets, buildSite };
