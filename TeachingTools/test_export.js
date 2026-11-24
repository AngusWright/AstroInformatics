"use strict";
const assert = require("node:assert/strict");
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const test = require("node:test");
const { tools, planAssets, buildSite } = require("../deploy.js");

const ROOT = path.resolve(__dirname, ".."), OUTPUT = path.join(ROOT, "_site");
function filesIn(directory) {
  return fs.readdirSync(directory, { withFileTypes: true }).flatMap(entry =>
    entry.isDirectory() ? filesIn(path.join(directory, entry.name)).map(name => `${entry.name}/${name}`) : entry.name
  ).sort();
}

test("the published artifact contains exactly the declared native assets, copied without changes", () => {
  const expected = new Set(["index.html", ".nojekyll", ".astroinformatics-html"]);
  for (const tool of tools) {
    for (const asset of tool.static_files) {
      const name = path.posix.normalize(`${tool.path}/${asset}`);
      expected.add(name);
      assert.deepEqual(fs.readFileSync(path.join(OUTPUT, name)), fs.readFileSync(path.join(ROOT, name)), name);
    }
  }
  const actual = filesIn(OUTPUT);
  assert.deepEqual(actual, [...expected].sort());
  assert.ok(actual.every(name => !/shinylive|webr|app\.json|(^|\/)edit\/|\.R$/i.test(name)));
});

test("gallery order, lectures, overview and all local links work beneath the GitHub project path", () => {
  const base = new URL("https://example.test/Astroinformatics/");
  const gallery = fs.readFileSync(path.join(OUTPUT, "index.html"), "utf8");
  const overview = fs.readFileSync(path.join(ROOT, "AstroInformatics.qmd"), "utf8");
  let previous = -1;
  for (const tool of tools) {
    const position = gallery.indexOf(`href="${tool.path}/"`);
    assert.ok(position > previous, `${tool.path}: missing or out-of-order gallery link`);
    previous = position;
    assert.ok(gallery.includes(tool.lecture), `${tool.path}: missing lecture`);
    assert.ok(overview.includes(`](${tool.path}/)`), `${tool.path}: missing overview link`);
  }
  for (const page of ["index.html", ...tools.map(tool => `${tool.path}/index.html`)]) {
    const html = fs.readFileSync(path.join(OUTPUT, page), "utf8");
    const document = new URL(page, base);
    for (const match of html.matchAll(/\b(?:href|src)="([^"]+)"/g)) {
      const target = new URL(match[1], document);
      if (target.origin !== base.origin) continue;
      assert.ok(target.pathname.startsWith(base.pathname), `${page}: link escapes the project: ${match[1]}`);
      let name = decodeURIComponent(target.pathname.slice(base.pathname.length));
      if (!name || name.endsWith("/")) name += "index.html";
      assert.ok(fs.statSync(path.join(OUTPUT, name)).isFile(), `${page}: missing ${match[1]}`);
    }
  }
});

test("the publisher rejects omitted dependencies, non-native payloads and escaping paths", () => {
  const missing = structuredClone(tools);
  const conjugate = missing.find(tool => tool.path === "Conjugate");
  conjugate.static_files = conjugate.static_files.filter(asset => asset !== "model.js");
  assert.throws(() => planAssets(missing), /Conjugate: model\.js is not in the static_files manifest/);
  for (const asset of ["../../outside.js", "../docs/shinylive/load-shinylive.js", "app.json", "https://example.test/app.js"]) {
    const invalid = structuredClone(tools);
    invalid[0].static_files.push(asset);
    assert.throws(() => planAssets(invalid), /invalid asset path|not a local HTML runtime asset/);
  }
});

test("rebuilds discard stale payloads but do not replace unrelated directories", () => {
  const temporary = fs.mkdtempSync(path.join(os.tmpdir(), "astro-html-export-"));
  try {
    const output = path.join(temporary, "site");
    buildSite(output);
    fs.mkdirSync(path.join(output, "shinylive"));
    fs.writeFileSync(path.join(output, "shinylive", "old.js"), "stale");
    fs.writeFileSync(path.join(output, "Conjugate", "app.json"), "stale");
    buildSite(output);
    assert.deepEqual(filesIn(output), filesIn(OUTPUT));
    assert.equal(fs.existsSync(path.join(output, "shinylive")), false);
    assert.equal(fs.existsSync(path.join(output, "Conjugate", "app.json")), false);

    const unrelated = path.join(temporary, "unrelated");
    fs.mkdirSync(unrelated);
    fs.writeFileSync(path.join(unrelated, "keep.txt"), "preserve");
    assert.throws(() => buildSite(unrelated), /not owned by the HTML build/);
    assert.equal(fs.readFileSync(path.join(unrelated, "keep.txt"), "utf8"), "preserve");
    assert.throws(() => buildSite(ROOT), /must not contain the source repository/);
  } finally {
    fs.rmSync(temporary, { recursive: true, force: true });
  }
});
