"use strict";
const assert = require("node:assert/strict");
const { spawn } = require("node:child_process");
const { mkdtemp, readFile, rm, mkdir, writeFile } = require("node:fs/promises");
const http = require("node:http");
const os = require("node:os");
const path = require("node:path");
const { pathToFileURL } = require("node:url");
const { setTimeout: sleep } = require("node:timers/promises");

async function runBrowser(page, suite) {
  const root = path.resolve(__dirname, "..", process.env.HTML_SITE_ROOT || ".");
  const directory = path.join(root, page);
  assert.ok(page === "" || /^[A-Za-z][A-Za-z0-9_-]*$/.test(page), "Invalid tool directory");
  const route = page ? `/${page}/` : "/", name = page || "Gallery";
  const html = await readFile(path.join(directory, "index.html"), "utf8");
  const files = new Map([[route, path.join(directory, "index.html")],
    [`${route}index.html`, path.join(directory, "index.html")]]);
  for (const match of html.matchAll(/(?:src|href)="([^"]+\.(?:js|css))"/g)) {
    const url = new URL(match[1], `http://localhost${route}`);
    assert.equal(url.origin, "http://localhost", "Runtime assets must be local");
    const file = path.resolve(root, "." + decodeURIComponent(url.pathname));
    assert.ok(file.startsWith(root + path.sep), "Runtime assets must stay within the repository");
    files.set(url.pathname, file);
  }
  const profile = await mkdtemp(path.join(os.tmpdir(), `${name.toLowerCase()}-chrome-`));
  const types = { ".html": "text/html", ".css": "text/css", ".js": "text/javascript" };
  const server = http.createServer(async (req, res) => {
    const pathname = new URL(req.url, "http://localhost").pathname;
    if (pathname === "/favicon.ico") { res.writeHead(204); res.end(); return; }
    if (!files.has(pathname)) { res.writeHead(404); res.end("Not found"); return; }
    try {
      const file = files.get(pathname), body = await readFile(file);
      res.writeHead(200, { "Content-Type": types[path.extname(file)] }); res.end(body);
    } catch (error) {
      res.writeHead(500); res.end(error.message);
    }
  });
  let chrome, socket, sequence = 0, stderr = "";
  const pending = new Map(), exceptions = [], requests = [];
  try {
    await new Promise((resolve, reject) => {
      server.once("error", reject); server.listen(0, "127.0.0.1", resolve);
    });
    const origin = `http://127.0.0.1:${server.address().port}`, base = `${origin}${route}`;
    const binary = process.env.CHROME_BIN ||
      (process.platform === "darwin" ? "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" : "google-chrome");
    chrome = spawn(binary, ["--headless", "--remote-debugging-port=0", `--user-data-dir=${profile}`,
      "--no-first-run", "--no-default-browser-check", "--disable-background-networking", "about:blank"],
    { stdio: ["ignore", "ignore", "pipe"] });
    const endpoint = await new Promise((resolve, reject) => {
      const timeout = setTimeout(() => reject(new Error(`Chrome startup timed out:\n${stderr}`)), 20000);
      const finish = (error, value) => { clearTimeout(timeout); error ? reject(error) : resolve(value); };
      chrome.once("error", error => finish(error));
      chrome.once("exit", code => finish(new Error(`Chrome exited (${code}):\n${stderr}`)));
      chrome.stderr.on("data", data => {
        stderr = (stderr + data.toString()).slice(-12000);
        const match = stderr.match(/DevTools listening on (ws:\/\/\S+)/);
        if (match) finish(null, match[1]);
      });
    });
    const targets = await (await fetch(`http://${new URL(endpoint).host}/json/list`)).json();
    socket = new WebSocket(targets.find(target => target.type === "page").webSocketDebuggerUrl);
    await new Promise((resolve, reject) => {
      socket.addEventListener("open", resolve, { once: true });
      socket.addEventListener("error", reject, { once: true });
    });
    socket.addEventListener("message", event => {
      const message = JSON.parse(event.data);
      if (message.id) {
        const request = pending.get(message.id);
        if (!request) return;
        clearTimeout(request.timeout); pending.delete(message.id);
        if (message.error) request.reject(new Error(JSON.stringify(message.error)));
        else request.resolve(message.result);
      } else if (message.method === "Runtime.exceptionThrown") {
        exceptions.push(message.params.exceptionDetails.exception?.description || message.params.exceptionDetails.text);
      } else if (message.method === "Network.requestWillBeSent") requests.push(message.params.request.url);
    });
    function call(method, params = {}) {
      return new Promise((resolve, reject) => {
        const id = ++sequence;
        const timeout = setTimeout(() => { pending.delete(id); reject(new Error(`Browser command timed out: ${method}`)); }, 15000);
        pending.set(id, { resolve, reject, timeout }); socket.send(JSON.stringify({ id, method, params }));
      });
    }
    async function evaluate(expression) {
      const result = await call("Runtime.evaluate", { expression, returnByValue: true, awaitPromise: true });
      if (result.exceptionDetails) throw new Error(result.exceptionDetails.exception?.description || result.exceptionDetails.text);
      return result.result.value;
    }
    async function waitFor(expression) {
      for (let i = 0; i < 200; i++) { if (await evaluate(expression)) return; await sleep(25); }
      throw new Error(`Page condition timed out: ${expression}`);
    }
    const click = id => evaluate(`document.getElementById(${JSON.stringify(id)}).click()`);
    const text = id => evaluate(`document.getElementById(${JSON.stringify(id)}).textContent`);
    const value = id => evaluate(`document.getElementById(${JSON.stringify(id)}).value`);
    const change = (id, value) => evaluate(`(() => {
      const input = document.getElementById(${JSON.stringify(id)});
      ${typeof value === "boolean" ? `input.checked = ${value};` : `input.value = ${JSON.stringify(String(value))};`}
      input.dispatchEvent(new Event("change", {bubbles:true}));
    })()`);
    async function healthy() {
      if (page) assert.equal(await evaluate('document.getElementById("error").hidden'), true, await text("error"));
      assert.deepEqual(exceptions, []);
    }
    async function navigate(url = base) {
      await evaluate("delete document.documentElement?.dataset.ready");
      const result = await call("Page.navigate", { url });
      if (result.errorText) throw new Error(result.errorText);
      await waitFor(page ? 'document.documentElement?.dataset.ready === "true" && document.readyState === "complete"' :
        `location.href === ${JSON.stringify(url)} && document.readyState === "complete"`);
      await healthy();
    }
    async function screenshot(name) {
      const destination = process.env.HTML_SCREENSHOT_DIR || process.env.GIBBS_SCREENSHOT_DIR;
      if (!destination) return;
      await mkdir(destination, { recursive: true });
      const image = await call("Page.captureScreenshot", { format: "png", captureBeyondViewport: true });
      await writeFile(path.join(destination, name), Buffer.from(image.data, "base64"));
    }
    async function hasColour(id, colour) {
      return evaluate(`(() => {
        const canvas = document.getElementById(${JSON.stringify(id)});
        const pixels = canvas.getContext("2d").getImageData(0,0,canvas.width,canvas.height).data;
        const colour = ${JSON.stringify(colour)};
        for(let i=0;i<pixels.length;i+=4)
          if(colour.every((channel,j)=>Math.abs(channel-pixels[i+j])<4)) return true;
        return false;
      })()`);
    }
    await Promise.all([call("Page.enable"), call("Runtime.enable"), call("Network.enable")]);
    await call("Emulation.setDeviceMetricsOverride", { width: 1440, height: 1100, deviceScaleFactor: 1, mobile: false });
    await navigate();
    await suite({ evaluate, click, change, text, value, waitFor, healthy, navigate,
      openFile: () => navigate(pathToFileURL(path.join(directory, "index.html")).href),
      screenshot, hasColour, call, base });
    await healthy();
    const allowed = new Set([origin + "/favicon.ico", ...[...files].flatMap(([pathname, file]) =>
      [origin + pathname, pathToFileURL(file).href])]);
    assert.ok(requests.every(url => allowed.has(url)), `Unexpected runtime request: ${requests.filter(url => !allowed.has(url)).join(", ")}`);
    console.log(`${name}: browser interactions, plots, and local-only runtime passed.`);
  } finally {
    for (const request of pending.values()) { clearTimeout(request.timeout); request.reject(new Error("Browser closed")); }
    if (socket) socket.close();
    if (chrome && chrome.pid !== undefined && chrome.exitCode === null && chrome.signalCode === null) {
      await new Promise(resolve => {
        const timeout = setTimeout(() => chrome.kill("SIGKILL"), 5000);
        chrome.once("exit", () => { clearTimeout(timeout); resolve(); }); chrome.kill("SIGTERM");
      });
    }
    await new Promise(resolve => server.close(resolve));
    await rm(profile, { recursive: true, force: true, maxRetries: 5, retryDelay: 100 });
  }
}
module.exports = { runBrowser };
