"use strict";
const assert = require("node:assert/strict");
const { runLessonBrowser } = require("../TeachingTools/lesson_browser_test.js");

runLessonBrowser("RobustStatistics", async b => {
  assert.equal(await b.evaluate("document.title"), "Mean Sensitivity: One Bad Star");
  assert.equal(await b.evaluate('document.querySelector(".lecture-link").href'),
    "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section1b_html.html");
  assert.ok(await b.hasColour("sample-plot", [213, 101, 52]));
  assert.ok(await b.hasColour("sample-plot", [80, 89, 166]));
  assert.ok(await b.hasColour("influence-plot", [0, 114, 178]));
  await b.evaluate("window.originalRobustData = lesson.result.data");
  const clean = await b.evaluate("lesson.result.clean"), jitter = await b.evaluate("lesson.result.jitter");
  const original = await b.evaluate("lesson.result.summaries");
  const image = await b.evaluate('document.getElementById("sample-plot").toDataURL()');
  await b.change("position", 12);
  const moved = await b.evaluate("lesson.result.summaries");
  assert.ok(Math.abs(moved.mean - original.mean - 4 / 40) < 1e-12);
  assert.equal(moved.median, original.median); assert.equal(moved.nmad, original.nmad);
  assert.deepEqual(await b.evaluate("lesson.result.clean"), clean);
  assert.deepEqual(await b.evaluate("lesson.result.jitter"), jitter);
  assert.equal(await b.evaluate("lesson.result.data === originalRobustData"), true);
  assert.notEqual(await b.evaluate('document.getElementById("sample-plot").toDataURL()'), image);
  await b.change("trim", 0);
  assert.equal(await b.text("trimmed"), await b.text("mean"));
  await b.change("trim", 0.4); await b.healthy();
  await b.change("contaminants", 0);
  assert.deepEqual(await b.evaluate("lesson.result.observed"), clean);
  await b.change("position", -12);
  assert.deepEqual(await b.evaluate("lesson.result.observed"), clean, "No effect when no values are replaced");
  for (const n of [20, 120]) {
    await b.change("n", n); await b.healthy();
    assert.equal(await b.evaluate("lesson.result.observed.length"), n);
  }
  await b.change("n", 20); await b.change("contaminants", 10); await b.change("position", 12);
  assert.ok(Number(await b.text("median")) > 4, "Half-sample contamination can move even the median");
  await b.click("reset"); await b.click("redraw");
  await b.evaluate("window.dragData = lesson.result.data; window.beforeDrag = lesson.result.observed.slice()");
  const handle = await b.evaluate(`(() => {
    const canvas = document.getElementById("sample-plot"), rect = canvas.getBoundingClientRect();
    const ratio = canvas.width / canvas.clientWidth;
    const image = canvas.getContext("2d").getImageData(0, 0, canvas.width, Math.ceil(85 * ratio));
    let sx = 0, sy = 0, count = 0;
    for (let y = 0; y < image.height; y++) for (let x = 0; x < image.width; x++) {
      const i = 4 * (y * image.width + x);
      if (image.data[i] === 213 && image.data[i+1] === 101 && image.data[i+2] === 52) {
        sx += x; sy += y; count++;
      }
    }
    return count ? { x: rect.left + sx/count/ratio, y: rect.top + sy/count/ratio, left: rect.left, right: rect.right } : null;
  })()`);
  assert.ok(handle, "The draggable contaminant handle must actually be painted");
  await b.call("Input.dispatchMouseEvent", { type: "mouseMoved", x: handle.x, y: handle.y });
  await b.call("Input.dispatchMouseEvent", { type: "mousePressed", x: handle.x, y: handle.y, button: "left", buttons: 1, clickCount: 1 });
  await b.call("Input.dispatchMouseEvent", { type: "mouseMoved", x: handle.left + 90, y: handle.y, button: "left", buttons: 1 });
  assert.ok(Number(await b.value("position")) < -8, "Dragging must update the real control");
  assert.ok(await b.evaluate("lesson.result.observed[0] < beforeDrag[0] - 15"), "Dragging moves the replacement in the statistical model");
  assert.equal(await b.evaluate("lesson.result.data === dragData && lesson.parameters.seed === 43"), true,
    "Dragging after New realisation must not reset its generation");
  assert.deepEqual(await b.evaluate("lesson.result.observed.slice(1)"), await b.evaluate("beforeDrag.slice(1)"));
  await b.call("Input.dispatchMouseEvent", { type: "mouseMoved", x: handle.left + 1, y: handle.y, button: "left", buttons: 1 });
  assert.equal(await b.value("position"), "-12");
  await b.call("Input.dispatchMouseEvent", { type: "mouseMoved", x: handle.right + 20, y: handle.y, button: "left", buttons: 1 });
  assert.equal(await b.value("position"), "12");
  await b.call("Input.dispatchMouseEvent", { type: "mouseReleased", x: handle.right + 20, y: handle.y, button: "left", buttons: 0, clickCount: 1 });
  await b.healthy();
  await b.change("seed", 0); await b.healthy();
  await b.change("seed", 2147483647); await b.healthy();
  await b.change("seed", 1.5);
  assert.equal(await b.evaluate('document.getElementById("error").hidden'), false);
  assert.equal(await b.text("mean"), "—");
  await b.click("reset"); await b.healthy();
}).catch(error => { console.error(error); process.exitCode = 1; });
