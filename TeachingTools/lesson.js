(function (root) {
  "use strict";
  const A = root.Astro, U = root.AstroUI;

  function mount({ compute, render, presets = {}, steps = 0 }) {
    const form = U.$("controls"), error = U.$("error");
    const fields = [...form.elements].filter(field => field.id && field.matches("input, select"));
    let generation = 0, result = null, parameters = null, frame = steps, playing = false, timer = null;
    const animation = ["play", "step", "rewind", "finish", "frame"].map(U.$).filter(Boolean);
    if (U.$("frame")) { U.$("frame").min = "0"; U.$("frame").max = String(steps); U.$("frame").step = "1"; }

    function pause() {
      playing = false;
      if (timer !== null) clearTimeout(timer);
      timer = null;
      if (U.$("play")) { U.$("play").textContent = "Play animation"; U.$("play").setAttribute("aria-pressed", "false"); }
    }
    function read() {
      const values = {};
      for (const field of fields) {
        if (field.type === "checkbox") { values[field.id] = field.checked; continue; }
        if (field.type === "number" || field.type === "range") {
          const value = field.valueAsNumber, name = field.labels[0]?.textContent.trim() || field.id;
          A.check(Number.isFinite(value), `${name} must be a finite number.`);
          A.check((field.min === "" || value >= Number(field.min)) && (field.max === "" || value <= Number(field.max)),
            `${name} must be between ${field.min || "-infinity"} and ${field.max || "infinity"}.`);
          A.check(!field.validity.stepMismatch, `${name} must use increments of ${field.step}.`);
          values[field.id] = value;
        } else {
          A.check([...field.options].some(option => option.value === field.value && !option.disabled),
            `Choose a valid ${field.id}.`);
          values[field.id] = field.value;
        }
      }
      if (values.seed !== undefined) {
        A.check(Number.isInteger(values.seed), "Seed must be an integer.");
        values.seed = (values.seed + generation) % 2147483648;
      }
      return values;
    }
    function attempt(action) {
      try { action(); }
      catch (failure) {
        if (!(failure instanceof A.ValidationError)) throw failure;
        pause(); result = null;
        error.hidden = false; error.textContent = failure.message;
        document.querySelectorAll("canvas").forEach(canvas => U.blank(canvas, "Adjust the settings to continue."));
        document.querySelectorAll("[data-result]").forEach(node => { node.textContent = "\u2014"; });
        document.querySelectorAll("tbody").forEach(node => node.replaceChildren());
        animation.forEach(control => { control.disabled = true; });
      }
    }
    function draw() {
      if (result === null) return;
      attempt(() => {
        if (U.$("frame")) U.$("frame").value = String(frame);
        if (U.$("frame-label")) U.$("frame-label").textContent = `${frame} / ${steps}`;
        render(result, parameters, frame);
      });
    }
    function refresh() {
      pause(); U.outputs();
      attempt(() => {
        parameters = read(); result = compute(parameters); frame = steps;
        A.check(result !== null && result !== undefined, "The model did not produce a result.");
        error.hidden = true; error.textContent = "";
        animation.forEach(control => { control.disabled = false; });
        draw();
      });
    }
    function set(values) {
      for (const [id, value] of Object.entries(values)) {
        const field = U.$(id);
        if (!field || !fields.includes(field)) throw new Error(`Unknown lesson control: ${id}`);
        if (field.type === "checkbox") field.checked = value;
        else field.value = String(value);
      }
      if (Object.hasOwn(values, "seed")) generation = 0;
      refresh();
    }
    function tick() {
      if (!playing) return;
      frame = Math.min(steps, frame + 1); draw();
      if (frame >= steps || result === null) pause();
      else timer = setTimeout(tick, Number(U.$("delay")?.value || 200));
    }
    form.addEventListener("submit", event => event.preventDefault());
    form.addEventListener("input", U.outputs);
    form.addEventListener("change", event => {
      if (event.target.id === "delay") return;
      if (event.target.id === "seed") generation = 0;
      refresh();
    });
    U.$("redraw")?.addEventListener("click", () => { generation++; refresh(); });
    U.$("reset")?.addEventListener("click", () => { form.reset(); generation = 0; refresh(); });
    document.querySelectorAll("[data-preset]").forEach(button => button.addEventListener("click", () => {
      const values = presets[button.dataset.preset];
      if (!values) throw new Error(`Unknown lesson preset: ${button.dataset.preset}`);
      set(values);
    }));
    U.$("play")?.addEventListener("click", () => {
      if (playing) { pause(); return; }
      if (frame >= steps) frame = 0;
      playing = true; U.$("play").textContent = "Pause"; U.$("play").setAttribute("aria-pressed", "true");
      tick();
    });
    for (const [id, move] of [["step", () => Math.min(steps, frame + 1)], ["rewind", () => 0], ["finish", () => steps]]) {
      U.$(id)?.addEventListener("click", () => { pause(); frame = move(); draw(); });
    }
    U.$("frame")?.addEventListener("input", event => { pause(); frame = Number(event.target.value); draw(); });

    const tabs = [...document.querySelectorAll('[role="tab"]')];
    function selectTab(selected, focus = false) {
      for (const tab of tabs) {
        const active = tab === selected;
        tab.setAttribute("aria-selected", String(active)); tab.tabIndex = active ? 0 : -1;
        U.$(tab.getAttribute("aria-controls")).hidden = !active;
      }
      if (focus) selected.focus();
      draw();
    }
    tabs.forEach((tab, i) => {
      tab.addEventListener("click", () => selectTab(tab));
      tab.addEventListener("keydown", event => {
        const index = event.key === "ArrowRight" ? (i + 1) % tabs.length :
          event.key === "ArrowLeft" ? (i + tabs.length - 1) % tabs.length :
          event.key === "Home" ? 0 : event.key === "End" ? tabs.length - 1 : -1;
        if (index >= 0) { event.preventDefault(); selectTab(tabs[index], true); }
      });
    });
    if (U.$("settings") && matchMedia("(max-width: 800px)").matches) U.$("settings").open = false;
    window.addEventListener("resize", draw);
    window.addEventListener("pagehide", pause);
    refresh();
    document.documentElement.dataset.ready = "true";
    return { refresh, set, redraw: draw, get result() { return result; }, get parameters() { return parameters; },
      get frame() { return frame; }, get generation() { return generation; } };
  }

  root.AstroLesson = { mount };
})(globalThis);
