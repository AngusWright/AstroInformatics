(function () {
  "use strict";
  const U = AstroUI, M = MarginalisationModel;
  const purple = "#5059A6", teal = "#087F8C", orange = "#D56534";
  const pct = value => `${U.number(100 * value)}%`;
  function render(result, p) {
    U.$("broad-mass").textContent = pct(p.weight);
    U.$("interval-mass").textContent = pct(result.probability);
    U.$("mode").textContent = result.mode.toFixed(2);
    U.$("window-mass").textContent = pct(result.mass);
    const joint = U.plot(U.$("joint"), [-5, 7], [-6, 8], "Parameter x", "Nuisance parameter y");
    const peak = Math.max(...result.joint.flat());
    joint.clip();
    U.heatmap(joint, result.x, result.y, result.joint, value => {
      const t = Math.sqrt(value / peak), white = [255, 255, 255], dark = [80, 89, 166];
      return `rgb(${white.map((v, i) => Math.round(v + t * (dark[i] - v))).join(",")})`;
    });
    joint.ctx.fillStyle = "#087F8C20";
    joint.ctx.fillRect(joint.xp(result.lower), joint.top, joint.xp(result.upper) - joint.xp(result.lower), joint.bottom - joint.top);
    U.line(joint, [[result.lower, -6], [result.lower, 8]], teal);
    U.line(joint, [[result.upper, -6], [result.upper, 8]], teal);
    joint.ctx.restore();
    const ymax = Math.max(...result.marginal, ...(p.show_profile ? result.profile : [0])) * 1.12;
    const marginal = U.plot(U.$("marginal"), [-5, 7], [0, ymax], "Parameter x", "Density in x");
    marginal.clip();
    const band = U.linspace(result.lower, result.upper, 150).map(x => [x, M.mixtureMarginal(x, p.weight, p.width)]);
    U.area(marginal, band, "#087F8C66");
    U.line(marginal, result.x.map((x, i) => [x, result.marginal[i]]), purple, 3);
    if (p.show_profile) U.line(marginal, result.x.map((x, i) => [x, result.profile[i]]), orange, 2, [6, 4]);
    marginal.ctx.restore();
    const peakFactor = 2 * Math.PI * Math.sqrt(1 - 0.45 ** 2);
    U.table("components", [
      ["Narrow", pct(1 - p.weight), U.number((1 - p.weight) / (peakFactor * 0.23 * 0.3)), "0.23, 0.30"],
      ["Broad", pct(p.weight), U.number(p.weight / (peakFactor * p.width ** 2)), `${U.number(p.width)}, ${U.number(p.width)}`]
    ]);
    U.$("interval-summary").textContent =
      `P(${result.lower.toFixed(2)} ≤ x ≤ ${result.upper.toFixed(2)}) ≈ ${pct(result.probability)}, integrating over all y. ` +
      `The image contains approximately ${pct(result.mass)} of the total joint probability; omitted tails are not renormalised.`;
  }
  window.lesson = AstroLesson.mount({
    compute: M.compute, render,
    presets: {
      tall: { weight: 0.7, width: 1.5, centre: 1.2, show_profile: false },
      narrow: { weight: 0.15, width: 1, centre: -1, show_profile: false },
      profile: { weight: 0.7, width: 1.5, centre: 1.2, show_profile: true }
    }
  });
})();
