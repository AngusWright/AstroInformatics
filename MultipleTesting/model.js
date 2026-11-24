(function (root) {
  "use strict";
  const A = typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro;

  function searchFalseProbability(threshold, trials, twoSided = false) {
    A.check(Number.isFinite(threshold) && threshold >= 0, "Threshold must be non-negative and finite.");
    A.check(Number.isSafeInteger(trials) && trials >= 0, "Trial count must be a non-negative integer.");
    A.check(typeof twoSided === "boolean", "Choose a one- or two-sided search.");
    const single = A.normalCdf(-threshold) * (twoSided ? 2 : 1);
    return { single, global: trials === 0 ? 0 : -Math.expm1(trials * Math.log1p(-single)), expected: trials * single };
  }

  function compute(p) {
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    A.check(Number.isInteger(p.side) && p.side >= 5 && p.side <= 60, "Side length must be an integer from 5 to 60.");
    A.check(Number.isFinite(p.threshold) && p.threshold >= 1.5 && p.threshold <= 5, "Threshold must be between 1.5 and 5 sigma.");
    A.check(["one", "two"].includes(p.tails), "Choose positive peaks or either sign.");
    A.check(Number.isFinite(p.signal_fraction) && p.signal_fraction >= 0 && p.signal_fraction <= 0.1, "Signal fraction must be between 0 and 0.1.");
    A.check(Number.isFinite(p.amplitude) && p.amplitude >= 0 && p.amplitude <= 8, "Signal amplitude must be between 0 and 8 sigma.");
    const rng = A.createRng(p.seed);
    // Fixed-size streams keep overlapping cells unchanged when the search area changes.
    const noise = Array.from({ length: 3600 }, () => rng.normal());
    const membership = Array.from({ length: 3600 }, () => rng.uniform());
    const cells = noise.slice(0, p.side ** 2).map((value, i) => {
      const truth = membership[i] < p.signal_fraction;
      const z = value + (truth ? p.amplitude : 0);
      return { noise: value, truth, z, hit: (p.tails === "two" ? Math.abs(z) : z) >= p.threshold };
    });
    return { cells, side: p.side, threshold: p.threshold, twoSided: p.tails === "two" };
  }

  function reveal(result, frame = 100) {
    A.check(Number.isFinite(frame) && frame >= 0 && frame <= 100, "Search progress must be between 0 and 100.");
    const tested = Math.floor(result.cells.length * frame / 100);
    let trueDiscoveries = 0, falseDiscoveries = 0, nullTrials = 0;
    for (let i = 0; i < tested; i++) {
      const cell = result.cells[i];
      if (!cell.truth) nullTrials++;
      if (cell.hit) cell.truth ? trueDiscoveries++ : falseDiscoveries++;
    }
    const discoveries = trueDiscoveries + falseDiscoveries;
    return { tested, nullTrials, trueDiscoveries, falseDiscoveries, discoveries,
      falseDiscoveryFraction: discoveries ? falseDiscoveries / discoveries : null,
      probability: searchFalseProbability(result.threshold, nullTrials, result.twoSided) };
  }

  const api = { searchFalseProbability, compute, reveal };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.MultipleTestingModel = api;
})(globalThis);
