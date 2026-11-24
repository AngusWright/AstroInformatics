(function () {
  'use strict';
  const A = window.Astro, U = window.AstroUI, M = window.MCMCModel, $ = U.$;
  const colours = ['#0072b2', '#b45309', '#008260', '#8d4eaa', '#c3415b', '#607d23', '#795548', '#487b91'];
  const tabNames = ['explore', 'diagnostics', 'dependence', 'explain'];
  let result = null, reference = null, cursor = 0, active = 'explore', generation = 0;
  let playing = false, timer = null, animation = null, animationFrame = null;
  let diagnosticCache = null;

  function numeric(id) {
    const raw = $(id).value.trim();
    A.check(raw !== '', `${$(id).labels[0]?.textContent || id} is required.`);
    const value = Number(raw);
    A.check(Number.isFinite(value), `${id} must be a finite number.`);
    return value;
  }
  function config() {
    const family = $('family').value;
    return {
      family, observations: numeric('observations'),
      theta: [numeric('theta1'), numeric('theta2')], dataSd: family === 'gaussian' ? numeric('data-sd') : 1.2,
      rho: family === 'gaussian' ? numeric('rho') : 0.8,
      bounds: [[numeric('lower1'), numeric('upper1')], [numeric('lower2'), numeric('upper2')]],
      prior: $('prior').value, proposal: $('proposal').value, proposalSd: numeric('proposal-sd'),
      globalWeight: numeric('global-weight'), iterations: numeric('iterations'),
      chains: numeric('chains'), seed: numeric('seed'), starts: null
    };
  }
  function display() {
    const burn = numeric('burn'), thin = numeric('thin'), focus = numeric('focus-chain');
    A.check(Number.isInteger(burn) && burn >= 0 && burn <= result.config.iterations,
      'Burn-in must be an integer from zero to the iteration count.');
    A.check(Number.isInteger(thin) && thin >= 1 && thin <= 10000, 'Thinning must be an integer from 1 to 10,000.');
    A.check(Number.isInteger(focus) && focus >= 0 && focus < result.chains.length, 'Choose an existing chain.');
    return { burn, thin, focus, trail: numeric('trail'), delay: numeric('delay'),
      lines: $('lines').checked, prior: $('show-prior').checked };
  }
  function pause() {
    playing = false;
    clearTimeout(timer);
    timer = null;
    $('play').textContent = 'Play';
    $('playback-state').textContent = 'Paused';
  }
  function stopAnimation() {
    animation = null;
    if (animationFrame !== null) cancelAnimationFrame(animationFrame);
    animationFrame = null;
  }
  function invalidate(message, discardModel) {
    pause();
    stopAnimation();
    if (discardModel) { result = null; reference = null; }
    diagnosticCache = null;
    $('error').textContent = message;
    $('error').hidden = false;
    document.querySelectorAll('canvas').forEach(canvas => U.blank(canvas, 'Correct the settings to continue.'));
    ['step-table', 'sample-summary', 'chain-summary'].forEach(id => U.table(id, []));
    ['step-description', 'decision-detail', 'sample-count', 'data-summary', 'progress-text'].forEach(id => { $(id).textContent = ''; });
    ['play', 'step', 'rewind', 'finish', 'diagnostic-finish', 'seek'].forEach(id => { $(id).disabled = true; });
  }
  function guarded(fn, discardModel = false) {
    try {
      fn();
      $('error').hidden = true;
      $('error').textContent = '';
      ['play', 'step', 'rewind', 'finish', 'diagnostic-finish', 'seek'].forEach(id => { $(id).disabled = !result; });
      return true;
    } catch (error) {
      if (!(error instanceof A.ValidationError)) throw error;
      invalidate(error.message, discardModel);
      return false;
    }
  }
  function familyControls(reset) {
    const family = $('family').value, defaults = M.defaults(family);
    if (reset) {
      $('theta1').value = defaults.theta[0]; $('theta2').value = defaults.theta[1];
      ['lower1', 'upper1', 'lower2', 'upper2'].forEach((id, i) => { $(id).value = defaults.bounds[Math.floor(i / 2)][i % 2]; });
    }
    const labels = family === 'gaussian' ? ['mean x', 'mean y'] : family === 'binomial'
      ? ['probability p', 'draw count N'] : family === 'beta' ? ['shape α', 'shape β'] : ['shape k', 'rate r'];
    for (let d = 0; d < 2; d++) {
      $(`theta${d + 1}-label`).textContent = `Generating ${labels[d]}`;
      $(`lower${d + 1}-label`).textContent = `${labels[d]}: lower bound`;
    }
    $('theta2').step = family === 'binomial' ? '1' : 'any';
    $('lower2').step = $('upper2').step = family === 'binomial' ? '1' : 'any';
    $('gaussian-settings').hidden = family !== 'gaussian';
    $('mixture-settings').hidden = $('proposal').value !== 'mixture';
    $('proposal-sd').disabled = $('proposal').value === 'independence';
    const mixed = family === 'binomial' ? ' N moves by −1, 0 or +1; this SD applies only to p.' : ' This is the SD in each coordinate, not the variance.';
    $('proposal-hint').textContent = $('proposal').value === 'randomwalk' ? `Symmetric random walk.${mixed}`
      : $('proposal').value === 'independence'
        ? 'Global normals use SD = prior width / 3; N uses a uniform integer proposal. The Hastings correction is required.'
        : `The proposal-density ratio includes both mixture components.${mixed}`;
    U.outputs();
  }
  function regenerate() {
    pause(); stopAnimation(); diagnosticCache = null;
    const built = guarded(() => {
      familyControls(false);
      const next = M.run(config()), nextReference = M.reference(next.target);
      result = next; reference = nextReference; cursor = 0; generation++;
      $('focus-chain').replaceChildren(...next.chains.map((_, i) => {
        const option = document.createElement('option');
        option.value = i; option.textContent = `Chain ${i + 1}`; return option;
      }));
      $('seek').max = $('progress').max = next.config.iterations;
      $('burn').max = next.config.iterations;
    }, true);
    if (built) guarded(render);
  }
  function dataFor(d) {
    const key = `${generation}/${cursor}/${d.burn}/${d.thin}`;
    if (!diagnosticCache || diagnosticCache.key !== key) {
      const retained = M.retained(result, cursor, d.burn, d.thin);
      diagnosticCache = { key, retained, stats: [0, 1].map(dim => M.diagnostics(retained, dim)) };
    }
    return diagnosticCache;
  }
  function fmt(value) {
    return value === -Infinity ? '−∞' : value === Infinity ? '∞' : U.number(value);
  }
  function pair(values) { return `(${values.map(fmt).join(', ')})`; }
  function marker(frame, point, colour, type = 'filled', radius = 5) {
    const { ctx } = frame, x = frame.xp(point[0]), y = frame.yp(point[1]);
    ctx.strokeStyle = colour; ctx.fillStyle = colour; ctx.lineWidth = 2;
    ctx.beginPath();
    if (type === 'cross' || type === 'plus') {
      const r = radius;
      if (type === 'cross') { ctx.moveTo(x - r, y - r); ctx.lineTo(x + r, y + r); ctx.moveTo(x - r, y + r); ctx.lineTo(x + r, y - r); }
      else { ctx.moveTo(x - r, y); ctx.lineTo(x + r, y); ctx.moveTo(x, y - r); ctx.lineTo(x, y + r); }
      ctx.stroke();
    } else {
      ctx.arc(x, y, radius, 0, 2 * Math.PI);
      if (type === 'open') { ctx.fillStyle = 'white'; ctx.fill(); ctx.stroke(); }
      else { ctx.fill(); ctx.strokeStyle = 'white'; ctx.lineWidth = 1; ctx.stroke(); }
    }
  }
  function focusMove(d) { return cursor ? result.chains[d.focus].moves[cursor - 1] : null; }
  function joint(d) {
    const move = focusMove(d), bounds = result.config.bounds, { axes, nx, ny, logs, maxLog } = reference;
    const ranges = bounds.map((b, dim) => U.extent([...b, result.config.theta[dim], ...(move ? [move.proposal[dim]] : []), ...(result.target.binomial && dim === 1 ? [b[0] - 0.5, b[1] + 0.5] : [])], 0.06));
    const frame = U.plot($('joint-plot'), ranges[0], ranges[1], result.target.labels[0], result.target.labels[1]);
    const { ctx } = frame;
    frame.clip();
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const drop = logs[j * nx + i] - maxLog;
      if (drop < -16) continue;
      ctx.fillStyle = `rgba(8,127,140,${0.06 + 0.67 * Math.exp(drop / 2)})`;
      const left = frame.xp(axes[0].edges[i]), right = frame.xp(axes[0].edges[i + 1]);
      const top = frame.yp(axes[1].edges[j + 1]), bottom = frame.yp(axes[1].edges[j]);
      ctx.fillRect(left, top, right - left + 0.4, bottom - top + 0.4);
    }
    if (d.prior) U.line(frame, [[bounds[0][0], bounds[1][0]], [bounds[0][1], bounds[1][0]],
      [bounds[0][1], bounds[1][1]], [bounds[0][0], bounds[1][1]], [bounds[0][0], bounds[1][0]]], '#7a8698', 1.5, [5, 4]);
    result.chains.forEach((chain, j) => {
      const first = Math.max(0, cursor - d.trail), colour = colours[j];
      const past = chain.states.slice(first, cursor + 1);
      if (d.lines && past.length > 1) {
        const burnEnd = Math.max(0, Math.min(cursor, d.burn) - first + 1);
        if (burnEnd > 1) U.line(frame, past.slice(0, burnEnd), '#8992a188', 1.2);
        U.line(frame, past.slice(Math.max(0, burnEnd - 1)), colour + '88', 1.2);
      }
      const samples = [];
      const count = Math.max(0, Math.ceil((cursor - d.burn) / d.thin)), stride = Math.max(1, Math.ceil(count / 1500));
      for (let t = d.burn + 1; t <= cursor; t += d.thin * stride) samples.push(chain.states[t]);
      U.points(frame, samples, colour + '65', 2);
      marker(frame, chain.states[0], colour, 'cross', 4);
      if (j !== d.focus) marker(frame, chain.states[cursor], colour, 'filled', 4);
    });
    marker(frame, result.config.theta, '#365b99', 'plus', 8);
    if (move) {
      marker(frame, move.before, '#25344b', 'open', 6);
      U.line(frame, [move.before, move.proposal], '#c75800', 2, [5, 4]);
      marker(frame, move.proposal, move.accepted ? '#c75800' : '#9e3f42', 'cross', 6);
      let shown = move.state;
      if (animation && animation.cursor === cursor) {
        const p = Math.min(1, (performance.now() - animation.start) / animation.duration);
        const travel = p < 0.7 ? p / 0.7 : move.accepted ? 1 : (1 - p) / 0.3;
        shown = move.before.map((v, dim) => v + travel * (move.proposal[dim] - v));
      }
      marker(frame, shown, move.accepted ? '#008260' : '#25344b', 'filled', 6);
    } else marker(frame, result.chains[d.focus].states[0], '#25344b', 'open', 6);
    ctx.restore();
  }
  function acceptance(d) {
    const move = focusMove(d);
    if (!move) { U.blank($('acceptance-plot'), 'Advance one iteration to propose and decide.'); return; }
    const alpha = Math.exp(move.logAlpha), u = Math.exp(move.logU);
    const frame = U.plot($('acceptance-plot'), [0, 1], [0, 1], 'Uniform draw U and acceptance probability α', '');
    const { ctx } = frame;
    frame.clip();
    ctx.fillStyle = '#edf1f7';
    ctx.fillRect(frame.xp(0), frame.yp(0.65), frame.xp(1) - frame.xp(0), frame.yp(0.35) - frame.yp(0.65));
    ctx.fillStyle = '#00826055';
    ctx.fillRect(frame.xp(0), frame.yp(0.65), frame.xp(alpha) - frame.xp(0), frame.yp(0.35) - frame.yp(0.65));
    U.line(frame, [[alpha, 0.2], [alpha, 0.8]], '#008260', 2, [4, 3]);
    U.line(frame, [[u, 0.12], [u, 0.88]], '#c75800', 3);
    ctx.fillStyle = move.accepted ? '#008260' : '#9e3f42';
    ctx.font = '600 16px system-ui, sans-serif'; ctx.textAlign = 'center';
    ctx.fillText(move.accepted ? 'ACCEPT' : 'REJECT · stay at previous state', (frame.left + frame.right) / 2, frame.yp(0.94));
    ctx.font = '12px system-ui, sans-serif'; ctx.fillStyle = '#25344b';
    ctx.fillText(`α = ${fmt(alpha)}     U = ${fmt(u)}`, (frame.left + frame.right) / 2, frame.yp(0.04));
    ctx.restore();
  }
  function exploreText(d) {
    const move = focusMove(d), chain = result.chains[d.focus];
    $('progress-text').textContent = `Iteration ${cursor.toLocaleString()} / ${result.config.iterations.toLocaleString()} · chain ${d.focus + 1}`;
    $('seek').value = $('progress').value = cursor;
    $('playback-state').textContent = playing ? 'Playing cached chain' : 'Paused';
    $('data-summary').textContent = `${result.config.observations.toLocaleString()} seeded ${result.config.family} observations. ${result.target.dataSummary}`;
    if (!move) {
      $('step-description').textContent = `Starting at ${pair(chain.states[0])}. No proposals attempted yet. Playback is paused until you press Play or Next iteration.`;
      $('decision-detail').textContent = 'Green shading will show the acceptance region [0, α]. The orange line is the cached uniform draw U. The draw is independent of the proposal.';
      U.table('step-table', [['Current state', pair(chain.states[0])], ['Log unnormalised posterior', fmt(result.target.logDensity(chain.states[0]))]]);
      return;
    }
    $('step-description').textContent = `${move.component}: proposed ${pair(move.proposal)} from ${pair(move.before)}. ${move.accepted ? 'Accepted' : 'Rejected: record the previous state again'}. State ${cursor <= d.burn ? 'is burn-in' : (cursor - d.burn - 1) % d.thin === 0 ? 'is retained' : 'is omitted by display thinning'}.`;
    $('decision-detail').textContent = `Accept if U < α, equivalently log U < min(0, log r). ${move.logTargetProposal === -Infinity ? 'The proposal is outside posterior support, so α = 0.' : result.config.proposal === 'randomwalk' ? 'The symmetric proposal terms cancel.' : 'The non-symmetric proposal terms are included in both directions.'} A rejection still contributes a repeated state.`;
    U.table('step-table', [
      ['Before → proposed', `${pair(move.before)} → ${pair(move.proposal)}`],
      ['Recorded state', pair(move.state)],
      ['log π(before) / log π(proposed)', `${fmt(move.logTargetBefore)} / ${fmt(move.logTargetProposal)}`],
      ['log q(proposed | before)', fmt(move.logQForward)],
      ['log q(before | proposed)', fmt(move.logQReverse)],
      ['Log Hastings ratio log r', fmt(move.logRatio)],
      ['Acceptance probability α', fmt(Math.exp(move.logAlpha))],
      ['Uniform U / log U', `${fmt(Math.exp(move.logU))} / ${fmt(move.logU)}`],
      ['Decision', move.accepted ? 'ACCEPT' : 'REJECT — repeat previous state']
    ]);
  }
  function trace(dim, d) {
    const range = U.extent(result.config.bounds[dim], 0.04);
    const frame = U.plot($(`trace-${dim}`), [0, result.config.iterations], range, 'Attempted iteration', result.target.labels[dim]);
    frame.clip();
    U.line(frame, [[d.burn, range[0]], [d.burn, range[1]]], '#7a8698', 1.2, [4, 4]);
    U.line(frame, [[0, reference.moments[dim].mean], [result.config.iterations, reference.moments[dim].mean]], '#087f8c', 1.2, [5, 3]);
    result.chains.forEach((chain, j) => {
      const points = chain.states.slice(0, cursor + 1).map((state, t) => [t, state[dim]]);
      if (d.lines) {
        U.line(frame, points.slice(0, Math.min(d.burn, cursor) + 1), '#8992a180', 1);
        U.line(frame, points.slice(Math.min(d.burn, cursor)), colours[j] + 'bb', 1.2);
      } else {
        U.points(frame, points.slice(0, Math.min(d.burn, cursor) + 1), '#8992a180', 1.4);
        U.points(frame, points.slice(Math.min(d.burn, cursor) + 1), colours[j] + '99', 1.4);
      }
    });
    frame.ctx.restore();
  }
  function marginal(dim, d, diagnostics) {
    const samples = diagnostics.retained.flatMap(chain => chain.map(x => x[dim]));
    const axis = reference.axes[dim], density = reference.marginal[dim], discrete = result.target.binomial && dim === 1;
    const range = discrete ? [axis.edges[0], axis.edges[axis.edges.length - 1]] : result.config.bounds[dim];
    const bins = discrete ? axis.centers.length : 36;
    const width = (range[1] - range[0]) / bins, histogram = new Array(bins).fill(0);
    for (const x of samples) histogram[Math.min(bins - 1, Math.max(0, Math.floor((x - range[0]) / width)))]++;
    const histDensity = histogram.map(v => samples.length ? v / samples.length / width : 0);
    const ymax = Math.max(...density, ...histDensity, 1e-6) * 1.12;
    const frame = U.plot($(`marginal-${dim}`), range, [0, ymax], result.target.labels[dim], discrete ? 'Probability' : 'Density');
    frame.clip();
    const { ctx } = frame;
    histDensity.forEach((value, i) => {
      ctx.fillStyle = '#365b9950';
      const left = frame.xp(range[0] + i * width), right = frame.xp(range[0] + (i + 1) * width);
      ctx.fillRect(left + 0.5, frame.yp(value), Math.max(0, right - left - 1), frame.yp(0) - frame.yp(value));
    });
    if (discrete) {
      axis.centers.forEach((x, i) => U.line(frame, [[x, 0], [x, density[i]]], '#087f8c', 2));
      U.points(frame, axis.centers.map((x, i) => [x, density[i]]), '#087f8c', 3);
    } else U.line(frame, axis.centers.map((x, i) => [x, density[i]]), '#087f8c', 2.5);
    frame.ctx.restore();
  }
  function inference(d) {
    const diagnostics = dataFor(d);
    for (let dim = 0; dim < 2; dim++) {
      $(`trace-title-${dim}`).textContent = `Trace: ${result.target.labels[dim]}`;
      $(`marginal-title-${dim}`).textContent = `Marginal: ${result.target.labels[dim]}`;
      trace(dim, d); marginal(dim, d, diagnostics);
    }
    const perChain = diagnostics.retained[0].length;
    $('sample-count').textContent = `${diagnostics.stats[0].count.toLocaleString()} retained states (${perChain.toLocaleString()} per chain) from ${cursor.toLocaleString()} revealed attempts per chain. Discard B = ${d.burn}; retain every k = ${d.thin}.`;
    U.table('sample-summary', diagnostics.stats.map((s, dim) => [
      result.target.labels[dim], fmt(s.mean), fmt(s.sd), fmt(reference.moments[dim].mean),
      fmt(reference.moments[dim].sd), fmt(s.rhat), fmt(s.ess)
    ]));
  }
  function dependence(d) {
    const diagnostics = dataFor(d);
    for (let dim = 0; dim < 2; dim++) {
      const acf = diagnostics.stats[dim].acf;
      if (!acf.length || !acf.every(Number.isFinite)) U.blank($(`acf-${dim}`), 'Need at least 3 nonconstant retained states in every chain.');
      else {
        const frame = U.plot($(`acf-${dim}`), [0, Math.max(1, acf.length - 1)], [-1, 1], 'Lag in retained states', 'ACF');
        frame.clip();
        U.line(frame, [[0, 0], [acf.length - 1, 0]], '#8992a1', 1);
        acf.forEach((value, lag) => U.line(frame, [[lag, 0], [lag, value]], '#087f8c', 2));
        frame.ctx.restore();
      }
      const range = U.extent(result.config.bounds[dim], 0.03);
      const frame = U.plot($(`lag-${dim}`), range, range, `${result.target.labels[dim]}[t]`, `${result.target.labels[dim]}[t + k]`);
      frame.clip();
      U.line(frame, [[range[0], range[0]], [range[1], range[1]]], '#8992a1', 1, [4, 4]);
      diagnostics.retained.forEach((chain, j) => {
        U.points(frame, chain.slice(1).map((x, i) => [chain[i][dim], x[dim]]), colours[j] + '80', 2);
      });
      frame.ctx.restore();
    }
    U.table('chain-summary', result.chains.map((chain, j) => {
      const accepted = chain.moves.slice(0, cursor).filter(m => m.accepted).length;
      return [`Chain ${j + 1}`, `${accepted} / ${cursor}`, cursor ? `${fmt(100 * accepted / cursor)}%` : '—', diagnostics.retained[j].length];
    }));
  }
  function explanation() {
    const c = result.config;
    const descriptions = {
      gaussian: 'Each observation is a bivariate Gaussian with unknown means (μx, μy), common known SD σ and known correlation ρ. The posterior is proportional to the product of those bivariate likelihoods and the two independent bounded priors.',
      binomial: 'Each observation is Binomial(N, p), with unknown probability p and integer draw count N. The likelihood includes the N-dependent binomial coefficient: dropping it would change the target. N must be at least the largest observed count.',
      beta: 'Observations lie strictly between 0 and 1 and follow Beta(α, β). Both shapes are inferred. log L = n[log Γ(α + β) − log Γ(α) − log Γ(β)] + (α − 1)Σ log z + (β − 1)Σ log(1 − z).',
      gamma: 'Positive observations follow Gamma(shape k, rate r). Both parameters are inferred. log L = n[k log r − log Γ(k)] + (k − 1)Σ log z − rΣz. The mean is k / r, not k × r.'
    };
    $('model-equations').textContent = `${descriptions[c.family]} The generating pair is ${pair(c.theta)}; it need not equal the posterior mean for this finite data set. Prior: ${c.prior === 'uniform' ? 'uniform on the specified box' : 'independent bounded Gaussians with midpoint means and width/4 SDs (discrete on N)'}.`;
  }
  function render() {
    if (!result) return;
    const d = display();
    exploreText(d);
    explanation();
    if (active === 'explore') { joint(d); acceptance(d); }
    else if (active === 'diagnostics') inference(d);
    else if (active === 'dependence') dependence(d);
  }
  function animate() {
    animationFrame = null;
    if (!animation || !result || !$('error').hidden) return;
    if (active === 'explore') guarded(() => { const d = display(); joint(d); acceptance(d); });
    if (!animation) return;
    if (performance.now() - animation.start < animation.duration) animationFrame = requestAnimationFrame(animate);
    else { animation = null; guarded(render); }
  }
  function seek(position, animated = false) {
    A.check(result && Number.isInteger(position) && position >= 0 && position <= result.config.iterations,
      'Seek position must be an integer inside the cached chain.');
    stopAnimation();
    cursor = position; diagnosticCache = null;
    if (animated && cursor) animation = { cursor, start: performance.now(), duration: display().delay * 0.85 };
    render();
    if (animation) animationFrame = requestAnimationFrame(animate);
  }
  function tick() {
    if (!playing || !result) return;
    guarded(() => {
      if (cursor >= result.config.iterations) { pause(); render(); return; }
      seek(cursor + 1, true);
      if (cursor >= result.config.iterations) { pause(); render(); }
      else timer = setTimeout(tick, display().delay);
    });
  }
  function activate(name, focus = false) {
    active = name;
    tabNames.forEach(tab => {
      $(`tab-${tab}`).setAttribute('aria-selected', String(tab === name));
      $(`tab-${tab}`).tabIndex = tab === name ? 0 : -1;
      $(tab).hidden = tab !== name;
    });
    if (focus) $(`tab-${name}`).focus();
    if ($('error').hidden) guarded(render);
  }
  $('controls').addEventListener('submit', event => event.preventDefault());
  $('display-controls').addEventListener('submit', event => event.preventDefault());
  $('controls').addEventListener('input', () => U.outputs());
  $('controls').addEventListener('change', event => {
    if (event.target.id === 'family') familyControls(true);
    regenerate();
  });
  $('display-controls').addEventListener('input', () => U.outputs());
  $('display-controls').addEventListener('change', () => {
    pause(); stopAnimation(); diagnosticCache = null;
    if (!result) { regenerate(); return; }
    guarded(render);
  });
  $('play').addEventListener('click', () => guarded(() => {
    if (playing) { pause(); stopAnimation(); render(); return; }
    if (cursor >= result.config.iterations) seek(0);
    playing = true; $('play').textContent = 'Pause'; tick();
  }));
  $('step').addEventListener('click', () => guarded(() => { pause(); seek(Math.min(cursor + 1, result.config.iterations), true); }));
  $('rewind').addEventListener('click', () => guarded(() => { pause(); seek(0); }));
  ['finish', 'diagnostic-finish'].forEach(id => $(id).addEventListener('click', () => guarded(() => {
    pause(); seek(result.config.iterations);
  })));
  $('seek').addEventListener('input', () => guarded(() => { pause(); seek(numeric('seek')); }));
  $('seek').addEventListener('change', () => guarded(() => { pause(); seek(numeric('seek')); }));
  tabNames.forEach((name, i) => {
    $(`tab-${name}`).addEventListener('click', () => activate(name));
    $(`tab-${name}`).addEventListener('keydown', event => {
      let next;
      if (event.key === 'ArrowRight' || event.key === 'ArrowDown') next = (i + 1) % tabNames.length;
      else if (event.key === 'ArrowLeft' || event.key === 'ArrowUp') next = (i + tabNames.length - 1) % tabNames.length;
      else if (event.key === 'Home') next = 0;
      else if (event.key === 'End') next = tabNames.length - 1;
      else return;
      event.preventDefault(); activate(tabNames[next], true);
    });
  });
  const resizeObserver = new ResizeObserver(() => {
    if (result && $('error').hidden) guarded(render);
  });
  resizeObserver.observe($('main'));
  window.addEventListener('resize', () => { if (result && $('error').hidden) guarded(render); });
  window.MCMCDemo = Object.freeze({
    snapshot() {
      const move = result && cursor ? result.chains[Number($('focus-chain').value)].moves[cursor - 1] : null;
      return {
        cursor, playing, generation, valid: Boolean(result) && $('error').hidden,
        family: result?.config.family ?? null, iterations: result?.config.iterations ?? null,
        state: result ? result.chains[Number($('focus-chain').value)].states[cursor].slice() : null,
        proposal: move ? move.proposal.slice() : null, accepted: move?.accepted ?? null,
        dataSummary: result?.target.dataSummary ?? null,
        signature: result ? JSON.stringify(result.chains.map(c => [c.states[0], c.states[c.states.length - 1]])) : null
      };
    }
  });
  if (window.matchMedia('(max-width: 800px)').matches) $('settings').open = false;
  try { familyControls(false); regenerate(); }
  finally { document.documentElement.dataset.ready = 'true'; }
})();
