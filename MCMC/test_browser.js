'use strict';
const assert = require('node:assert/strict');
const { setTimeout: sleep } = require('node:timers/promises');
const { mkdir, rmdir } = require('node:fs/promises');
const path = require('node:path');
const { runBrowser } = require('../TeachingTools/browser_test.js');

async function suite(b) {
  const snapshot = () => b.evaluate('MCMCDemo.snapshot()');
  assert.equal(await b.evaluate('document.title'), 'Metropolis-Hastings MCMC');
  let initial = await snapshot();
  assert.equal(initial.cursor, 0);
  assert.equal(initial.playing, false);
  assert.equal(initial.valid, true);
  assert.ok(await b.hasColour('joint-plot', [37, 52, 75]), 'Initial state / axes are drawn');
  await b.click('step');
  const first = await snapshot();
  assert.equal(first.cursor, 1);
  assert.ok(first.proposal);
  assert.match(await b.text('step-description'), /proposed/);
  await b.click('rewind');
  await b.click('step');
  const replay = await snapshot();
  assert.deepEqual(replay.state, first.state);
  assert.deepEqual(replay.proposal, first.proposal);
  assert.equal(replay.accepted, first.accepted);
  assert.equal(replay.signature, initial.signature);
  await b.click('play');
  await b.waitFor('MCMCDemo.snapshot().cursor >= 3');
  await b.click('play');
  const paused = (await snapshot()).cursor;
  await sleep(750);
  assert.equal((await snapshot()).cursor, paused);
  assert.equal((await snapshot()).playing, false);
  await b.click('finish');
  assert.equal((await snapshot()).cursor, 1000);
  await b.change('seek', 17);
  assert.equal((await snapshot()).cursor, 17);
  const beforeDisplay = await snapshot();
  await b.change('burn', 4);
  await b.change('thin', 3);
  await b.change('lines', false);
  await b.change('show-prior', false);
  assert.equal((await snapshot()).generation, beforeDisplay.generation);
  assert.equal((await snapshot()).signature, beforeDisplay.signature);
  assert.equal((await snapshot()).cursor, 17);
  await b.click('tab-diagnostics');
  assert.match(await b.text('sample-count'), /20 retained states \(5 per chain\)/);
  assert.ok(await b.hasColour('marginal-0', [8, 127, 140]), 'Full posterior marginal is drawn');
  assert.equal(await b.evaluate('document.getElementById("sample-summary").rows.length'), 2);
  await b.click('tab-dependence');
  assert.equal(await b.evaluate('document.getElementById("chain-summary").rows.length'), 4);
  assert.match(await b.text('chain-summary'), /Chain 1/);
  await b.evaluate('document.getElementById("tab-dependence").dispatchEvent(new KeyboardEvent("keydown",{key:"End",bubbles:true}))');
  assert.equal(await b.evaluate('document.getElementById("tab-explain").getAttribute("aria-selected")'), 'true');
  await b.evaluate('document.getElementById("tab-explain").dispatchEvent(new KeyboardEvent("keydown",{key:"Home",bubbles:true}))');
  assert.equal(await b.evaluate('document.getElementById("tab-explore").getAttribute("aria-selected")'), 'true');
  await b.change('iterations', 120);
  assert.equal((await snapshot()).cursor, 0);
  assert.equal((await snapshot()).playing, false);
  assert.equal((await snapshot()).dataSummary, initial.dataSummary);
  for (const family of ['binomial', 'beta', 'gamma', 'gaussian']) {
    await b.change('family', family);
    await b.healthy();
    assert.equal((await snapshot()).family, family);
    await b.change('proposal', 'mixture');
    await b.change('prior', 'gaussian');
    await b.click('finish');
    await b.click('tab-diagnostics');
    await b.healthy();
    assert.match(await b.text('sample-count'), /retained states/);
    await b.click('tab-dependence');
    await b.healthy();
    await b.click('tab-explore');
  }
  await b.change('proposal', 'independence');
  assert.equal(await b.evaluate('document.getElementById("proposal-sd").disabled'), true);
  await b.click('finish');
  await b.healthy();
  const beforeInvalidDisplay = await snapshot();
  await b.change('thin', 0);
  assert.equal(await b.evaluate('document.getElementById("error").hidden'), false);
  assert.equal(await b.evaluate('document.getElementById("play").disabled'), true);
  await b.change('thin', 1);
  await b.healthy();
  assert.equal((await snapshot()).generation, beforeInvalidDisplay.generation);
  assert.equal((await snapshot()).signature, beforeInvalidDisplay.signature);
  await b.change('data-sd', 0);
  assert.equal(await b.evaluate('document.getElementById("error").hidden'), false);
  assert.equal((await snapshot()).state, null);
  assert.equal(await b.text('step-table'), '');
  await b.change('data-sd', 1.2);
  await b.healthy();
  await b.change('chains', 1);
  await b.change('burn', 0);
  await b.click('step');
  await b.click('tab-diagnostics');
  await b.healthy();
  assert.match(await b.text('sample-count'), /1 retained states/);
  await b.call('Emulation.setDeviceMetricsOverride', { width: 390, height: 844, deviceScaleFactor: 2, mobile: true });
  await b.waitFor('document.getElementById("marginal-0").width >= document.getElementById("marginal-0").clientWidth * 1.9');
  assert.equal(await b.evaluate('document.documentElement.scrollWidth <= window.innerWidth + 1'), true);
  await b.navigate();
  assert.equal(await b.evaluate('document.getElementById("settings").open'), false);
  await b.click('step');
  await b.screenshot('mcmc-mobile.png');
  await b.openFile();
  assert.equal((await snapshot()).cursor, 0);
  assert.equal((await snapshot()).playing, false);
  await b.click('step');
  await b.healthy();
  await b.call('Emulation.setDeviceMetricsOverride', { width: 1440, height: 1100, deviceScaleFactor: 1, mobile: false });
  await b.click('finish');
  await b.screenshot('mcmc-desktop.png');
}
async function main() {
  const profileRoot = path.join(__dirname, '.browser-test');
  const created = await mkdir(profileRoot, { recursive: true });
  const previous = process.env.TMPDIR;
  process.env.TMPDIR = profileRoot;
  try {
    await runBrowser('MCMC', suite);
  } finally {
    if (previous === undefined) delete process.env.TMPDIR;
    else process.env.TMPDIR = previous;
    if (created) {
      try { await rmdir(profileRoot); }
      catch (error) { if (error.code !== 'ENOTEMPTY' && error.code !== 'ENOENT') throw error; }
    }
  }
}
main().catch(error => { console.error(error); process.exitCode = 1; });
