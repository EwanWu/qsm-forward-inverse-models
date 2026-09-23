'use strict';
const assert = require('node:assert/strict');
const Q = require('../src/qsm-core.js');
let count = 0;
function test(name, fn) { fn(); count++; console.log('PASS ' + name); }
const maxError = (a, b) => a.reduce((m, v, i) => Math.max(m, Math.abs(v - b[i])), 0);
const dot = (a, b) => a.reduce((s, v, i) => s + v * b[i], 0);
const near = (a, b, eps = 1e-10) => assert.ok(Math.abs(a - b) < eps, `${a} vs ${b}`);
test('3-D complex FFT round trip preserves every voxel', () => {
  const re = Q.gaussianNoise(Q.SIZE, 1), im = Q.gaussianNoise(Q.SIZE, 2), a = re.slice(), b = im.slice();
  Q.fft3(a, b); Q.fft3(a, b, true);
  assert.ok(maxError(a, re) < 4e-15); assert.ok(maxError(b, im) < 4e-15);
});
test('Dipole has zero DC, analytic axial limits, and exact magic-cone zeros', () => {
  const d = Q.dipole(); near(d[0], 0); near(d[Q.idx(1, 0, 0)], 1 / 3); near(d[Q.idx(0, 0, 1)], -2 / 3); near(d[Q.idx(1, 1, 1)], 0);
  const c = new Float64Array(Q.SIZE).fill(2); assert.ok(Math.max(...Q.apply(c, d)) < 1e-14);
});
test('Tilted 3-D dipole remains real and self-adjoint, including Nyquist bins', () => {
  const d = Q.dipole(37), a = Q.gaussianNoise(Q.SIZE, 3), b = Q.gaussianNoise(Q.SIZE, 4);
  for (let z = 0; z < Q.N; z++) for (let y = 0; y < Q.N; y++) for (let x = 0; x < Q.N; x++) near(d[Q.idx(x, y, z)], d[Q.idx((Q.N - x) % Q.N, (Q.N - y) % Q.N, (Q.N - z) % Q.N)]);
  near(dot(a, Q.apply(b, d)), dot(Q.apply(a, d), b), 1e-10);
  const f = Q.transform(a); for (let i = 0; i < Q.SIZE; i++) { f.re[i] *= d[i]; f.im[i] *= d[i]; } Q.fft3(f.re, f.im, true);
  assert.ok(Math.max(...f.im.map(Math.abs)) < 1e-14);
});
test('Unit dipole response equals the 3-D IFFT kernel and sphere center obeys cubic symmetry', () => {
  const impulse = new Float64Array(Q.SIZE); impulse[0] = 1;
  const d = Q.dipole(); assert.ok(maxError(Q.apply(impulse, d), Q.kernel(d)) < 1e-14);
  near(Q.run('forward').stats.centerField, 0, 1e-14);
});
test('Orthonormal 2-D DCT round trip matches Neumann spectral convention', () => {
  const a = Q.gaussianNoise(16 * 16); assert.ok(maxError(a, Q.dct2(Q.dct2(a, 16), 16, true)) < 2e-14);
});
test('Wrapped-gradient LS recovers an Itoh-valid nonperiodic phase up to gauge', () => {
  const n = 32, truth = Float64Array.from({ length: n * n }, (_, i) => 0.35 * (i % n) - 0.24 * Math.floor(i / n) + 0.4 * Math.sin((i % n) / 7));
  const wrapped = Float64Array.from(truth, Q.wrap), unwrapped = Q.alignTo(Q.leastSquaresUnwrap(wrapped, n), truth);
  assert.ok(maxError(unwrapped, truth) < 2e-12);
});
test('Shared noisy phase comparison has real path failures and finite LS errors', () => {
  const r = Q.run('unwrap', { te: 40, noise: 0.7 });
  assert.ok(r.stats.pathRmse > 1); assert.ok(r.stats.leastSquaresRmse < r.stats.pathRmse);
  const clean = Q.run('unwrap', { noise: 0 }); assert.ok(clean.stats.pathRmse < 1e-12); assert.ok(clean.stats.leastSquaresRmse < 1e-12);
});
test('Unwrap local response is mean-gauged identity when no edge branch changes', () => {
  const r = Q.run('unwrap', { noise: 0 }); assert.equal(r.stats.branchChanged, false);
  for (const c of r.curves.filter(c => c.unit === '1')) { near(c.data[32], 1 - 1 / 4096, 1e-9); near(c.data[0], -1 / 4096, 1e-9); }
  // Place a known center edge half an epsilon above a -π branch. The fixed
  // positive central perturbation must then change that edge's branch.
  const phase = r.maps.find(m => m.key === 'truePhase').data, eps = Q.gaussianNoise(64 * 64), c = 32 * 64 + 32;
  let chosenNoise = null;
  for (const j of [c - 1, c + 1, c - 64, c + 64]) for (let k = -2; k <= 2; k++) {
    const n = (-Math.PI + r.stats.epsilon / 2 + 2 * Math.PI * k - phase[j] + phase[c]) / (eps[j] - eps[c]);
    if (n > 0 && n < 1.4) chosenNoise = n;
  }
  assert.notEqual(chosenNoise, null);
  const crossing = Q.run('unwrap', { noise: chosenNoise }); assert.equal(crossing.stats.branchChanged, true); assert.ok(crossing.stats.branchChanges >= 1);
});
test('Phase simulation uses the documented negative receive sign and correct ppm/time units', () => {
  const r = Q.run('unwrap', { te: 30, noise: 0 }), n = 64, x = 20, z = 28;
  const X = (x - (n - 1) / 2) / (n / 2), Z = (z - (n - 1) / 2) / (n / 2);
  const b = 0.24 * X + 0.14 * Z + 0.16 * Math.exp(-((X + 0.3) ** 2 + (Z - 0.1) ** 2) / 0.13) - 0.13 * Math.exp(-((X - 0.3) ** 2 + (Z + 0.2) ** 2) / 0.07);
  near(r.maps.find(m => m.key === 'truePhase').data[z * n + x], -2 * Math.PI * 42.58e6 * 3 * 0.03 * b * 1e-6);
  assert.match(r.meta.model, /Negative receive-phase/);
});
test('Discrete 3-D ball is normalized and annihilates harmonic polynomials on eroded mask', () => {
  const smv = Q.sphericalMean(3), mask = Q.sphereMask(13), valid = Q.sphereMask(10), field = new Float64Array(Q.SIZE);
  near(smv.spatial.reduce((s, v) => s + v, 0), 1);
  for (let z = 0; z < Q.N; z++) for (let y = 0; y < Q.N; y++) for (let x = 0; x < Q.N; x++) { const X = (x - 16) / 13, Y = (y - 16) / 13, Z = (z - 16) / 13; field[Q.idx(x, y, z)] = mask[Q.idx(x, y, z)] * (X * X - Y * Y + Z * (2 * Z * Z - 3 * X * X - 3 * Y * Y)); }
  const hp = Q.apply(field, smv.highpass); for (let i = 0; i < Q.SIZE; i++) if (valid[i]) near(hp[i], 0, 1e-14);
});
test('Harmonic projection recovers all supported quadratic harmonics exactly', () => {
  const mask = Q.sphereMask(13), inner = Q.sphereMask(11), shell = Float64Array.from(mask, (v, i) => v - inner[i]), field = new Float64Array(Q.SIZE);
  for (let z = 0; z < Q.N; z++) for (let y = 0; y < Q.N; y++) for (let x = 0; x < Q.N; x++) field[Q.idx(x, y, z)] = dot([0.2, -0.1, 0.3, 0.4, 0.2, -0.3, 0.1, 0.15, -0.3], Q.harmonicBasis((x - 16) / 13, (y - 16) / 13, (z - 16) / 13));
  const local = Q.harmonicProjection(field, shell).local; assert.ok(Math.max(...local.map(Math.abs)) < 1e-13);
});
test('Masked regularized SMV inversion suppresses pure harmonic field', () => {
  const mask = Q.sphereMask(13), valid = Q.sphereMask(10), field = new Float64Array(Q.SIZE);
  for (let z = 0; z < Q.N; z++) for (let y = 0; y < Q.N; y++) for (let x = 0; x < Q.N; x++) field[Q.idx(x, y, z)] = (x - 16) * (z - 16) / 100;
  const r = Q.sharpSolve(field, mask, valid, Q.sphericalMean(3).highpass, 0.05); assert.ok(Math.max(...r.x.map(Math.abs)) < 1e-12);
});
test('Background comparison uses a common eroded mask and a converged solve', () => {
  const r = Q.run('background'); assert.ok(r.stats.relativeResidual < 1e-6);
  assert.ok(r.stats.retainedVoxels < r.stats.originalVoxels); assert.equal(r.curves.length, 2);
  assert.ok(r.stats.sharpRmse < 0.04); assert.ok(r.stats.harmonicRmse < 0.04);
});
test('Each inverse kernel predicts its noiseless reconstruction exactly', () => {
  const d = Q.dipole(23), chi = Q.phantom(), field = Q.apply(chi, d), filters = Q.inverseFilters(d, 0.15, 0.04);
  for (const g of Object.values(filters)) {
    const transfer = Float64Array.from(d, (v, i) => v * g[i]);
    assert.ok(maxError(Q.apply(field, g), Q.apply(chi, transfer)) < 5e-15);
    near(g[0], 0); near(Q.kernel(transfer).reduce((s, v) => s + v, 0), 0, 1e-13);
  }
});
test('Regularization lowers white-noise amplification with identical input realization', () => {
  const low = Q.run('inversion', { threshold: 0.08, lambda: 0.02 }), high = Q.run('inversion', { threshold: 0.25, lambda: 0.15 });
  assert.ok(maxError(low.maps.find(m => m.key === 'noisyField').data, high.maps.find(m => m.key === 'noisyField').data) === 0);
  for (const name of ['tkd', 'tikhonov', 'gradient']) assert.ok(high.stats[name + 'NoiseGain'] < low.stats[name + 'NoiseGain']);
});
console.log(`${count} numerical checks passed.`);
