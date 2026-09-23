/* QSM teaching laboratory. Genuine 3-D periodic Fourier operators; no dependencies. */
(function (root, factory) {
  'use strict';
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.QSM = api;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';
  const N = 32, SIZE = N * N * N, TAU = 2 * Math.PI;
  const idx = (x, y, z, n = N) => (z * n + y) * n + x;
  const centered = (i, n = N) => i < n / 2 ? i : i - n;
  const wrap = x => Math.atan2(Math.sin(x), Math.cos(x));
  const sum = a => a.reduce((s, x) => s + x, 0);
  const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
  const fftPlans = new Map();
  function plan(n) {
    if (fftPlans.has(n)) return fftPlans.get(n);
    if (n < 2 || (n & (n - 1))) throw Error('FFT length must be a power of two.');
    const rev = new Uint32Array(n), c = new Float64Array(n / 2), s = new Float64Array(n / 2);
    const bits = Math.round(Math.log2(n));
    for (let i = 0; i < n; i++) { let a = i, b = 0; for (let j = 0; j < bits; j++) { b = b * 2 + (a & 1); a >>>= 1; } rev[i] = b; }
    for (let i = 0; i < n / 2; i++) { c[i] = Math.cos(TAU * i / n); s[i] = Math.sin(TAU * i / n); }
    const p = { rev, c, s }; fftPlans.set(n, p); return p;
  }
  function fft1(re, im, inverse = false) {
    const n = re.length, p = plan(n);
    for (let i = 0; i < n; i++) { const j = p.rev[i]; if (j > i) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; } }
    for (let len = 2; len <= n; len *= 2) {
      const half = len / 2, step = n / len;
      for (let start = 0; start < n; start += len) for (let j = 0; j < half; j++) {
        const a = start + j, b = a + half, k = j * step, wi = inverse ? p.s[k] : -p.s[k];
        const tr = p.c[k] * re[b] - wi * im[b], ti = p.c[k] * im[b] + wi * re[b];
        re[b] = re[a] - tr; im[b] = im[a] - ti; re[a] += tr; im[a] += ti;
      }
    }
    if (inverse) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
  }
  function fft3(re, im, inverse = false, n = N) {
    if (re.length !== n ** 3 || im.length !== re.length) throw Error('Invalid FFT volume.');
    const r = new Float64Array(n), q = new Float64Array(n);
    for (let axis = 0; axis < 3; axis++) for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) {
      for (let j = 0; j < n; j++) { const i = axis === 0 ? idx(j, a, b, n) : axis === 1 ? idx(a, j, b, n) : idx(a, b, j, n); r[j] = re[i]; q[j] = im[i]; }
      fft1(r, q, inverse);
      for (let j = 0; j < n; j++) { const i = axis === 0 ? idx(j, a, b, n) : axis === 1 ? idx(a, j, b, n) : idx(a, b, j, n); re[i] = r[j]; im[i] = q[j]; }
    }
    return { re, im };
  }
  function transform(a) { return fft3(Float64Array.from(a), new Float64Array(a.length)); }
  function apply(a, multiplier) {
    const f = transform(a);
    for (let i = 0; i < a.length; i++) { f.re[i] *= multiplier[i]; f.im[i] *= multiplier[i]; }
    return fft3(f.re, f.im, true).re;
  }
  function kernel(multiplier) { return fft3(Float64Array.from(multiplier), new Float64Array(SIZE), true).re; }
  function dipole(angle = 0) {
    const a = angle * Math.PI / 180, bx = Math.sin(a), bz = Math.cos(a), d = new Float64Array(SIZE);
    // At an even-grid Nyquist component, +N/2 and -N/2 are identical bins.
    // Average the continuum symbol over these equivalent signs. This preserves
    // a real self-adjoint discrete operator for tilted B0, including that plane.
    for (let z = 0; z < N; z++) for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      const kx = centered(x), ky = centered(y), kz = centered(z), k2 = kx * kx + ky * ky + kz * kz;
      const cross = x === N / 2 || z === N / 2 ? 0 : 2 * bx * bz * kx * kz;
      d[idx(x, y, z)] = k2 ? 1 / 3 - (bx * bx * kx * kx + bz * bz * kz * kz + cross) / k2 : 0;
    }
    return d;
  }
  function gaussianNoise(size, seed = 7183) {
    let state = seed >>> 0;
    const u = () => { state = (1664525 * state + 1013904223) >>> 0; return (state + 0.5) / 4294967296; };
    const a = new Float64Array(size);
    for (let i = 0; i < size; i += 2) { const r = Math.sqrt(-2 * Math.log(u())), t = TAU * u(); a[i] = r * Math.cos(t); if (i + 1 < size) a[i + 1] = r * Math.sin(t); }
    return a;
  }
  function phantom(chi = 0.6, single = false) {
    const a = new Float64Array(SIZE);
    for (let z = 0; z < N; z++) for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      const X = x - N / 2, Y = y - N / 2, Z = z - N / 2;
      let v = 0;
      if (single) { if (X * X + Y * Y + Z * Z <= 4 ** 2) v = chi; }
      else {
        if (((X + 4) / 4.5) ** 2 + (Y / 5) ** 2 + ((Z + 1) / 6) ** 2 <= 1) v += chi;
        if (((X - 5) / 2.5) ** 2 + (Y / 3.5) ** 2 + ((Z - 3) / 3.5) ** 2 <= 1) v -= chi * 0.6;
        if (((X - 3) / 2) ** 2 + ((Y + 1) / 2) ** 2 + ((Z + 6) / 2) ** 2 <= 1) v += chi * 0.35;
      }
      a[idx(x, y, z)] = v;
    }
    return a;
  }
  function sphereMask(radius) {
    const a = new Float64Array(SIZE);
    for (let z = 0; z < N; z++) for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) a[idx(x, y, z)] = Math.hypot(x - N / 2, y - N / 2, z - N / 2) <= radius ? 1 : 0;
    return a;
  }
  function slice(a, shifted = false, mask = null) {
    const b = new Float64Array(N * N);
    for (let z = 0; z < N; z++) for (let x = 0; x < N; x++) {
      const i = shifted ? idx((x + N / 2) % N, 0, (z + N / 2) % N) : idx(x, N / 2, z);
      b[z * N + x] = mask && !mask[i] ? NaN : a[i];
    }
    return b;
  }
  function profile(a) { return Float64Array.from({ length: N }, (_, j) => a[idx(0, 0, (j + N / 2) % N)]); }
  function map(key, a, unit = 'ppm', range = null, shifted = false, mask = null) { return { key, data: slice(a, shifted, mask), width: N, height: N, unit, range }; }
  function error(a, truth, mask = null, align = true) {
    let offset = 0, count = 0;
    for (let i = 0; i < a.length; i++) if (!mask || mask[i]) { offset += a[i] - truth[i]; count++; }
    offset = align ? offset / count : 0;
    let s = 0;
    for (let i = 0; i < a.length; i++) if (!mask || mask[i]) s += (a[i] - truth[i] - offset) ** 2;
    return { rmse: Math.sqrt(s / count), offset, count };
  }
  function alignTo(a, truth, mask = null) { const e = error(a, truth, mask); return Float64Array.from(a, x => x - e.offset); }
  function forward(p) {
    const chi = p.chi ?? 0.6, angle = p.angle ?? 0, object = phantom(chi, true), d = dipole(angle), field = apply(object, d), h = kernel(d);
    return {
      maps: [map('susceptibility', object, 'ppm', [-chi, chi]), map('localField', field, 'ppm', [-chi / 3, chi / 3]), map('dipoleSpectrum', d, '1', [-2 / 3, 1 / 3], true), map('dipoleKernel', h, '1', [-0.12, 0.12], true)],
      curves: [{ key: 'dipolePsf', data: profile(h), x: Float64Array.from({ length: N }, (_, i) => i - N / 2), unit: '1' }],
      stats: { chi, angle, meanField: sum(field) / SIZE, centerField: field[idx(16, 16, 16)], n: N },
      meta: { boundary: 'periodic 32³ cube; central x–z slices; B0 lies in the x–z plane', normalization: 'forward FFT unnormalized, inverse FFT / N³; D(0)=0; susceptibility and relative field both ppm', nyquist: 'Tilted symbols are averaged over equivalent signs of even-grid Nyquist components to preserve Hermitian symmetry.' }
    };
  }
  const dctPlans = new Map();
  function dctPlan(n) {
    if (dctPlans.has(n)) return dctPlans.get(n);
    const c = new Float64Array(n * n);
    for (let k = 0; k < n; k++) for (let x = 0; x < n; x++) c[k * n + x] = Math.sqrt((k ? 2 : 1) / n) * Math.cos(Math.PI * (x + 0.5) * k / n);
    dctPlans.set(n, c); return c;
  }
  function dct2(a, n, inverse = false) {
    const c = dctPlan(n), t = new Float64Array(a.length), b = new Float64Array(a.length);
    for (let y = 0; y < n; y++) for (let k = 0; k < n; k++) for (let x = 0; x < n; x++) t[y * n + k] += a[y * n + x] * (inverse ? c[x * n + k] : c[k * n + x]);
    for (let x = 0; x < n; x++) for (let k = 0; k < n; k++) for (let y = 0; y < n; y++) b[k * n + x] += t[y * n + x] * (inverse ? c[y * n + k] : c[k * n + y]);
    return b;
  }
  function leastSquaresUnwrap(w, n) {
    const rhs = new Float64Array(w.length);
    // A is the forward edge difference. rhs=Aᵀ g, no edges across the boundary.
    for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
      const i = y * n + x;
      if (x + 1 < n) { const g = wrap(w[i + 1] - w[i]); rhs[i] -= g; rhs[i + 1] += g; }
      if (y + 1 < n) { const g = wrap(w[i + n] - w[i]); rhs[i] -= g; rhs[i + n] += g; }
    }
    const f = dct2(rhs, n);
    for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
      const eigenvalue = 4 - 2 * Math.cos(Math.PI * x / n) - 2 * Math.cos(Math.PI * y / n);
      f[y * n + x] = eigenvalue > 1e-14 ? f[y * n + x] / eigenvalue : 0;
    }
    return dct2(f, n, true);
  }
  function pathUnwrap(w, n) {
    const a = new Float64Array(w.length); a[0] = w[0];
    for (let y = 0; y < n; y++) {
      if (y) a[y * n] = a[(y - 1) * n] + wrap(w[y * n] - w[(y - 1) * n]);
      for (let x = 1; x < n; x++) { const i = y * n + x; a[i] = a[i - 1] + wrap(w[i] - w[i - 1]); }
    }
    return a;
  }
  function unwrapRun(p) {
    const n = 64, te = p.te ?? 30, noise = p.noise ?? 0.25, b0 = 3, a = new Float64Array(n * n), w = new Float64Array(n * n), eps = gaussianNoise(n * n);
    // Smooth synthetic frequency field, in ppm, sampled on a nonperiodic square.
    for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
      const X = (x - (n - 1) / 2) / (n / 2), Z = (y - (n - 1) / 2) / (n / 2), i = y * n + x;
      const field = 0.24 * X + 0.14 * Z + 0.16 * Math.exp(-((X + 0.3) ** 2 + (Z - 0.1) ** 2) / 0.13) - 0.13 * Math.exp(-((X - 0.3) ** 2 + (Z + 0.2) ** 2) / 0.07);
      a[i] = -TAU * 0.04258 * b0 * te * field; w[i] = wrap(a[i] + noise * eps[i]);
    }
    const row = alignTo(pathUnwrap(w, n), a), ls = alignTo(leastSquaresUnwrap(w, n), a), er = error(row, a), el = error(ls, a);
    const epsilon = 0.01, center = (n / 2) * n + n / 2, perturbed = w.slice(); perturbed[center] = wrap(perturbed[center] + epsilon);
    const rowPerturbed = pathUnwrap(perturbed, n), lsPerturbed = leastSquaresUnwrap(perturbed, n);
    const differential = (after, before) => { const gaugeChange = (sum(after) - sum(before)) / (n * n); return Float64Array.from(after, (v, i) => (v - before[i] - gaugeChange) / epsilon); };
    const rowResponse = differential(rowPerturbed, row), lsResponse = differential(lsPerturbed, ls);
    let branchChanges = 0;
    for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
      const i = y * n + x;
      for (const j of [x + 1 < n ? i + 1 : -1, y + 1 < n ? i + n : -1]) if (j >= 0) {
        const expected = epsilon * ((j === center ? 1 : 0) - (i === center ? 1 : 0));
        if (Math.abs(wrap(perturbed[j] - perturbed[i]) - wrap(w[j] - w[i]) - expected) > 1e-6) branchChanges++;
      }
    }
    let max = 1; for (const v of a) max = Math.max(max, Math.abs(v));
    const mm = (key, data, range) => ({ key, data, width: n, height: n, unit: 'rad', range });
    return { maps: [mm('truePhase', a, [-max, max]), mm('wrappedPhase', w, [-Math.PI, Math.PI]), mm('pathPhase', row, [-max, max]), mm('leastSquaresPhase', ls, [-max, max]), mm('pathError', Float64Array.from(row, (v, i) => v - a[i]), [-Math.PI, Math.PI]), mm('leastSquaresError', Float64Array.from(ls, (v, i) => v - a[i]), [-Math.PI, Math.PI])],
      curves: [{ key: 'truePhase', data: a.slice(n * n / 2, n * n / 2 + n), unit: 'rad' }, { key: 'pathPhase', data: row.slice(n * n / 2, n * n / 2 + n), unit: 'rad' }, { key: 'leastSquaresPhase', data: ls.slice(n * n / 2, n * n / 2 + n), unit: 'rad' }, { key: 'pathResponse', data: rowResponse.slice(n * n / 2, n * n / 2 + n), x: Float64Array.from({ length: n }, (_, i) => i - n / 2), unit: '1' }, { key: 'leastSquaresResponse', data: lsResponse.slice(n * n / 2, n * n / 2 + n), x: Float64Array.from({ length: n }, (_, i) => i - n / 2), unit: '1' }],
      stats: { pathRmse: er.rmse, leastSquaresRmse: el.rmse, te, noise, b0, epsilon, branchChanges, branchChanged: branchChanges > 0 }, meta: { boundary: 'Nonperiodic square, Neumann graph Laplacian (no boundary-crossing edges), DCT diagonalization; gauges aligned to truth for evaluation only.', response: 'Local finite-difference response about the displayed wrapped phase/noise: add 0.01 rad at center before wrapping, repeat estimator, remove the mean from each estimate, and divide their difference by 0.01. Central row shown. It is a local derivative only if the edge branch pattern stays unchanged; branchChanges counts changed edges. Phase unwrapping has no global PSF.', model: 'Negative receive-phase convention: phi=-2π gamma_bar B0 TE b; additive phase Gaussian noise before wrapping.' } };
  }
  function harmonicBasis(x, y, z) { return [1, x, y, z, x * y, x * z, y * z, x * x - y * y, 2 * z * z - x * x - y * y]; }
  function solveDense(matrix, rhs) {
    const n = rhs.length, a = matrix.map((r, i) => [...r, rhs[i]]);
    for (let k = 0; k < n; k++) {
      let p = k; for (let j = k + 1; j < n; j++) if (Math.abs(a[j][k]) > Math.abs(a[p][k])) p = j;
      [a[k], a[p]] = [a[p], a[k]];
      if (Math.abs(a[k][k]) < 1e-14) throw Error('Singular harmonic basis.');
      const pivot = a[k][k]; for (let j = k; j <= n; j++) a[k][j] /= pivot;
      for (let i = 0; i < n; i++) if (i !== k) { const f = a[i][k]; for (let j = k; j <= n; j++) a[i][j] -= f * a[k][j]; }
    }
    return a.map(row => row[n]);
  }
  function harmonicProjection(field, fitMask) {
    const G = Array.from({ length: 9 }, () => new Float64Array(9)), r = new Float64Array(9);
    for (let z = 0; z < N; z++) for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      const i = idx(x, y, z); if (!fitMask[i]) continue;
      const v = harmonicBasis((x - 16) / 13, (y - 16) / 13, (z - 16) / 13);
      for (let j = 0; j < 9; j++) { r[j] += v[j] * field[i]; for (let k = 0; k < 9; k++) G[j][k] += v[j] * v[k]; }
    }
    const c = solveDense(G, r), bg = new Float64Array(SIZE);
    for (let z = 0; z < N; z++) for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) bg[idx(x, y, z)] = dot(c, harmonicBasis((x - 16) / 13, (y - 16) / 13, (z - 16) / 13));
    return { background: bg, local: Float64Array.from(field, (v, i) => v - bg[i]), coefficients: c };
  }
  function sphericalMean(radius) {
    const s = new Float64Array(SIZE); let count = 0;
    for (let z = -radius; z <= radius; z++) for (let y = -radius; y <= radius; y++) for (let x = -radius; x <= radius; x++) if (x * x + y * y + z * z <= radius * radius) { s[idx((x + N) % N, (y + N) % N, (z + N) % N)] = 1; count++; }
    for (let i = 0; i < SIZE; i++) s[i] /= count;
    const f = transform(s), h = Float64Array.from(f.re, v => 1 - v); h[0] = 0;
    return { spatial: s, highpass: h, count };
  }
  function cg(operator, rhs, maxIterations = 120, tolerance = 1e-6) {
    const x = new Float64Array(rhs.length), r = Float64Array.from(rhs), p = Float64Array.from(rhs); let rr = dot(r, r), initial = rr, iteration = 0;
    if (initial === 0) return { x, iterations: 0, relativeResidual: 0 };
    for (; iteration < maxIterations; iteration++) {
      const Ap = operator(p), denominator = dot(p, Ap); if (denominator <= 1e-25) break;
      const alpha = rr / denominator;
      for (let i = 0; i < x.length; i++) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
      const next = dot(r, r);
      if (next / initial < tolerance ** 2) { rr = next; iteration++; break; }
      const beta = next / rr; for (let i = 0; i < x.length; i++) p[i] = r[i] + beta * p[i]; rr = next;
    }
    return { x, iterations: iteration, relativeResidual: Math.sqrt(rr / initial) };
  }
  function sharpSolve(field, mask, valid, highpass, alpha) {
    const masked = (a, m) => Float64Array.from(a, (v, i) => v * m[i]);
    const A = a => masked(apply(masked(a, mask), highpass), valid);
    const AT = a => masked(apply(masked(a, valid), highpass), mask);
    const rhs = AT(A(field));
    const solution = cg(a => { const b = AT(A(a)); for (let i = 0; i < SIZE; i++) b[i] += alpha * alpha * a[i]; return b; }, rhs);
    return { ...solution, highpassed: A(field) };
  }
  function background(p) {
    const strength = p.strength ?? 1, radius = Math.max(2, Math.min(5, Math.round(p.radius ?? 3))), alpha = Math.max(0.005, p.alpha ?? 0.05);
    const mask = sphereMask(13), valid = sphereMask(13 - radius), shell = Float64Array.from(mask, (v, i) => v && !sphereMaskInnerCache[i] ? 1 : 0);
    const object = phantom(0.5), local = apply(object, dipole(0)), bg = new Float64Array(SIZE), total = new Float64Array(SIZE);
    for (let z = 0; z < N; z++) for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      const X = (x - 16) / 13, Y = (y - 16) / 13, Z = (z - 16) / 13, i = idx(x, y, z);
      bg[i] = strength * (0.23 * Z + 0.12 * X + 0.18 * (X * X - Y * Y) + 0.12 * Z * (2 * Z * Z - 3 * X * X - 3 * Y * Y));
      total[i] = local[i] + bg[i];
    }
    const poly = harmonicProjection(total, shell), smv = sphericalMean(radius), sharp = sharpSolve(total, mask, valid, smv.highpass, alpha);
    const impulse = new Float64Array(SIZE); impulse[idx(16, 16, 16)] = 1;
    const polyImpulse = harmonicProjection(impulse, shell).local, sharpImpulse = sharpSolve(impulse, mask, valid, smv.highpass, alpha).x;
    const spatialProfile = a => Float64Array.from({ length: N }, (_, z) => valid[idx(16, 16, z)] ? a[idx(16, 16, z)] : NaN);
    const polyAligned = alignTo(poly.local, local, valid), sharpAligned = alignTo(sharp.x, local, valid);
    const pe = error(poly.local, local, valid), se = error(sharp.x, local, valid), range = [-0.3, 0.3];
    return { maps: [map('totalField', total, 'ppm', [-0.6, 0.6], false, mask), map('trueBackground', bg, 'ppm', [-0.6, 0.6], false, mask), map('trueLocalField', local, 'ppm', range, false, valid), map('harmonicLocalField', polyAligned, 'ppm', range, false, valid), map('sharpLocalField', sharpAligned, 'ppm', range, false, valid), map('sharpHighpass', sharp.highpassed, 'ppm', range, false, valid)],
      curves: [{ key: 'harmonicResponse', data: spatialProfile(polyImpulse), x: Float64Array.from({ length: N }, (_, i) => i - 16), unit: '1' }, { key: 'sharpResponse', data: spatialProfile(sharpImpulse), x: Float64Array.from({ length: N }, (_, i) => i - 16), unit: '1' }], stats: { harmonicRmse: pe.rmse, sharpRmse: se.rmse, retainedVoxels: se.count, originalVoxels: sum(mask), radius, alpha, iterations: sharp.iterations, relativeResidual: sharp.relativeResidual },
      meta: { harmonic: 'Least-squares fit of all nine degree ≤2 solid harmonics to the 11<r≤13 voxel outer shell. True background also contains a cubic harmonic; local field can contaminate this fit.', sharp: 'Educational RESHARP-like masked solve: min_x ||Mv (I-Sr) M x - Mv (I-Sr) M b_total||² + alpha² ||x||². Sr is the normalized discrete 3-D ball mean; Mv is the radius-eroded mask. Conjugate gradients solves the normal equations. Displayed results and RMSE share Mv, with only a constant reference aligned to truth.', response: 'Curves are central z profiles of each raw linear estimator applied to a unit field impulse at the geometric center. Masks break translation invariance; these are center-specific responses, not global PSFs. Curves are shown before truth-based constant alignment.', boundary: 'Periodic FFT used to apply the ball convolution; radius-eroded data never sample outside the original spherical mask. Reconstruction is supported on M. This is a finite-grid teaching model, not a clinical SHARP implementation.' } };
  }
  const sphereMaskInnerCache = sphereMask(11);
  function inverseFilters(d, threshold, lambda) {
    const tkd = new Float64Array(SIZE), zero = new Float64Array(SIZE), grad = new Float64Array(SIZE);
    for (let z = 0; z < N; z++) for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      const i = idx(x, y, z), v = d[i], L = 4 * (Math.sin(Math.PI * x / N) ** 2 + Math.sin(Math.PI * y / N) ** 2 + Math.sin(Math.PI * z / N) ** 2);
      tkd[i] = Math.abs(v) > 1e-14 ? Math.sign(v) / Math.max(Math.abs(v), threshold) : 0;
      zero[i] = v / (v * v + lambda * lambda);
      grad[i] = i ? v / (v * v + lambda * lambda * L) : 0;
    }
    return { tkd, tikhonov: zero, gradient: grad };
  }
  function inversion(p) {
    const threshold = Math.max(0.01, p.threshold ?? 0.15), lambda = Math.max(0.001, p.lambda ?? 0.04), noise = p.noise ?? 0.005, angle = p.angle ?? 0;
    const object = phantom(0.6), d = dipole(angle), field = apply(object, d), eps = gaussianNoise(SIZE), observed = Float64Array.from(field, (v, i) => v + noise * eps[i]), filters = inverseFilters(d, threshold, lambda), mask = sphereMask(13);
    const meanChi = sum(object) / SIZE, truth = Float64Array.from(object, v => v - meanChi), maps = [map('trueSusceptibility', truth, 'ppm', [-0.4, 0.6]), map('noisyField', observed, 'ppm', [-0.25, 0.25])], curves = [], stats = { threshold, lambda, noise, angle, meanChi };
    for (const name of ['tkd', 'tikhonov', 'gradient']) {
      const rec = apply(observed, filters[name]), transfer = Float64Array.from(d, (v, i) => v * filters[name][i]), h = kernel(transfer), residual = apply(rec, d);
      maps.push(map(name, rec, 'ppm', [-0.4, 0.6]));
      curves.push({ key: name + 'Psf', data: profile(h), x: Float64Array.from({ length: N }, (_, i) => i - N / 2), unit: '1' });
      stats[name + 'Rmse'] = error(rec, truth, mask, false).rmse;
      stats[name + 'Residual'] = error(residual, observed, null, false).rmse;
      stats[name + 'NoiseGain'] = Math.sqrt(dot(filters[name], filters[name]) / SIZE);
    }
    return { maps, curves, stats, meta: { objective: 'TKD: G=sign(D)/max(|D|,threshold), G=0 at exact zeros. Zero-order: G=D/(D²+lambda²). Gradient: G=D/(D²+lambda² L), L=4 sum_j sin²(pi k_j/N). All DC bins fixed to zero.', psf: 'Exact global linear resolution kernel h=IFFT3(GD), central z line (x=y=0). Curves retain absolute amplitudes. Reconstructed maps are central x–z slices.', evaluation: 'All methods use the identical 3-D field/noise realization. Truth is referenced to zero mean over the full periodic cube; RMSE uses the same r≤13 sphere without additional method-specific offset fitting. NoiseGain is the exact global white-noise amplification sqrt(mean(G²)).', boundary: 'Periodic 32³ box; no brain masking or data weighting inside these diagonal inverse models.' } };
  }
  function run(name, p = {}) {
    if (name === 'forward') return forward(p);
    if (name === 'unwrap') return unwrapRun(p);
    if (name === 'background') return background(p);
    if (name === 'inversion') return inversion(p);
    throw Error('Unknown QSM experiment: ' + name);
  }
  return { N, SIZE, idx, centered, wrap, fft1, fft3, transform, apply, kernel, dipole, gaussianNoise, phantom, sphereMask, slice, profile, error, alignTo, dct2, leastSquaresUnwrap, pathUnwrap, harmonicBasis, harmonicProjection, sphericalMean, sharpSolve, inverseFilters, run };
});
