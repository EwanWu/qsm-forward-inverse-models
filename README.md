# Susceptibility-field modelling and QSM inversion

Inspectable numerical models by Yue Wu for the quantitative susceptibility mapping (QSM) measurement chain: susceptibility source, local magnetic field, field estimation, background-field removal and inversion. This repository is **synthetic methodological work**, not the clinical CADASIL processing pipeline or a full implementation of named QSM toolboxes.

## What is implemented

| File | Language | Evidence |
| --- | --- | --- |
| `src/qsm-core.js` | JavaScript | A 3D periodic dipole forward operator, Fourier transforms, synthetic phantoms, phase and inverse-model calculations used by the [interactive QSM chapters](https://mri.ewanarc.com/qsm/). |
| `src/qsm_method_comparisons.py` | Python | Controlled 32 × 32 × 32 comparisons of background-field approaches, structural priors and field orientations, with explicit noise and support conditions. |
| `src/qsm-comparisons.js` | JavaScript | Browser controls for selected comparisons; it expects the EwanArc page structure. |
| `tests/numerical-checks.cjs` | JavaScript | Fifteen numerical checks, including dipole symmetry, FFT round trip, self-adjointness, background-removal and inverse-filter limits. |

Run the short, dependency-free example and checks with Node.js:

```bash
node examples/quick_demo.cjs
node tests/numerical-checks.cjs
```

For the larger synthetic Python comparison, install NumPy and SciPy and run:

```bash
python src/qsm_method_comparisons.py
```

That script creates arrays and summaries in `src/output/`; it can take several minutes. The browser script is provided for code inspection; the live site supplies its page controls.

## Interpretation

These models use periodic boundaries and deliberately controlled conditions. The comparisons illustrate how orientation, noise, finite masks and regularisation alter estimates. “PDF-like”, “RESHARP-like” and morphology-weighted examples are **limited educational implementations**; they must not be cited as exact reproductions of PDF, RESHARP, MEDI, iLSQR, FANSI or clinical pipelines. Neither cellular magnetometry nor histological validation is claimed.

Original code is visible for research-software evaluation; no general reuse licence is granted. Contact Yue Wu at ewan.wu7@gmail.com for adaptation or project-specific materials.
