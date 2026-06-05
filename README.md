# SusLifePath Paper 0 — Code and Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20554998.svg)](https://doi.org/10.5281/zenodo.20554998)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Companion repository for the manuscript

> Toson, F. *A simplified engineering algorithm for collision risk assessment and classification of LEO satellites.* Acta Astronautica, submitted 2026.

The MATLAB toolchain reproduces all numerical results, figures and tables of the paper from a May 2026 Celestrak General Perturbations catalogue snapshot, including the N=100 statistical validation and the snapshot-ECOB proxy benchmark added in the camera-ready revision.

## Contents

```
suslifepath-paper0/
├── code/                              MATLAB toolchain
│   ├── databasecreator_real.m         build resident population from Celestrak CSV
│   ├── individual_probability_flux.m  P_ind = rho_all * sigma * v_rel * T_op
│   ├── collective_probability.m       P_col = rho_op * sigma * v_rel * T_tot * f_frag
│   ├── ecob_proxy.m                   snapshot-ECOB proxy: mass-based fragmentation,
│   │                                  200-yr horizon, same density snapshot
│   ├── risk_index.m                   normalise + classify (5-class ECSS-compatible)
│   ├── paper0_casestudies.m           Table 2: three case-study indices
│   ├── paper0_figures.m               Figure 2 (altitude distribution) + Figure 3 (asymmetry bars)
│   ├── paper0_validation.m            Section 4.2: statistical validation on N=100 random LEO,
│   │                                  Figure 4 + Table 3 (top-10 asymmetry)
│   ├── paper0_parametric.m            Section 4.3: Figure 5 + Table 4 (cost thresholds C*(h))
│   ├── paper0_ecob_benchmark.m        Section 5.1: ECOB-proxy benchmark on case studies +
│   │                                  N=100 sample, Figure 6 + Table 5 (Spearman rho)
│   └── paper0_sensitivity.m           Figure 7 (f_frag exponent sweep, class invariance)
├── data/                              Celestrak GP catalogue snapshot, May 2026
│   ├── celestrak_active.csv           14 727 active LEO satellites
│   └── celestrak_debris.csv           2 597 tracked debris (Fengyun-1C, Iridium-33,
│                                      Cosmos-2251, Cosmos-1408)
├── figures/                           created at runtime by the scripts
├── LICENSE                            MIT
└── README.md
```

## Reproducing the paper

Tested on MATLAB R2026a. Statistics and Machine Learning Toolbox is required only by `paper0_validation.m` and `paper0_ecob_benchmark.m` (uses `lognrnd` and `corr`); the core algorithm runs on base MATLAB. From the `code/` directory:

| Run                      | Reproduces                                                     | Console output                  |
|--------------------------|----------------------------------------------------------------|---------------------------------|
| `paper0_casestudies`     | Table 2 (case-study individual + collective indices)           | Numeric table                   |
| `paper0_figures`         | Figure 2 + Figure 3 (altitude distribution, asymmetry bars)    | PDFs to `../figures/`           |
| `paper0_validation`      | Section 4.2: Figure 4 + Table 3 (top-10 asymmetry, N=100)      | Stats per tier + LaTeX snippet  |
| `paper0_parametric`      | Section 4.3: Figure 5 + Table 4 (cost thresholds C*(h))        | PDF + threshold table           |
| `paper0_ecob_benchmark`  | Section 5.1: Figure 6 + Table 5 (Spearman rho, ECOB benchmark) | Rank correlation + LaTeX snippet|
| `paper0_sensitivity`     | Figure 7 (f_frag exponent sweep)                               | PDF + class invariance check    |

All paths are auto-discovered from `mfilename('fullpath')`. No manual configuration needed; just `cd` to `code/` and run.

Expected output for `paper0_casestudies`:

```
ENVISAT           | 1018.16  class 5 (Very High) |  903.19  class 5 (Very High)
Sentinel-6        |    0.02  class 2 (Low)       |    0.05  class 2 (Low)
Starlink V2 Mini  |    0.09  class 2 (Low)       |   29.95  class 5 (Very High)
```

Expected output for `paper0_validation` (N=100, `rng(20260605, 'twister')`):

```
Asymmetry R = C_col / C_ind
  median   = 73.41
  P(R > 1) = 100%
  p99      = 384.64
Top-10 all in 465-546 km (constellation band, Starlink-dominated)
```

Expected output for `paper0_ecob_benchmark` (same seed):

```
ENVISAT           C_col=  903.19  C_eco=  275.01  ratio= 0.30
Sentinel-6        C_col=    0.05  C_eco=    0.07  ratio= 1.56
Starlink V2 Mini  C_col=   29.95  C_eco=   96.39  ratio= 3.22
Spearman rho (joint n=103) = 0.906   Pearson on log10 = 0.943
```

## Data source

The two CSV files in `data/` are snapshots of the Celestrak GP (General Perturbations) catalogue downloaded from `https://celestrak.org/NORAD/elements/` on 29 May 2026, restricted to active LEO satellites and to the four principal tracked fragmentation clouds (Fengyun-1C 2007, Iridium-33 / Cosmos-2251 2009, Cosmos-1408 2021).

## Notes on the ECOB proxy

`ecob_proxy.m` is a **snapshot-only** proxy of the ECOB family of orbital-sustainability indices (Letizia, Colombo, Lewis and Krag). It applies the same kinetic-flux density as `collective_probability.m` but with three ECOB-aligned modifications: mass-based fragmentation weight `(M/M_ref)^0.75` (NASA Standard Breakup Model), fixed 200-yr horizon, no debris-evolution propagation. **It is not a re-implementation of the full ECOB**, which requires long-term Monte Carlo evolution via MASTER, MOCAT-MC or equivalent. Its only purpose is the internal benchmarking in Section 5.1 of the paper, where it produces a Spearman rank correlation of 0.91 with `collective_probability.m` on a joint sample of 103 catalogued LEO objects.

## Citation

If you use this code or data, please cite the paper and this software release:

> Toson, F. (2026). *SusLifePath Paper 0: simplified collision-risk classifier for LEO satellites — code and data* (v2.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.20554998

## Licence

MIT (see `LICENSE`). You may use, modify and redistribute the code freely, including for commercial purposes, provided the copyright notice is retained.

## Contact

Federico Toson — `federico.toson@unipd.it`
CISAS "G. Colombo", University of Padova
