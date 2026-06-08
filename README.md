# SusLifePath Paper 0 — Code and Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20554998.svg)](https://doi.org/10.5281/zenodo.20554998)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Companion repository for the manuscript

> Toson, F. *A simplified engineering algorithm for collision risk assessment and classification of LEO satellites.* Manuscript in preparation for submission to Acta Astronautica, 2026.

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
│   ├── paper0_seed_sweep.m            Section 4.2: 100-seed sweep — median R range,
│   │                                  Clopper-Pearson lower bound on P(R>1)
│   ├── paper0_robustness.m            Section 4.2: temporal robustness on past-2025 snapshot
│   ├── paper0_parametric.m            Section 4.3: Figure 5 + Table 4 (cost thresholds C*(h))
│   ├── paper0_ecob_benchmark.m        Section 5.1: ECOB-proxy benchmark on case studies +
│   │                                  N=100 sample, Figure 6 + Table 5 (Spearman rho)
│   ├── paper0_literature_anchor.m     Section 5.1: Table 6 — anchor against published
│   │                                  ECOB values for MetOp-A and Sentinel-2
│   └── paper0_sensitivity.m           Figure 7 (f_frag exponent sweep, class invariance)
├── data/                              Celestrak GP catalogue snapshots
│   ├── celestrak_active.csv           14 727 active LEO satellites, 29 May 2026 snapshot
│   ├── celestrak_active_past2025.csv  12 962 active LEO satellites, synthetic past-2025
│   │                                  snapshot (live June 2026 with 2026 launches removed)
│   └── celestrak_debris.csv           2 597 tracked debris (Fengyun-1C, Iridium-33,
│                                      Cosmos-2251, Cosmos-1408)
├── figures/                           created at runtime by the scripts
├── LICENSE                            MIT
└── README.md
```

## Reproducing the paper

Tested on MATLAB R2026a. Statistics and Machine Learning Toolbox is required only by `paper0_validation.m` and `paper0_ecob_benchmark.m` (uses `lognrnd` and `corr`); the core algorithm runs on base MATLAB. From the `code/` directory:

| Run                        | Reproduces                                                     | Console output                  |
|----------------------------|----------------------------------------------------------------|---------------------------------|
| `paper0_casestudies`       | Table 2 (case-study individual + collective indices)           | Numeric table                   |
| `paper0_figures`           | Figure 2 + Figure 3 (altitude distribution, asymmetry bars)    | PDFs to `../figures/`           |
| `paper0_validation`        | Section 4.2: Figure 4 + Table 3 (top-10 asymmetry, N=100)      | Stats per tier + LaTeX snippet  |
| `paper0_seed_sweep`        | Section 4.2: 100-seed sweep statistics                         | Sweep ranges + Clopper-Pearson  |
| `paper0_robustness`        | Section 4.2: temporal robustness on past-2025 snapshot         | Side-by-side stats + LaTeX      |
| `paper0_parametric`        | Section 4.3: Figure 5 + Table 4 (cost thresholds C*(h))        | PDF + threshold table           |
| `paper0_ecob_benchmark`    | Section 5.1: Figure 6 + Table 5 (Spearman rho, ECOB benchmark) | Rank correlation + LaTeX snippet|
| `paper0_literature_anchor` | Section 5.1: Table 6 (anchor vs published ECOB values)         | Ratios + LaTeX snippet          |
| `paper0_sensitivity`       | Figure 7 (f_frag exponent sweep)                               | PDF + class invariance check    |

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

Expected output for `paper0_seed_sweep` (100 seeds, ~5 min runtime):

```
Statistic              median        p05        p95         SD
median R                77.86      67.70      90.93       6.47
P(R>1)                  0.990      0.970      1.000      0.010
Spearman rho            0.915      0.885      0.938      0.017

Clopper-Pearson 95% lower bound for 100/100: P(R>1) >= 96.38%
```

Expected output for `paper0_robustness` (May 2026 vs synthetic past 2025):

```
                      May 2026     Past 2025
Active LEO objects    14 727       12 962
Median R              73.4         74.1
P(R>1)                100%         99%
Spearman rho          0.906        0.922
```

Expected output for `paper0_literature_anchor` (vs Letizia et al. 2017):

```
Mission        Ours/ENVISAT  Letizia/ENVISAT
ENVISAT        1.000         1.000
MetOp-A        0.144         0.180
Sentinel-2     0.025         0.030
(matches published ECOB values to within factor of ~2)
```

## Data source

The CSV files in `data/` are snapshots of the Celestrak GP (General Perturbations) catalogue downloaded from `https://celestrak.org/NORAD/elements/`:
- `celestrak_active.csv` and `celestrak_debris.csv`: 29 May 2026, restricted to active LEO satellites and to the four principal tracked fragmentation clouds (Fengyun-1C 2007, Iridium-33 / Cosmos-2251 2009, Cosmos-1408 2021).
- `celestrak_active_past2025.csv`: derived from the live June 2026 active-satellite catalogue by removing all objects launched in 2026 (the filter `OBJECT_ID NOT LIKE '2026-%'` removes about 1 900 predominantly Starlink V2 Mini elements). This synthetic "past-2025" snapshot is used by `paper0_robustness.m` for the temporal-robustness check of Section 4.2.

## Notes on the ECOB proxy

`ecob_proxy.m` is a **snapshot-only** proxy of the ECOB family of orbital-sustainability indices (Letizia, Colombo, Lewis and Krag). It applies the same kinetic-flux density as `collective_probability.m` but with three ECOB-aligned modifications: mass-based fragmentation weight `(M/M_ref)^0.75` (NASA Standard Breakup Model), fixed 200-yr horizon, no debris-evolution propagation. **It is not a re-implementation of the full ECOB**, which requires long-term Monte Carlo evolution via MASTER, MOCAT-MC or equivalent. Its only purpose is the internal benchmarking in Section 5.1 of the paper, where it produces a Spearman rank correlation of 0.91 with `collective_probability.m` on a joint sample of 103 catalogued LEO objects.

## AI disclosure

Consistent with Elsevier policy, the manuscript that this repository accompanies includes a "Declaration of generative AI" section disclosing that an AI assistant was used during manuscript preparation for language refinement, computational scripting, and statistical sensitivity analysis. All scientific content, numerical results, and methodological choices are the responsibility of the sole author, who verified each computation independently against the open-source toolchain in this repository.

## Citation

If you use this code or data, please cite the paper and this software release:

> Toson, F. (2026). *SusLifePath Paper 0: simplified collision-risk classifier for LEO satellites — code and data* (v2.0.0) [Software]. Zenodo.

## Licence

MIT (see `LICENSE`). You may use, modify and redistribute the code freely, including for commercial purposes, provided the copyright notice is retained.

## Contact

Federico Toson — `federico.toson@unipd.it`
CISAS "G. Colombo", University of Padova
