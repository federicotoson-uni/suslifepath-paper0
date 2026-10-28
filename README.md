# SusLifePath Paper 0: code and data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20554997.svg)](https://doi.org/10.5281/zenodo.20554997)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Companion repository for the manuscript

> Toson, F. *A simplified engineering algorithm for collision risk assessment and classification of LEO satellites.* Advances in Space Research, 2026, under revision. Pre-print: [arXiv:2606.17947](https://doi.org/10.48550/arXiv.2606.17947).

The MATLAB toolchain reproduces all numerical results, figures and tables of the paper from a May 2026 Celestrak General Perturbations catalogue snapshot. Version 2.1.0 is the revision for Advances in Space Research: a catalogued mission is removed from its own shell count (`mission.in_catalogue = true`), the two indices are evaluated for all 14 727 active LEO objects in place of the N = 100 sample, and two analyses are added (Sentinel-6 end-of-life options, collision avoidance). The scripts of the submitted version are kept and marked as such.

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
│   ├── paper0_sentinel6_eol.m         Section 4.1: Sentinel-6 collective index under the screening
│   │                                  convention, a graveyard orbit and a 25-year descent
│   ├── paper0_parametric.m            Section 4.3: Figure 5 + Table 4 (cost thresholds C*(h))
│   ├── paper0_cubesat_corollary.m     Section 4.3: educational smallsat corollary
│   │                                  (3U CubeSat / 2P PocketCube at the Starlink shell)
│   ├── paper0_ecob_benchmark.m        Section 5.1: Table 5 (ECOB proxy on the three cases);
│   │                                  the N=100 benchmark of the submitted version
│   ├── paper0_literature_anchor.m     Section 5.1: Table 6, anchor against published
│   │                                  ECOB values for MetOp-A and Sentinel-2
│   ├── paper0_sensitivity.m           Figure 7 (f_frag exponent sweep, class invariance)
│   ├── paper0_catalogue.m             Sections 4.2 and 5.1: all active LEO objects. Table 3,
│   │                                  Figure 4, Figure 6, Spearman rho against the ECOB proxy,
│   │                                  sensitivities A2 (target fragments), A4 (eccentric
│   │                                  orbits), A5 (collision avoidance), past-2025 snapshot
│   ├── paper0_validation.m            submitted version: N=100 random sample (superseded by
│   │                                  paper0_catalogue.m)
│   ├── paper0_seed_sweep.m            submitted version: 100-seed sweep on the N=100 sample
│   └── paper0_robustness.m            submitted version: temporal robustness on the N=100 sample
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

Tested on MATLAB R2026a. Statistics and Machine Learning Toolbox is required by `paper0_catalogue.m`, `paper0_validation.m` and `paper0_ecob_benchmark.m` (`lognrnd`, `corr`, `prctile`); the core algorithm runs on base MATLAB. From the `code/` directory:

| Run                        | Reproduces                                                     | Console output                  |
|----------------------------|----------------------------------------------------------------|---------------------------------|
| `paper0_casestudies`       | Table 2 (case-study individual + collective indices)           | Numeric table                   |
| `paper0_figures`           | Figure 2 + Figure 3 (altitude distribution, asymmetry bars)    | PDFs to `../figures/`           |
| `paper0_sentinel6_eol`     | Section 4.1: Sentinel-6 end-of-life options (convention, graveyard, 25-yr descent with and without avoidance) | Normalised collective index and class per option |
| `paper0_parametric`        | Section 4.3: Figure 5 + Table 4 (cost thresholds C*(h))        | PDF + threshold table           |
| `paper0_cubesat_corollary` | Section 4.3: educational smallsat corollary (R~43 at 1 MUSD, R~458 at 50 kUSD) | Numeric table + C* per geometry |
| `paper0_ecob_benchmark`    | Section 5.1: Table 5 (ECOB proxy on the three cases); the N=100 benchmark of the submitted version | Case table + rank correlation   |
| `paper0_literature_anchor` | Section 5.1: Table 6 (anchor vs published ECOB values)         | Ratios + LaTeX snippet          |
| `paper0_sensitivity`       | Figure 7 (f_frag exponent sweep)                               | PDF + class invariance check    |
| `paper0_catalogue`         | Sections 4.2 and 5.1: Table 3, Figure 4 (`Figure_6_validation.pdf`), Figure 6 (`Figure_7_ecob_benchmark.pdf`), Spearman rho on the whole catalogue, A2, A4, A5, past-2025 snapshot. Run it last: it writes the figure files used in the paper | Group statistics + sensitivities |
| `paper0_validation`        | Submitted version: N=100 sample statistics (`Figure_6_validation_submitted.pdf`) | Stats per tier + LaTeX snippet  |
| `paper0_seed_sweep`        | Submitted version: 100-seed sweep statistics                   | Sweep ranges + Clopper-Pearson  |
| `paper0_robustness`        | Submitted version: temporal robustness on past-2025 snapshot   | Side-by-side stats + LaTeX      |

All paths are auto-discovered from `mfilename('fullpath')`. No manual configuration needed; just `cd` to `code/` and run. The whole set runs in about ten minutes; `paper0_catalogue` alone takes under a minute.

Expected output for `paper0_casestudies` (Table 2 of the paper prints three decimals: 0.011, 0.023, 0.094, 29.9):

```
ENVISAT            | 1018.16  class 5 (Very High) | 903.19  class 5 (Very High)
Sentinel-6         |  0.01  class 2 (Low      ) |  0.02  class 2 (Low      )
Starlink V2 Mini   |  0.09  class 2 (Low      ) | 29.93  class 5 (Very High)
```

Expected output for `paper0_sentinel6_eol`:

```
Sentinel-6, normalised collective index (operational objects in the shell without 6A: 1)
convention: 25 yr at the 1336 km density               0.023  L
graveyard in the band, 200-yr horizon                  0.157  M
25-yr descent, perigee 350 km (density x869)          16.594  VH
25-yr descent, perigee 450 km (density x1441)         27.494  VH
25-yr descent, perigee 550 km (density x510)           9.745  H
descent with manoeuvres 90%, perigee 350 km (x241)     4.608  H
descent with manoeuvres 90%, perigee 450 km (x349)     6.656  H
descent with manoeuvres 90%, perigee 550 km (x174)     3.315  H
descent with manoeuvres 99%, perigee 350 km (x178)     3.410  H
descent with manoeuvres 99%, perigee 450 km (x239)     4.573  H
descent with manoeuvres 99%, perigee 550 km (x140)     2.672  H
Reference disposed on 700 x 550 km: density x10.7, reference collective x8.6, Starlink R 37
Reference disposed on 700 x 600 km: density x3.0, reference collective x2.6, Starlink R 124
```

Expected output for `paper0_catalogue` (D = 100 draws of the priors, `rng(20260605, 'twister')`; brackets are the minimum and maximum over the draws):

```
group                           n | median R [min max]     | share R>1 [min max]    | share R>100 [min max]
all                         14719 |   84.1 [  82.9   85.5] | 0.994 [0.992 0.995] | 0.381 [0.370 0.389]
constellations (by name)    11684 |   88.0 [  86.8   89.8] | 1.000 [1.000 1.000] | 0.408 [0.395 0.418]
others, below 600 km         2124 |   87.0 [  84.6   90.2] | 1.000 [1.000 1.000] | 0.398 [0.373 0.429]
others, 600-1500 km           908 |    4.4 [   4.1    4.9] | 0.900 [0.870 0.925] | 0.000 [0.000 0.002]
others, above 1500 km           3 |    3.4 [   1.5    7.3] | 1.000 [0.667 1.000] | 0.000 [0.000 0.000]
Starlink V2 Mini  ind    0.0943  col   29.9328  R   317.30
525 km shell: 2389 operational, 1821 Starlink (76%); against the other operators R 75, collective 7.1
Reference at 550 km: operational  96% of the shell, Starlink R    62, catalogue median   16.4
Reference at 700 km: operational  19% of the shell, Starlink R   317, catalogue median   84.1
Spearman R_col vs ECOB proxy (n = 14722): median 0.944 [0.943 0.947]
A2: Spearman with R_col 0.9990 [0.9989 0.9991]; median R 64.0 (was 84.1); share R>1 0.990
A4: 120 objects with apogee - perigee > 100 km
    operational density, orbit-averaged/perigee: median 0.60 [0.00 777.88]; effect on R: median 0.96 [0.17 4.14]
A5 eta 0.90, manoeuvrable targets by name  : Starlink R    81.1; constellation median R    17.3, share R>100 0.005; all, share R>1 0.993
A5 eta 0.99, manoeuvrable targets by name  : Starlink R    50.9; constellation median R     7.7, share R>100 0.003; all, share R>1 0.988
Snapshot without the 2026 launches: 12962 objects; median R 87.0 [85.7 88.0]; share R>1 0.993
```

Expected output for `paper0_ecob_benchmark` (Table 5; the Spearman value refers to the N=100 joint sample of the submitted version):

```
Case                 M [kg]        C_ind        C_col        C_eco
ENVISAT                8211     1018.158      903.186      275.013
Sentinel-6             1192        0.011        0.023        0.036
Starlink V2 Mini        800        0.094       29.933       96.355
  Spearman rho = 0.9058   (joint sample n=103)
```

Expected output for `paper0_validation` (submitted version, N=100, `rng(20260605, 'twister')`; its case-study block prints the submitted values, Sentinel-6 0.022 and 0.047, without the self-count correction):

```
Asymmetry R = C_col / C_ind
  median   = 73.41
  P(R > 1) = 100%
  p99      = 384.64
Top-10 all in 465-546 km (constellation band, Starlink-dominated)
```

Expected output for `paper0_seed_sweep` (submitted version, 100 seeds, ~5 min runtime):

```
Statistic              median        p05        p95         SD
median R                77.86      67.70      90.93       6.47
P(R>1)                  0.990      0.970      1.000      0.010
Spearman rho            0.915      0.885      0.938      0.017

Clopper-Pearson 95% lower bound for 100/100: P(R>1) >= 96.38%
```

Expected output for `paper0_robustness` (submitted version, May 2026 vs synthetic past 2025):

```
                      May 2026     Past 2025
Active LEO objects    14 727       12 962
Median R              73.4         74.1
P(R>1)                100%         99%
Spearman rho          0.906        0.922
```

Expected output for `paper0_literature_anchor` (Table 6; the proxy column of the table is `C_eco` divided by its ENVISAT value, 0.144 and 0.025):

```
Mission             C_col        C_eco Letizia ECOB   Ours/ENVISAT Letizia/ENVISAT
ENVISAT           903.186      275.013        1.000          1.000          1.000
MetOp-A            20.747       39.693        0.180          0.023          0.180
Sentinel-2          4.729        6.819        0.030          0.005          0.030
```

## Data source

The CSV files in `data/` are snapshots of the Celestrak GP (General Perturbations) catalogue downloaded from `https://celestrak.org/NORAD/elements/`:
- `celestrak_active.csv` and `celestrak_debris.csv`: 29 May 2026, restricted to active LEO satellites and to the four principal tracked fragmentation clouds (Fengyun-1C 2007, Iridium-33 / Cosmos-2251 2009, Cosmos-1408 2021). Rocket bodies and defunct payloads are not part of the snapshot; the paper says so in Section 5.2.
- `celestrak_active_past2025.csv`: derived from the live June 2026 active-satellite catalogue by removing all objects launched in 2026 (the filter `OBJECT_ID NOT LIKE '2026-%'` removes about 1 900 predominantly Starlink V2 Mini elements). This synthetic "past-2025" snapshot is used by `paper0_catalogue.m` and `paper0_robustness.m` for the temporal-robustness check of Section 4.2.

## Notes on the ECOB proxy

`ecob_proxy.m` is a **snapshot-only** proxy of the ECOB family of orbital-sustainability indices (Letizia, Colombo, Lewis and Krag). It applies the same kinetic-flux density as `collective_probability.m` but with three ECOB-aligned modifications: mass-based fragmentation weight `(M/M_ref)^0.75` (NASA Standard Breakup Model), fixed 200-yr horizon, no debris-evolution propagation. **It is not a re-implementation of the full ECOB**, which requires long-term Monte Carlo evolution via MASTER, MOCAT-MC or equivalent. Its only purpose is the internal sensitivity check in Section 5.1 of the paper, where it produces a Spearman rank correlation of 0.94 with `collective_probability.m` on the whole catalogue (0.91 on the joint sample of 103 objects of the submitted version).

## AI disclosure

Consistent with Elsevier policy, the manuscript that this repository accompanies includes a "Declaration of generative AI" section disclosing that an AI assistant was used during manuscript preparation for language refinement, computational scripting, and statistical sensitivity analysis. All scientific content, numerical results, and methodological choices are the responsibility of the sole author, who verified each computation independently against the open-source toolchain in this repository.

## Citation

If you use this code or data, please cite the paper and this software release:

> Toson, F. (2026). *suslifepath-paper0: collision-risk index toolchain and catalogue snapshot* (v2.1.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.20554997
>
> The concept DOI [10.5281/zenodo.20554997](https://doi.org/10.5281/zenodo.20554997) always resolves to the latest version; each release also carries its own version DOI on Zenodo.

## Licence

MIT (see `LICENSE`). You may use, modify and redistribute the code freely, including for commercial purposes, provided the copyright notice is retained.

## Contact

Federico Toson, `federico.toson@unipd.it`
CISAS "G. Colombo", University of Padova
