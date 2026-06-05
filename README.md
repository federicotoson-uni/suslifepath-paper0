# SusLifePath Paper 0 — Code and Data

Companion repository for the manuscript

> Toson, F. *A simplified engineering algorithm for collision risk assessment and classification of LEO satellites.* Acta Astronautica, submitted 2026.

The MATLAB toolchain reproduces all numerical results, figures and tables of the paper from a May 2026 Celestrak General Perturbations catalogue snapshot.

## Contents

```
suslifepath-paper0/
├── code/                              MATLAB toolchain
│   ├── databasecreator_real.m         build resident population from Celestrak CSV
│   ├── individual_probability_flux.m  P_ind = rho_all * sigma * v_rel * T_op
│   ├── collective_probability.m       P_col = rho_op * sigma * v_rel * T_tot * f_frag
│   ├── risk_index.m                   normalise + classify (5-class ECSS-compatible)
│   ├── paper0_casestudies.m           Table 2: three case-study indices
│   ├── paper0_figures.m               Figure 2 (altitude distribution) + Figure 3 (asymmetry bars)
│   ├── paper0_parametric.m            Figure 4 + Table 3 (cost thresholds)
│   └── paper0_sensitivity.m           Figure 5 (f_frag exponent sensitivity)
├── data/                              Celestrak GP catalogue snapshot, May 2026
│   ├── celestrak_active.csv           14 727 active LEO satellites
│   └── celestrak_debris.csv           2 597 tracked debris (Fengyun-1C, Iridium-33,
│                                      Cosmos-2251, Cosmos-1408)
├── figures/                           created at runtime by the scripts
├── LICENSE                            MIT
└── README.md
```

## Reproducing the paper

Tested on MATLAB R2026a. No toolboxes beyond base MATLAB. From the `code/` directory:

| Run            | Reproduces                                  | Console output |
|----------------|---------------------------------------------|----------------|
| `paper0_casestudies` | Table 2 (case-study individual + collective indices) | Numeric table |
| `paper0_figures`     | Figure 2, Figure 3 + Table 2 sanity print  | PDFs to `../figures/` |
| `paper0_parametric`  | Figure 4 + Table 3 (cost thresholds C*(h)) | PDF + threshold table |
| `paper0_sensitivity` | Figure 5 (f_frag exponent sweep)           | PDF + class invariance check |

All paths are auto-discovered from `mfilename('fullpath')`. No manual configuration needed; just `cd` to `code/` and run.

Expected output for `paper0_casestudies`:

```
ENVISAT           | 1018.16  class 5 (Very High) |  903.19  class 5 (Very High)
Sentinel-6        |    0.02  class 2 (Low)       |    0.05  class 2 (Low)
Starlink V2 Mini  |    0.09  class 2 (Low)       |   29.95  class 5 (Very High)
```

## Data source

The two CSV files in `data/` are snapshots of the Celestrak GP (General Perturbations) catalogue downloaded from `https://celestrak.org/NORAD/elements/` on 29 May 2026, restricted to active LEO satellites and to the four principal tracked fragmentation clouds (Fengyun-1C 2007, Iridium-33 / Cosmos-2251 2009, Cosmos-1408 2021).

## Citation

If you use this code or data, please cite the paper and this software release. Once the Zenodo DOI is minted, the citation becomes:

> Toson, F. (2026). *SusLifePath Paper 0: simplified collision-risk classifier for LEO satellites — code and data.* Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

## Licence

MIT (see `LICENSE`). You may use, modify and redistribute the code freely, including for commercial purposes, provided the copyright notice is retained.

## Contact

Federico Toson — `federico.toson@unipd.it`
CISAS "G. Colombo", University of Padova
