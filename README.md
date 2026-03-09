# FBL Wiretap Channel — Average Information Leakage (AIL) Analysis

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![IEEE Journal](https://img.shields.io/badge/IEEE-Transactions%20on%20Wireless%20Communications-red.svg)](https://ieeexplore.ieee.org)

> **Paper:** "Performance Analysis of Finite Blocklength Transmissions Over Wiretap Fading Channels: An Average Information Leakage Perspective"  
> **Authors:** Milad Tatar Mamaghani, Xiangyun Zhou, Nan Yang, A. Lee Swindlehurst, H. Vincent Poor  
> **Affiliations:** Australian National University · UC Irvine · Princeton University  
> **Venue:** IEEE Transactions on Wireless Communications (Under Review / Accepted)

---

## 📖 Research Summary

This repository contains the complete MATLAB simulation code for our study on **physical-layer security (PLS)** of **finite blocklength (FBL)** wireless transmissions over wiretap fading channels.

Traditional PLS relies on the concept of *secrecy capacity* achieved with infinite-length codewords. However, next-generation wireless systems — particularly **uRLLC** and **mMTC** in beyond-5G networks — require short FBL packets, which introduce a fundamental tension between **reliability** and **secrecy**.

We address this using a new metric: the **Average Information Leakage (AIL)**, denoted δ̄. Key aspects:

- We derive **exact and approximate closed-form AIL expressions** for general fading channels with Gaussian signaling
- Case studies use **Artificial Noise (AN) beamforming** in a MISO wiretap channel under **Rayleigh** and **Rician fading**
- We uncover a clean analytical relationship between AIL (FBL regime) and the classical **Secrecy Outage Probability** (SOP, IBL regime)
- We formulate and solve an **Average Secrecy Throughput (AST) maximization** via both non-adaptive (statistical CSI) and adaptive (instantaneous CSI) strategies

### Key Findings

- Allowing a *small* amount of AIL can yield **significant reliability improvements**
- The **Rician LoS component** reduces AIL compared to Rayleigh fading
- AN beamforming achieves **asymptotically zero AIL** at high SNR; MRT does not
- Both **blocklength design** and **AN power allocation** critically impact AST

---

## 🏗️ Repository Structure

```
FBL-WiretapChannel-AIL/
│
├── README.md                   ← This file
├── LICENSE                     ← MIT License
├── .gitignore
│
├── src/                        ← Core reusable MATLAB functions
│   ├── beamforming.m           ← AN / MRT beamforming matrix generation
│   ├── channel_dispersion.m    ← Channel dispersion V(γ) computation
│   ├── eval_point.m            ← Evaluation point x₀ calculation (Laplace approx.)
│   ├── ail_approx.m            ← Approximate AIL (Proposition 1, all cases)
│   ├── ail_exact.m             ← Exact AIL via numerical integration
│   ├── ail_simulate.m          ← Monte-Carlo AIL simulation
│   ├── dist_estimate.m         ← Gamma / Rician parameter estimation for γ_e
│   ├── generalized_marcum_q.m  ← Generalized Marcum Q-function
│   ├── gamma_hb_pdf.m          ← PDF of main channel power gain
│   └── smooth_surf.m           ← 2D Gaussian surface smoother (for AST plots)
│
├── config/
│   └── sys_params.m            ← Centralised system parameter configuration
│
├── simulations/                ← Experiment scripts (reproduce paper figures)
│   ├── run_all.m               ← Master runner: launches all sims + plots
│   ├── sim_rayleigh_fading.m   ← AIL sims over Rayleigh fading (Figs. 3–4)
│   ├── sim_rician_fading.m     ← AIL sims over Rician fading  (Figs. 3–4)
│   ├── fig_cdf_comparison.m    ← CDF validation (Fig. 2)
│   ├── fig_ail_vs_epsilon_alpha.m          ← AIL vs ε and α  (Figs. 5–6)
│   ├── fig_ail_vs_blocklength_kfactor.m    ← AIL vs N, K, k  (Figs. 7–8)
│   ├── fig_ail_sop_comparison.m            ← AIL vs SOP      (Fig. 9)
│   ├── fig_kfactor_impact.m                ← K-factor impact (Fig. 10)
│   ├── fig_ast_optimization_surface.m      ← AST surface     (Fig. 11)
│   └── optimization/
│       ├── optim_non_adaptive.m    ← Non-adaptive AST maximisation (Sec. V-A)
│       ├── optim_adaptive.m        ← Adaptive AST maximisation (Sec. V-B)
│       └── optim_adaptive_ao.m     ← Closed-form AO variant
│
├── data/
│   ├── raw/                    ← (empty) Placeholder for raw channel data
│   └── processed/              ← Generated .mat result files (see note below)
│       ├── FinalResults.mat         ← ⚠️  MISSING — generate with run_all.m
│       └── Results_optimalAST.mat   ← ⚠️  MISSING — generate with fig_ast_optimization_surface.m
│
├── results/
│   └── figures/                ← Auto-saved PDF/PNG figures
│
├── docs/
│   ├── system_model.md         ← System model description
│   └── missing_mat_files.md    ← Guidance on regenerating .mat files
│
└── examples/
    └── quick_start.m           ← Minimal runnable demo
```

---

## ⚙️ System Model

Alice (k-antenna AP) → Bob (single-antenna IoT device) in the presence of passive Eve.

```
Alice  ──[main link, gb]──►  Bob
  │
  └──[wiretap link, ge]──►  Eve (passive, CSI unknown to Alice/Bob)
```

**Beamforming:** Alice uses either:
- **AN beamforming** (power split α to info signal, 1-α to artificial noise in null-space of hb)
- **MRT beamforming** (α = 1, all power to info signal)

**FBL secrecy rate** (normal approximation):
```
R*_s ≈ [C_s(γ_b, γ_e) − √(V_b/N)·Q⁻¹(ε) − √(V_e/N)·Q⁻¹(δ)]+
```

**AIL (Proposition 1 — Laplace approximation):**
```
δ̄ ≈ 1 − F_{γ_e}(x₀),    where  x₀ = (1+γ_b)/2^{R₀} − 1
```

---

## 🔧 Installation & Dependencies

### Requirements

| Requirement | Version | Purpose |
|---|---|---|
| MATLAB | R2021a or later | Core simulation engine |
| Optimization Toolbox | Any | `fmincon`, `solve`, `optimvar` |
| Statistics and Machine Learning Toolbox | Any | `fitdist`, `makedist`, `cdf` |
| Global Optimization Toolbox | Any | `ga` (for adaptive scheme) |
| Symbolic Math Toolbox | Any | `sym`, `vpasolve`, `int` |

### Setup

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/FBL-WiretapChannel-AIL.git
cd FBL-WiretapChannel-AIL
```

Then in MATLAB:
```matlab
% Add all source folders to path
addpath(genpath('.'));
```

---

## 🚀 How to Run

### Option 1: Run Everything (Full Reproduction)

```matlab
cd simulations
run_all    % Runs Rayleigh + Rician sims, saves results, plots all figures
```

> ⏱️ **Estimated runtime:** 2–8 hours depending on hardware (Monte-Carlo loops with `simNum=1e4`).

### Option 2: Quick Demo

```matlab
cd examples
quick_start   % Minimal demo: computes AIL for a single parameter set
```

### Option 3: Reproduce Individual Figures

| Figure | Script | Description |
|---|---|---|
| Fig. 2 | `simulations/fig_cdf_comparison.m` | CDF validation of γ_e approximation |
| Figs. 3–4 | `simulations/run_all.m` | AIL vs SNR (Rayleigh + Rician) |
| Figs. 5–6 | `simulations/fig_ail_vs_epsilon_alpha.m` | AIL vs ε and α |
| Figs. 7–8 | `simulations/fig_ail_vs_blocklength_kfactor.m` | AIL vs N, K, k |
| Fig. 9 | `simulations/fig_ail_sop_comparison.m` | AIL vs SOP |
| Fig. 10 | `simulations/fig_kfactor_impact.m` | K-factor impact |
| Figs. 11–12 | `simulations/fig_ast_optimization_surface.m` | AST optimisation surface |

---

## ⚠️ Missing `.mat` Files

Two pre-computed result files are **not included** in this repository due to size, but they are regenerated automatically by the simulation scripts:

### `data/processed/FinalResults.mat`

**Used by:** `simulations/run_all.m` (plotting section), and any post-processing scripts.

**Contains:**
```
deltaBar_Rayl_MRT            (3 × 30)   Exact AIL, Rayleigh MRT
deltaBar_Rayl_MRT_Approx     (3 × 30)   Approx AIL, Rayleigh MRT
deltaBar_Rayl_MRT_Approx_HighSNR (3×30) High-SNR approx, Rayleigh MRT
deltaBar_Rayl_ANI            (3 × 30)   Exact AIL, Rayleigh AN beamforming
deltaBar_Rayl_ANI_Approx     (3 × 30)   Approx AIL, Rayleigh AN
deltaBar_Rayl_ANI_Approx_HighSNR (3×30) High-SNR approx, Rayleigh AN
deltaBar_Rice_MRT            (3 × 30)   Exact AIL, Rician MRT
deltaBar_Rice_MRT_APPROX     (3 × 30)   Approx AIL, Rician MRT
deltaBar_Rice_ANI            (3 × 30)   Exact AIL, Rician AN beamforming
deltaBar_Rice_ANI_Approx     (3 × 30)   Approx AIL, Rician AN
... (and SIM counterparts)
```

**To regenerate:**
```matlab
% In MATLAB, from repo root:
addpath(genpath('.'));
cd simulations
NVec = [5e2, 5e3, 5e4];
SNRVec_dB = linspace(-20, 10, 30);
sim_rayleigh_fading   % saves to data/processed/FinalResults.mat
sim_rician_fading     % appends to data/processed/FinalResults.mat
```

---

### `data/processed/Results_optimalAST.mat`

**Used by:** `simulations/fig_ast_optimization_surface.m`

**Contains:**
```
AST_Rayl     (20 × 100)   Average Secrecy Throughput surface (α × N grid)
```

**To regenerate:**
```matlab
addpath(genpath('.'));
cd simulations
fig_ast_optimization_surface   % saves to data/processed/Results_optimalAST.mat
```

See `docs/missing_mat_files.md` for full technical details.

---

## 🔑 Key Parameters (default values in `config/sys_params.m`)

| Parameter | Symbol | Default | Description |
|---|---|---|---|
| `k` | k | 4 | Number of Alice's transmit antennas |
| `epsilon` | ε | 1e-3 | Target decoding error probability |
| `M` | m | 100 | Number of information bits per packet |
| `N` | N | 400 | Coding blocklength (channel uses) |
| `alpha` | α | 0.7 | AN power allocation fraction |
| `mu_B` | β_b | 3 | Alice-to-Bob large-scale channel gain |
| `mu_E` | β_e | 1 | Alice-to-Eve large-scale channel gain |
| `K_factor` | K | 5 | Rician K-factor (both links) |
| `phi` | φ | 1e-4 | AIL constraint for optimisation |
| `simNum` | — | 1e4 | Monte-Carlo channel realisations |

---

## 📈 Results Overview

Representative outputs from the paper:

- **AIL vs SNR**: Both exact and approximate expressions closely match Monte-Carlo simulation. AN beamforming achieves zero AIL at high SNR; MRT does not.
- **AIL vs ε**: Relaxing reliability (larger ε) monotonically reduces AIL — a fundamental reliability–secrecy trade-off.
- **AIL vs N**: Longer blocklengths always improve AIL; the rate of improvement depends on fading type and K-factor.
- **AST Optimisation**: A unique global maximum exists in the (α, N) space; the optimal N scales approximately linearly with m.

---

## 📚 Citation

If you use this code, please cite our paper:

```bibtex
@article{tatarmamaghani2024fbl,
  title   = {Performance Analysis of Finite Blocklength Transmissions Over
             Wiretap Fading Channels: An Average Information Leakage Perspective},
  author  = {Tatar Mamaghani, Milad and Zhou, Xiangyun and Yang, Nan and
             Swindlehurst, A. Lee and Poor, H. Vincent},
  journal = {IEEE Transactions on Wireless Communications},
  year    = {2024},
  note    = {Preliminary version presented at IEEE ICC 2024, Denver, CO, USA}
}
```

**Preliminary conference version (IEEE ICC 2024):**
```bibtex
@inproceedings{tatarmamaghani2024fbl_icc,
  title     = {Average Information Leakage of Finite Blocklength Transmissions
               Over Wiretap Channels},
  author    = {Tatar Mamaghani, Milad and Zhou, Xiangyun and Yang, Nan and
               Swindlehurst, A. Lee and Poor, H. Vincent},
  booktitle = {IEEE International Conference on Communications (ICC)},
  year      = {2024},
  address   = {Denver, CO, USA}
}
```

**Supported by:**
- Australian Research Council Discovery Projects (DP220101318)
- US National Science Foundation (CNS-2128448, ECCS-2335876)

---

## 📬 Contact

| Author | Affiliation | Email |
|---|---|---|
| Milad Tatar Mamaghani | Australian National University | milad.tatarmamaghani@anu.edu.au |
| Xiangyun Zhou | Australian National University | xiangyun.zhou@anu.edu.au |
| Nan Yang | Australian National University | nan.yang@anu.edu.au |
| A. Lee Swindlehurst | University of California, Irvine | swindle@uci.edu |
| H. Vincent Poor | Princeton University | poor@princeton.edu |

---

## 📄 License

This project is released under the [MIT License](LICENSE). See `LICENSE` for details.

> **Note:** This code is provided for academic and research purposes. Please cite the paper if you use it.
