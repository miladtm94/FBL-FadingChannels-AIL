# System Model

## Overview

This document describes the system model for the FBL wiretap channel study.

```
           k antennas
Alice (AP) ─────────────────────────────────► Bob (IoT device)
  │          Main channel: g_b = √β_b · h_b          │
  │                                                   │
  └─────────────────────────────────────────► Eve (passive eavesdropper)
             Wiretap channel: g_e = √β_e · h_e
```

## Channel Model

The received signals at Bob and Eve are:

```
y_r = √P · g_r · x + ξ_r,    r ∈ {b, e}
```

where:
- `P` — transmit power
- `g_r = √β_r · h_r` — complex channel vector (large-scale × small-scale)
- `β_r = β₀ · d_r^{-η}` — path loss with distance `d_r` and exponent `η`
- `h_r ∈ ℂ^{1×k}` — normalised small-scale fading vector
- `ξ_r ~ CN(0, σ²)` — AWGN
- `ρ = P/σ²` — transmit SNR

## Fading Models

### Rayleigh Fading
```
h_r ~ CN(0, I_k)     (i.i.d. entries)
```
- Bob's gain `γ̃_b = ρβ_b‖h_b‖²` ~ Gamma(k, ρβ_b)
- Eve's gain with MRT: `γ_e ~ Exp(ρβ_e)`
- Eve's gain with AN: has PDF given by eq. (19)

### Rician Fading
```
h_r = √(K/(1+K)) · ζ_r + √(1/(1+K)) · h̃_r
```
where `ζ_r` is the LoS component and `h̃_r ~ CN(0, I_k)`.

- Bob's gain `γ̃_b` ~ scaled non-central chi-squared (eq. 30–31)
- Eve's gain with MRT: scaled non-central chi-squared (eq. 39)
- Eve's gain with AN: approximated by Gamma(a_z, b_z) (eq. 34–35)

## Beamforming

### Artificial Noise (AN) Beamforming
Alice sends: `x = √α · w · s + √((1-α)/(k-1)) · U · v`

- `w` — principal eigenvector of `h_b^†h_b` (points towards Bob)
- `U` — null-space of `h_b` (orthonormal complement)
- `s` — unit-power information symbol
- `v ~ CN(0, I_{k-1})` — AN vector
- `α ∈ (0,1)` — power split factor

Result: Bob receives no interference from AN (`h_b U = 0`). Eve's channel is degraded.

### MRT Beamforming (α = 1)
Special case: `w = h_b^† / ‖h_b‖`, no AN. Maximises Bob's SNR.

## FBL Secrecy Rate

The achievable FBL secrecy coding rate (normal approximation):
```
R*_s ≈ [C_s(γ_b, γ_e) − √(V_b/N)·Q⁻¹(ε) − √(V_e/N)·Q⁻¹(δ)]+
```

where:
- `C_s = log₂(1+γ_b) − log₂(1+γ_e)` — IBL secrecy capacity
- `V_r = log₂²(e) · γ_r(γ_r+2)/(γ_r+1)²` — channel dispersion
- `N` — coding blocklength
- `ε` — decoding error probability at Bob
- `δ` — information leakage to Eve

## AIL Definition

Alice fixes ε and transmits at rate R_s = R*_s. The resulting leakage per realisation:
```
δ(γ_e) = Q(√(N/V_e) · [log₂((1+γ_b)/(1+γ_e)) − R₀])
```

The **Average Information Leakage (AIL)** averages over Eve's channel:
```
δ̄ = E{δ | h_b} = ∫₀^∞ Q(√(N/V_e(x)) · [log₂((1+γ_b)/(1+x)) − R₀]) · f_{γ_e}(x) dx
```

**Laplace approximation (Proposition 1):**
```
δ̄ ≈ 1 − F_{γ_e}(x₀),    x₀ = (1+γ_b)/2^{R₀} − 1
```

This equals the IBL Secrecy Outage Probability (SOP) with a specific redundancy rate — establishing the AIL–SOP equivalence (Corollary 1).
