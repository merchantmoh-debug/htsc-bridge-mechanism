# Cross-Channel Vertex Corrections in Multi-Boson Superconductors

**A Coupling-Level Bridge Mechanism for High-Temperature Superconductivity**

[![arXiv](https://img.shields.io/badge/arXiv-2026.xxxxx-b31b1b.svg)](https://arxiv.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Abstract

We derive a coupling-level correction to the electron-boson interaction in superconductors with coexisting phonon and spin-fluctuation pairing channels. The total coupling constant acquires a cross-channel term:

```
λ_total = λ_ph + λ_sf + γ · λ_ph · λ_sf
```

where γ = ω_sf / E_F is the first-order Migdal vertex correction, computed entirely from measurable material properties. This term is fed into the **unmodified** Allen-Dynes equation (PRB 12, 905, 1975) — no changes to the Tc solver.

We validate against three independent systems:

| System | Prediction | Measured | Error |
|--------|-----------|----------|-------|
| FeSe/SrTiO₃ Tc | 65.4 K | 65 K | 0.6% |
| YBa₂Cu₃O₇ isotope α_O | 0.021 | 0.017–0.025 | Within range |
| 2Δ/kTc (FeSe) | 5.31 | 5.36 | 0.9% |

We predict the oxygen isotope exponent for pressurized La₃Ni₂O₇: **α_O = 0.053** (ΔTc = 0.50 K), testable via high-pressure diamond anvil cell experiments.

## Quick Start

```bash
cd code/
python bridge_tc.py              # Self-test
python validate_fese.py          # FeSe/STO prediction
python validate_ybco.py          # YBCO isotope exponent
python validate_lno.py           # La₃Ni₂O₇ prediction
python validate_gap.py           # ω_sf convergence check
```

Requires only `numpy` (no exotic dependencies).

## Repository Structure

```
├── code/
│   ├── bridge_tc.py          # Core Tc calculator (the formula + solver)
│   ├── validate_fese.py      # FeSe/STO Tc enhancement
│   ├── validate_ybco.py      # YBCO oxygen isotope exponent
│   ├── validate_lno.py       # La₃Ni₂O₇ self-consistency + prediction
│   └── validate_gap.py       # Independent ω_sf convergence check
├── paper/
│   └── bridge_mechanism.tex  # Manuscript (LaTeX)
├── data/                     # Raw computational data (DFPT, when available)
└── README.md
```

## Key Results

### The Physics

The cross-coupling term `γ · λ_ph · λ_sf` is the first-order vertex correction in the Migdal expansion. When E_F is large (>200 meV), γ is negligible and the channels are effectively additive. When E_F drops to ~15 meV (as in monolayer FeSe on SrTiO₃), γ jumps to ~1.7 and the cross-term becomes the **dominant contribution** to the Tc enhancement.

This explains why the FeSe/STO interface shows an anomalously large Tc jump (8 K → 65 K) that neither phonon-only nor spin-only theories can account for.

### Universality

| System | E_F (meV) | γ | Cross-coupling |
|--------|-----------|---|----------------|
| FeSe/STO | 15 | 1.67 | **55% of Tc jump** |
| FeSe bulk | 30 | 0.13 | 1.5% |
| YBCO | 250 | 0.17 | 1.1% |
| La₃Ni₂O₇ | 500 | 0.14 | 2.0% |
| H₃S | 3000 | ~0 | 0% (BCS limit) |

The mechanism is universal but only experimentally dominant in the BCS-BEC crossover regime.

## Citation

```bibtex
@article{alzawahreh2026bridge,
  title={Cross-Channel Vertex Corrections in Multi-Boson Superconductors: 
         A Coupling-Level Bridge Mechanism for High-Temperature Superconductivity},
  author={Al-Zawahreh, Mohamad},
  year={2026},
  journal={arXiv preprint}
}
```

## License

Code: MIT License  
Paper: CC-BY 4.0
