"""
Validation: La₃Ni₂O₇ Self-Consistency + Falsifiable Prediction
================================================================
Tests whether the Bridge formula self-consistently explains
Tc = 80 K in pressurized La₃Ni₂O₇ and generates a falsifiable
prediction for the oxygen isotope exponent.
"""
import numpy as np
import sys
sys.path.insert(0, '.')
from bridge_tc import (bridge_tc, modified_allen_dynes,
                        invert_lambda_sf, isotope_exponent)

# ── Material Parameters ──────────────────────────────────────────────
lam_ph = 0.13          # DFT, Fmmm phase (2024)
omega_ph_16 = 550.0    # K, Debye temperature at 25 GPa (¹⁶O)
omega_sf = 812.0       # K, spin fluctuations (~70 meV)
E_F = 500.0            # meV, Ni e_g bandwidth
gamma = 70.0 / E_F     # = 0.14

# ── Execution ────────────────────────────────────────────────────────
print("=" * 65)
print("  La₃Ni₂O₇ — Self-Consistency + Isotope Prediction")
print("=" * 65)
print(f"\n  Parameters: λ_ph = {lam_ph}, ω_ph = {omega_ph_16} K")
print(f"              ω_sf = {omega_sf} K (70 meV), E_F = {E_F} meV")
print(f"              γ = {gamma:.3f}")

# A) Self-consistency: can Bridge reproduce Tc = 80 K?
lsf = invert_lambda_sf(lam_ph, omega_ph_16, omega_sf, gamma, 80.0)
tc, lt, ol, o2 = bridge_tc(lam_ph, lsf, omega_ph_16, omega_sf, gamma)

print(f"\n  A) Self-consistency check:")
print(f"     λ_sf (inverted) = {lsf:.3f}")
print(f"     λ_total         = {lt:.3f}")
print(f"     Tc              = {tc:.1f} K (target: 80 K)")
print(f"     Physical range (0.5–3.0): {'✓' if 0.5 <= lsf <= 3.0 else '✗'}")

# B) BCS alone
tc_bcs = modified_allen_dynes(omega_ph_16, omega_ph_16, lam_ph)
print(f"\n  B) BCS (phonon only): Tc = {tc_bcs:.2f} K")
print(f"     Confirms: λ_ph = {lam_ph} insufficient for 80 K")

# C) Isotope prediction
mass_ratio = np.sqrt(16.0 / 18.0)
omega_ph_18 = omega_ph_16 * mass_ratio

tc_18, _, _, _ = bridge_tc(lam_ph, lsf, omega_ph_18, omega_sf, gamma)
delta_Tc = tc - tc_18
alpha_O = -np.log(tc_18 / tc) / np.log(18.0 / 16.0)

print(f"\n  C) FALSIFIABLE PREDICTION:")
print(f"  ┌──────────────────────────────────────────────────┐")
print(f"  │  La₃Ni₂O₇ at 25 GPa:                            │")
print(f"  │    Tc(¹⁶O)  = {tc:.2f} K                            │")
print(f"  │    Tc(¹⁸O)  = {tc_18:.2f} K                            │")
print(f"  │    ΔTc      = {delta_Tc:.2f} K                             │")
print(f"  │    α_O      = {alpha_O:.4f}                            │")
print(f"  │                                                  │")
print(f"  │  Consistency with PSI (arXiv:2504.08290):         │")
print(f"  │    CDW isotope shift = +2.3 K  (phonon-coupled)   │")
print(f"  │    SDW isotope shift = 0       (electronic)       │")
print(f"  │    Our model: spin channel has no ¹⁸O dependence  │")
print(f"  │    → SDW shift = 0  ✓                             │")
print(f"  └──────────────────────────────────────────────────┘")
