"""
Validation: FeSe/STO Tc Enhancement
====================================
Reproduces the blind prediction of Tc = 65.4 K for monolayer FeSe on SrTiO₃.

Protocol:
  1. Invert λ_sf from the FeSe/SVO/STO buffer baseline (Tc = 48 K)
  2. Add the STO interface phonon coupling (λ_ph = 0.19)
  3. Predict Tc without re-fitting any parameters

All energy scales from published ARPES/RIXS/DFT data.
"""
import numpy as np
import sys
sys.path.insert(0, '.')
from bridge_tc import bridge_tc, invert_lambda_sf, isotope_exponent

# ── Material Parameters ──────────────────────────────────────────────
lam_ph_bulk = 0.17          # FeSe intrinsic e-ph coupling (DFT)
lam_ph_interface = 0.19     # STO Fuchs-Kliewer phonon coupling (ARPES)
lam_ph_total = lam_ph_bulk + lam_ph_interface

omega_ph_bulk = 230.0       # K, FeSe Debye temperature
omega_STO_FK = 1160.0       # K, STO FK phonon (100 meV)
omega_sf_eff = 290.0        # K, effective spin-fluctuation frequency (25 meV)

E_F = 15.0                  # meV, Fermi energy (Coldea et al.)
gamma = 25.0 / E_F          # Migdal vertex ratio = 1.667

# Effective phonon frequency (geometric weighted mean)
w_b = lam_ph_bulk / lam_ph_total
w_i = lam_ph_interface / lam_ph_total
omega_ph_eff = np.exp(w_b * np.log(omega_ph_bulk) + w_i * np.log(omega_STO_FK))

# ── Execution ────────────────────────────────────────────────────────
print("=" * 65)
print("  FeSe/STO Tc Enhancement — Bridge Mechanism Validation")
print("=" * 65)

# Step 1: Invert λ_sf from buffer baseline
print(f"\n  Step 1: Calibrate λ_sf from FeSe/SVO/STO (Tc = 48 K)")
print(f"  Parameters: λ_ph = {lam_ph_bulk}, ω_ph = {omega_ph_bulk} K,")
print(f"              ω_sf = {omega_sf_eff} K (25 meV), γ = {gamma:.3f}")

lsf = invert_lambda_sf(lam_ph_bulk, omega_ph_bulk, omega_sf_eff, gamma, 48.0)
tc_cal, lt_cal, ol_cal, o2_cal = bridge_tc(
    lam_ph_bulk, lsf, omega_ph_bulk, omega_sf_eff, gamma
)

print(f"\n  Result: λ_sf = {lsf:.3f}")
print(f"          λ_total = {lt_cal:.3f}")
print(f"          Tc = {tc_cal:.1f} K (target: 48.0 K)")

# Step 2: Blind prediction — add interface phonons
print(f"\n  Step 2: BLIND PREDICTION — add STO interface phonons")
print(f"  Added: λ_ph(interface) = {lam_ph_interface}")
print(f"         ω_ph(effective) = {omega_ph_eff:.0f} K")
print(f"  λ_sf held fixed at {lsf:.3f} (NOT re-fitted)")

tc_pred, lt_pred, ol_pred, o2_pred = bridge_tc(
    lam_ph_total, lsf, omega_ph_eff, omega_sf_eff, gamma
)

err = abs(tc_pred - 65.0) / 65.0 * 100

print(f"\n  ┌─────────────────────────────────────────┐")
print(f"  │  Tc (predicted) = {tc_pred:>6.1f} K               │")
print(f"  │  Tc (measured)  = {65.0:>6.1f} K               │")
print(f"  │  Error          = {err:>6.1f} %               │")
print(f"  └─────────────────────────────────────────┘")

# Step 3: Cross-coupling contribution
tc_no_cross, _, _, _ = bridge_tc(lam_ph_total, lsf, omega_ph_eff, omega_sf_eff, 0.0)
cross_K = tc_pred - tc_no_cross
cross_pct = cross_K / (tc_pred - tc_cal) * 100

print(f"\n  Cross-coupling analysis:")
print(f"    Without γ (additive only): Tc = {tc_no_cross:.1f} K")
print(f"    Cross-coupling adds:       {cross_K:.1f} K")
print(f"    Share of 48→65 K jump:     {cross_pct:.0f}%")
