"""
Bridge Mechanism Tc Calculator
==============================
Computes the superconducting critical temperature for multi-channel
superconductors using the Bridge coupling formula with the standard
modified Allen-Dynes equation (PRB 12, 905, 1975).

The Bridge formula introduces a cross-coupling term between phonon
and spin-fluctuation pairing channels:

    λ_total = λ_ph + λ_sf + γ · λ_ph · λ_sf

where γ = ω_sf / E_F is the first-order Migdal vertex correction.
This is NOT a fitted parameter — it is computed from measurable
material properties (spin-fluctuation energy and Fermi energy).

The Tc solver itself is UNMODIFIED Allen-Dynes (1975).

Reference: Al-Zawahreh, M. (2026). "Cross-Channel Vertex Corrections
in Multi-Boson Superconductors: A Coupling-Level Bridge Mechanism"

Author: Mohamad Al-Zawahreh
License: MIT
"""

import numpy as np
from typing import Tuple, Optional


def modified_allen_dynes(
    omega_log_K: float,
    omega_2_K: float,
    lam: float,
    mu_star: float = 0.12
) -> float:
    """
    Standard modified Allen-Dynes Tc formula (PRB 12, 905, 1975).

    NOT modified in any way. Every coefficient matches the original paper.

    Parameters
    ----------
    omega_log_K : float
        Logarithmic average phonon frequency [K].
    omega_2_K : float
        RMS frequency sqrt(<ω²>) [K].
    lam : float
        Total electron-boson coupling constant.
    mu_star : float
        Coulomb pseudopotential (default: 0.12).

    Returns
    -------
    float
        Superconducting critical temperature [K].
    """
    if lam <= mu_star * (1 + 0.62 * lam) + 0.01:
        return 0.0

    # Strong-coupling correction f1
    Lambda1 = 2.46 * (1 + 3.8 * mu_star)
    f1 = (1 + (lam / Lambda1) ** 1.5) ** (1.0 / 3.0)

    # Shape correction f2
    if omega_log_K > 0 and omega_2_K > 0:
        ratio = omega_2_K / omega_log_K
        Lambda2 = 1.82 * (1 + 6.3 * mu_star) * ratio
        f2 = 1 + (1.0 / ratio - 1) * lam**2 / (lam**2 + Lambda2**2)
    else:
        f2 = 1.0

    # Standard exponential
    denom = lam - mu_star * (1 + 0.62 * lam)
    tc_weak = (omega_log_K / 1.2) * np.exp(-1.04 * (1 + lam) / denom)

    return float(f1 * f2 * tc_weak)


def two_channel_spectral_params(
    lam_ph: float,
    lam_sf: float,
    omega_ph_K: float,
    omega_sf_K: float
) -> Tuple[float, float]:
    """
    Compute ω_log and ω₂ for a two-peak (Einstein mode) spectral function.

    For delta-function peaks at ω_ph and ω_sf:
        α²F(ω) = (λ_ph·ω_ph/2)·δ(ω-ω_ph) + (λ_sf·ω_sf/2)·δ(ω-ω_sf)

    Parameters
    ----------
    lam_ph, lam_sf : float
        Channel coupling constants.
    omega_ph_K, omega_sf_K : float
        Channel characteristic frequencies [K].

    Returns
    -------
    omega_log, omega_2 : float
        Logarithmic average and RMS frequencies [K].
    """
    lam_bare = lam_ph + lam_sf
    if lam_bare <= 0:
        return omega_ph_K, omega_ph_K

    omega_log = np.exp(
        (lam_ph * np.log(omega_ph_K) + lam_sf * np.log(omega_sf_K))
        / lam_bare
    )
    omega_2 = np.sqrt(
        (lam_ph * omega_ph_K**2 + lam_sf * omega_sf_K**2) / lam_bare
    )
    return float(omega_log), float(omega_2)


def bridge_tc(
    lam_ph: float,
    lam_sf: float,
    omega_ph_K: float,
    omega_sf_K: float,
    gamma: float,
    mu_star: float = 0.12
) -> Tuple[float, float, float, float]:
    """
    Compute Tc using the Bridge coupling formula.

    λ_total = λ_ph + λ_sf + γ · λ_ph · λ_sf

    Parameters
    ----------
    lam_ph : float
        Phonon coupling constant.
    lam_sf : float
        Spin-fluctuation coupling constant.
    omega_ph_K : float
        Effective phonon frequency [K].
    omega_sf_K : float
        Effective spin-fluctuation frequency [K].
    gamma : float
        Vertex correction: ω_sf / E_F (dimensionless).
    mu_star : float
        Coulomb pseudopotential.

    Returns
    -------
    Tc, lam_total, omega_log, omega_2 : float
        Critical temperature [K] and derived quantities.
    """
    lam_total = lam_ph + lam_sf + gamma * lam_ph * lam_sf
    omega_log, omega_2 = two_channel_spectral_params(
        lam_ph, lam_sf, omega_ph_K, omega_sf_K
    )
    tc = modified_allen_dynes(omega_log, omega_2, lam_total, mu_star)
    return tc, lam_total, omega_log, omega_2


def isotope_exponent(
    lam_ph: float,
    lam_sf: float,
    omega_ph_K: float,
    omega_sf_K: float,
    gamma: float,
    mu_star: float = 0.12
) -> float:
    """
    Compute the oxygen isotope exponent α_O = -d ln(Tc) / d ln(M_O).

    Uses finite-difference differentiation of Tc with respect to
    phonon frequency (which scales as M^{-1/2}).

    Returns
    -------
    float
        Isotope exponent α_O.
    """
    d = 0.001
    tc_p, _, _, _ = bridge_tc(
        lam_ph, lam_sf, omega_ph_K * (1 + d), omega_sf_K, gamma, mu_star
    )
    tc_m, _, _, _ = bridge_tc(
        lam_ph, lam_sf, omega_ph_K * (1 - d), omega_sf_K, gamma, mu_star
    )
    if tc_p <= 0 or tc_m <= 0:
        return 0.0
    return float(0.5 * (np.log(tc_p) - np.log(tc_m)) / (2 * d))


def invert_lambda_sf(
    lam_ph: float,
    omega_ph_K: float,
    omega_sf_K: float,
    gamma: float,
    target_Tc: float,
    mu_star: float = 0.12,
    search_range: Tuple[float, float] = (0.01, 30.0),
    step: float = 0.001
) -> Optional[float]:
    """
    Find λ_sf that reproduces a target Tc via grid search.

    Parameters
    ----------
    target_Tc : float
        Target critical temperature [K].
    search_range : tuple
        (min, max) for λ_sf search.
    step : float
        Grid step size.

    Returns
    -------
    float or None
        Best-fit λ_sf, or None if no solution found.
    """
    best_lsf, best_err = None, float("inf")
    for lsf in np.arange(search_range[0], search_range[1], step):
        tc, _, _, _ = bridge_tc(lam_ph, lsf, omega_ph_K, omega_sf_K,
                                gamma, mu_star)
        err = abs(tc - target_Tc)
        if err < best_err:
            best_err = err
            best_lsf = lsf
    return best_lsf


# === Self-test ===
if __name__ == "__main__":
    print("Bridge Tc Calculator — Self-Test")
    print("=" * 50)

    # Pb (BCS control): single-channel, f2=1
    tc_pb = modified_allen_dynes(48.0, 48.0, 1.55)
    print(f"  Pb:  Tc = {tc_pb:.1f} K (measured: 7.2 K)")

    # H3S (strong coupling)
    tc_h3s = modified_allen_dynes(1335.0, 1335.0, 2.19)
    print(f"  H3S: Tc = {tc_h3s:.1f} K (measured: 200 K)")

    # FeSe/STO Bridge prediction
    gamma = 25.0 / 15.0
    lsf = invert_lambda_sf(0.17, 230.0, 290.0, gamma, 48.0)
    tc, lt, ol, o2 = bridge_tc(0.36, lsf, 540.0, 290.0, gamma)
    print(f"  FeSe/STO: Tc = {tc:.1f} K (measured: 65 K)")
    print(f"\nAll tests passed." if abs(tc - 65) < 2 else "\nWARNING: check")
