"""
AD QSP Model - 67 Year Steady State Simulation
Based on: Madrasi et al. (2021) Alzheimer's & Dementia
"""

import numpy as np
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Parameters from 01_parameters.R
# =============================================================================

# Volumes (L)
V_plasma = 3.0
V_peripheral = 3.0
V_csf = 0.139
V_bisf = 0.261

# Kinetic rate constants
konPP = 1e-3
kcatBACE = 7.2e-3
kcatGamma = 1.2e-3
koffBACE = 120
koffGamma = 0.4
kcleave = 2e-7

# Aggregation
kM2G = 1.4e-5
kG2M_plasma = 1.4e5
kG2M_per = 1.4e5
kG2M_csf = 2.8e-3
kG2M_bisf = 1.4e-8
kG2P = 7.0e-8
kP2G_bisf = 7e-11

# Synthesis (nmol/s)
ksynthAPP_per = 7.28e-2
ksynthAPP_bisf = 2e-3
ksynthBACE_per = 9.0e-4
ksynthBACE_bisf = 1.0e-4
ksynthGamma_per = 6.0e-4
ksynthGamma_bisf = 5.2e-5
ksynthFcR_bisf = 2.5e-5

# Clearance
kclearAPP = 6.93e-5
kclearBACE = 1.2e-5
kclearGamma = 8.0e-6
kclearFcR = 9.6e-5
kclearAbeta_plasma = 1.9e-4
kclearAbeta_bisf = 5.5e-5
kclearAolig_plasma = 1.9e-4
kclearAolig_bisf = 2.2e-8
kclearAplaq = 8.0e-9
kclearBACEs = kclearBACE
kclearCTFb = kclearAPP

# Transport
k12Abeta = 3e-4
k21Abeta = 2.0e-5
k12Aoli = 3e-4
k21Aoli = 2.0e-5
k12BACE = 3e-4
k21BACE = 2.0e-5
k13Abeta = 1.72e-9
k31Abeta = 4.5e-5
k13Aoli = 1.72e-9
k31Aoli = 4.5e-5
k13BACE = 1.72e-9
k31BACE = 4.5e-5
k14Abeta = 1.48e-7
k41Abeta = 1.48e-8
k14Aoli = 1.48e-7
k41Aoli = 1.48e-8
k14BACE = 2.0e-8
k41BACE = 8.0e-5
k43Abeta = 7.5e-5
k43Aoli = 2.25e-6
k43BACE = 1.5e-5

# =============================================================================
# ODE System
# =============================================================================

def ode_abeta(t, y):
    # Unpack state variables
    (APP_per, BACE_per, BACEs_per, Gamma_per, APP_BACE_per, CTFb_per, CTFb_Gamma_per, Abeta_per, Aolig_per,
     APP_bisf, BACE_bisf, BACEs_bisf, Gamma_bisf, APP_BACE_bisf, CTFb_bisf, CTFb_Gamma_bisf, Abeta_bisf, Aolig_bisf, Aplaq_bisf, FcR_bisf,
     BACEs_plasma, Abeta_plasma, Aolig_plasma,
     BACEs_csf, Abeta_csf, Aolig_csf) = y

    # Volume ratios
    r_pla_per = V_plasma / V_peripheral
    r_per_pla = V_peripheral / V_plasma
    r_pla_csf = V_plasma / V_csf
    r_csf_pla = V_csf / V_plasma
    r_pla_bisf = V_plasma / V_bisf
    r_bisf_pla = V_bisf / V_plasma
    r_bisf_csf = V_bisf / V_csf

    # PERIPHERAL
    dAPP_per = ksynthAPP_per/V_peripheral - kclearAPP*APP_per - konPP*APP_per*BACE_per + koffBACE*APP_BACE_per
    dBACE_per = ksynthBACE_per/V_peripheral - kclearBACE*BACE_per - kcleave*BACE_per - konPP*APP_per*BACE_per + koffBACE*APP_BACE_per + kcatBACE*APP_BACE_per
    dBACEs_per = kcleave*BACE_per - kclearBACEs*BACEs_per - k21BACE*BACEs_per + k12BACE*BACEs_plasma*r_pla_per
    dGamma_per = ksynthGamma_per/V_peripheral - kclearGamma*Gamma_per - konPP*CTFb_per*Gamma_per + koffGamma*CTFb_Gamma_per + kcatGamma*CTFb_Gamma_per
    dAPP_BACE_per = konPP*APP_per*BACE_per - koffBACE*APP_BACE_per - kcatBACE*APP_BACE_per
    dCTFb_per = kcatBACE*APP_BACE_per - kclearCTFb*CTFb_per - konPP*CTFb_per*Gamma_per + koffGamma*CTFb_Gamma_per
    dCTFb_Gamma_per = konPP*CTFb_per*Gamma_per - koffGamma*CTFb_Gamma_per - kcatGamma*CTFb_Gamma_per
    dAbeta_per = kcatGamma*CTFb_Gamma_per - kclearAbeta_plasma*Abeta_per - kM2G*Abeta_per + kG2M_per*Aolig_per - k21Abeta*Abeta_per + k12Abeta*Abeta_plasma*r_pla_per
    dAolig_per = kM2G*Abeta_per - kG2M_per*Aolig_per - kclearAolig_plasma*Aolig_per - k21Aoli*Aolig_per + k12Aoli*Aolig_plasma*r_pla_per

    # BRAIN ISF
    dAPP_bisf = ksynthAPP_bisf/V_bisf - kclearAPP*APP_bisf - konPP*APP_bisf*BACE_bisf + koffBACE*APP_BACE_bisf
    dBACE_bisf = ksynthBACE_bisf/V_bisf - kclearBACE*BACE_bisf - kcleave*BACE_bisf - konPP*APP_bisf*BACE_bisf + koffBACE*APP_BACE_bisf + kcatBACE*APP_BACE_bisf
    dBACEs_bisf = kcleave*BACE_bisf - kclearBACEs*BACEs_bisf - k41BACE*BACEs_bisf + k14BACE*BACEs_plasma*r_pla_bisf - k43BACE*BACEs_bisf
    dGamma_bisf = ksynthGamma_bisf/V_bisf - kclearGamma*Gamma_bisf - konPP*CTFb_bisf*Gamma_bisf + koffGamma*CTFb_Gamma_bisf + kcatGamma*CTFb_Gamma_bisf
    dAPP_BACE_bisf = konPP*APP_bisf*BACE_bisf - koffBACE*APP_BACE_bisf - kcatBACE*APP_BACE_bisf
    dCTFb_bisf = kcatBACE*APP_BACE_bisf - kclearCTFb*CTFb_bisf - konPP*CTFb_bisf*Gamma_bisf + koffGamma*CTFb_Gamma_bisf
    dCTFb_Gamma_bisf = konPP*CTFb_bisf*Gamma_bisf - koffGamma*CTFb_Gamma_bisf - kcatGamma*CTFb_Gamma_bisf
    dAbeta_bisf = kcatGamma*CTFb_Gamma_bisf - kclearAbeta_bisf*Abeta_bisf - kM2G*Abeta_bisf + kG2M_bisf*Aolig_bisf - k41Abeta*Abeta_bisf + k14Abeta*Abeta_plasma*r_pla_bisf - k43Abeta*Abeta_bisf
    dAolig_bisf = kM2G*Abeta_bisf - kG2M_bisf*Aolig_bisf - kclearAolig_bisf*Aolig_bisf - kG2P*Aolig_bisf + kP2G_bisf*Aplaq_bisf - k41Aoli*Aolig_bisf + k14Aoli*Aolig_plasma*r_pla_bisf - k43Aoli*Aolig_bisf
    dAplaq_bisf = kG2P*Aolig_bisf - kP2G_bisf*Aplaq_bisf - kclearAplaq*Aplaq_bisf
    dFcR_bisf = ksynthFcR_bisf/V_bisf - kclearFcR*FcR_bisf

    # PLASMA
    dBACEs_plasma = -kclearBACEs*BACEs_plasma - k12BACE*BACEs_plasma + k21BACE*BACEs_per*r_per_pla - k13BACE*BACEs_plasma + k31BACE*BACEs_csf*r_csf_pla - k14BACE*BACEs_plasma + k41BACE*BACEs_bisf*r_bisf_pla
    dAbeta_plasma = -kclearAbeta_plasma*Abeta_plasma - kM2G*Abeta_plasma + kG2M_plasma*Aolig_plasma - k12Abeta*Abeta_plasma + k21Abeta*Abeta_per*r_per_pla - k13Abeta*Abeta_plasma + k31Abeta*Abeta_csf*r_csf_pla - k14Abeta*Abeta_plasma + k41Abeta*Abeta_bisf*r_bisf_pla
    dAolig_plasma = kM2G*Abeta_plasma - kG2M_plasma*Aolig_plasma - kclearAolig_plasma*Aolig_plasma - k12Aoli*Aolig_plasma + k21Aoli*Aolig_per*r_per_pla - k13Aoli*Aolig_plasma + k31Aoli*Aolig_csf*r_csf_pla - k14Aoli*Aolig_plasma + k41Aoli*Aolig_bisf*r_bisf_pla

    # CSF
    dBACEs_csf = k13BACE*BACEs_plasma*r_pla_csf - k31BACE*BACEs_csf + k43BACE*BACEs_bisf*r_bisf_csf
    dAbeta_csf = -kM2G*Abeta_csf + kG2M_csf*Aolig_csf + k13Abeta*Abeta_plasma*r_pla_csf - k31Abeta*Abeta_csf + k43Abeta*Abeta_bisf*r_bisf_csf
    dAolig_csf = kM2G*Abeta_csf - kG2M_csf*Aolig_csf + k13Aoli*Aolig_plasma*r_pla_csf - k31Aoli*Aolig_csf + k43Aoli*Aolig_bisf*r_bisf_csf

    return [dAPP_per, dBACE_per, dBACEs_per, dGamma_per, dAPP_BACE_per, dCTFb_per, dCTFb_Gamma_per, dAbeta_per, dAolig_per,
            dAPP_bisf, dBACE_bisf, dBACEs_bisf, dGamma_bisf, dAPP_BACE_bisf, dCTFb_bisf, dCTFb_Gamma_bisf, dAbeta_bisf, dAolig_bisf, dAplaq_bisf, dFcR_bisf,
            dBACEs_plasma, dAbeta_plasma, dAolig_plasma,
            dBACEs_csf, dAbeta_csf, dAolig_csf]

# =============================================================================
# Run Simulation
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("       AD QSP Model - 67 Year Steady State Simulation")
    print("=" * 70)
    print()

    # Initial conditions (all zeros)
    y0 = np.zeros(26)

    # Time span: 67 years in seconds
    t_end = 67 * 365 * 24 * 3600
    print(f"Simulation duration: 67 years ({t_end:.2e} seconds)")
    print("Running solver (this may take a moment)...")

    # Solve ODE
    sol = solve_ivp(ode_abeta, [0, t_end], y0, method='LSODA', rtol=1e-8, atol=1e-10)

    print(f"Solver completed with {len(sol.t)} time points")
    print()

    # Extract final values
    final = sol.y[:, -1]

    # State variable indices
    idx = {
        'Abeta_per': 7, 'Aolig_per': 8,
        'Abeta_bisf': 16, 'Aolig_bisf': 17, 'Aplaq_bisf': 18,
        'Abeta_plasma': 21, 'Aolig_plasma': 22,
        'Abeta_csf': 24, 'Aolig_csf': 25
    }

    # Expected values from Table S2
    expected = {
        'Abeta_plasma': 0.066,
        'Abeta_csf': 2.8,
        'Abeta_bisf': 0.75,
        'Aolig_plasma': 6.9e-12,
        'Aolig_csf': 0.02,
        'Aolig_bisf': 4.4,
        'Aplaq_bisf': 39
    }

    # Print comparison table
    print("=" * 70)
    print("       STEADY STATE VALIDATION (Table S2 Comparison)")
    print("=" * 70)
    print()
    header = f"{'Species':<15} {'Compart.':<10} {'Expected':<12} {'Simulated':<12} {'Fold':<10} {'Status':<8}"
    print(header)
    print("-" * 70)

    results = []
    for species in ['Abeta_plasma', 'Abeta_csf', 'Abeta_bisf', 'Aolig_plasma', 'Aolig_csf', 'Aolig_bisf', 'Aplaq_bisf']:
        exp_val = expected[species]
        sim_val = final[idx[species]]
        fold = sim_val / exp_val if exp_val != 0 else float('inf')

        # Determine compartment
        if 'plasma' in species:
            comp = 'Plasma'
        elif 'csf' in species:
            comp = 'CSF'
        else:
            comp = 'Brain ISF'

        # Check if within 2-fold
        status = 'OK' if 0.5 <= fold <= 2 else 'ADJUST'

        results.append((species, comp, exp_val, sim_val, fold, status))
        print(f"{species:<15} {comp:<10} {exp_val:<12.2e} {sim_val:<12.2e} {fold:<10.2f} {status:<8}")

    print("-" * 70)

    # Summary
    n_ok = sum(1 for r in results if r[5] == 'OK')
    print(f"\nValidation: {n_ok}/7 species within 2-fold tolerance")
    print()

    # Plaque half-life check
    plaque_t_half = np.log(2) / kclearAplaq / 3600 / 24 / 365
    print(f"Plaque clearance half-life: {plaque_t_half:.2f} years (target: ~2.75 years)")
