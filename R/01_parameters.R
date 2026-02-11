# =============================================================================
# 01_parameters.R
# All model parameters for AD QSP Model
# Based on: Madrasi et al. (2021) Alzheimer's & Dementia - Table S3
# Units: Time in seconds (s), Concentration in nM
# =============================================================================

# =============================================================================
# 1. Compartment Volumes (L)
# =============================================================================

V_plasma     <- 3.0      # Plasma volume (Pearson et al., 2008)
V_peripheral <- 3.0      # Peripheral tissue volume (assumed same as plasma)
V_csf        <- 0.139    # CSF volume (Nau et al., 2010)
V_bisf       <- 0.261    # Brain ISF volume (Shah & Betts, 2012)

# =============================================================================
# 2. Kinetic Rate Constants for Enzyme Reactions
# =============================================================================

# Association rate constants
konPP    <- 1e-3         # APP-BACE and CTFb-Gamma association (1/(nM*s)) [Schlosshauer & Baker, 2004]

# Turnover numbers (catalytic rate constants)
kcatBACE  <- 7.2e-3      # BACE turnover number (1/s) [Ahmed et al., 2010]
kcatGamma <- 1.2e-3      # Gamma-secretase turnover number (1/s) [Kamp et al., 2015]

# Dissociation rate constants
koffBACE  <- 120         # APP-BACE dissociation (1/s) [Baranello et al., 2015; Lin et al., 2000]
koffGamma <- 0.4         # CTFb-Gamma dissociation (1/s) [Chavez-Gutierrez et al., 2012]

# Soluble BACE formation
kcleave <- 2e-7          # Rate constant for sBACE formation from BACE (1/s) [Fitted]

# =============================================================================
# 3. Aggregation Rate Constants
# =============================================================================

# Monomer to Oligomer
kM2G <- 1.4e-5           # Monomer aggregation to oligomer (1/s) [Fitted]

# Oligomer to Monomer (compartment-specific)
kG2M_plasma <- 1.4e5     # Oligomer dissociation in plasma (1/s) [Fitted]
kG2M_per    <- 1.4e5     # Oligomer dissociation in peripheral (1/s) [Fitted]
kG2M_csf    <- 2.8e-3    # Oligomer dissociation in CSF (1/s) [Fitted]
kG2M_bisf   <- 1.4e-8    # Oligomer dissociation in brain ISF (1/s) [Fitted]

# Oligomer to Plaque (Brain ISF only)
kG2P      <- 7.0e-8      # Oligomer to plaque formation (1/s) [Fitted]
kP2G_bisf <- 7e-11       # Plaque to oligomer dissociation (1/s) [Assumed 1000x slower than association]

# =============================================================================
# 4. Synthesis Rate Constants
# =============================================================================

# APP synthesis
ksynthAPP_per  <- 7.28e-2    # APP synthesis in peripheral (nmol/s) [Fitted]
ksynthAPP_bisf <- 2e-3       # APP synthesis in brain ISF (nmol/s) [Fitted]

# BACE synthesis
ksynthBACE_per  <- 9.0e-4    # BACE synthesis in peripheral (nmol/s) [Fitted]
ksynthBACE_bisf <- 1.0e-4    # BACE synthesis in brain ISF (nmol/s) [Fitted]

# Gamma-secretase synthesis
ksynthGamma_per  <- 6.0e-4   # Gamma synthesis in peripheral (nmol/s) [Fitted]
ksynthGamma_bisf <- 5.2e-5   # Gamma synthesis in brain ISF (nmol/s) [Fitted]

# Fc receptor synthesis (Brain ISF only)
ksynthFcR_bisf <- 2.5e-5     # FcR synthesis in brain ISF (nmol/s) [Fitted]

# =============================================================================
# 5. Clearance Rate Constants
# =============================================================================

# Protein clearance
kclearAPP   <- 6.93e-5       # APP degradation (1/s) [Ring et al., 2007]
kclearBACE  <- 1.2e-5        # BACE degradation (1/s) [Huse et al., 2000]
kclearGamma <- 8.0e-6        # Gamma-secretase degradation (1/s) [Dries & Yu, 2008]
kclearFcR   <- 9.6e-5        # FcR degradation in brain ISF (1/s) [Mellman et al., 1983]

# Abeta monomer clearance
kclearAbeta_plasma <- 1.9e-4 # Abeta monomer clearance in plasma (1/s) [Ovod et al., 2017]
kclearAbeta_bisf   <- 5.5e-5 # Abeta monomer clearance in brain ISF (1/s) [Cirrito et al., 2003]

# Abeta oligomer clearance
kclearAolig_plasma <- 1.9e-4 # Abeta oligomer clearance in plasma (1/s) [Assumed same as monomer]
kclearAolig_bisf   <- 2.2e-8 # Abeta oligomer clearance in brain ISF (1/s) [Fitted]

# Abeta plaque clearance
kclearAplaq <- 8.0e-9        # Abeta plaque clearance in brain ISF (1/s) [Fitted]

# Soluble BACE clearance
kclearBACEs <- kclearBACE    # sBACE clearance (assumed same as BACE)

# CTFb clearance
kclearCTFb <- kclearAPP      # CTFb clearance (assumed same as APP)

# =============================================================================
# 6. Transport Rate Constants (Abeta and sBACE)
# =============================================================================

# --- Plasma <-> Peripheral ---
k12Abeta <- 3e-4             # Plasma to peripheral (Abeta) (1/s) [Fitted]
k21Abeta <- 2.0e-5           # Peripheral to plasma (Abeta) (1/s) [Fitted]

k12Aoli  <- 3e-4             # Plasma to peripheral (Aolig) (1/s) [Assumed same as Abeta]
k21Aoli  <- 2.0e-5           # Peripheral to plasma (Aolig) (1/s) [Assumed same as Abeta]

k12BACE  <- 3e-4             # Plasma to peripheral (sBACE) (1/s) [Assumed same as Abeta]
k21BACE  <- 2.0e-5           # Peripheral to plasma (sBACE) (1/s) [Fitted]

# --- Plasma <-> CSF ---
k13Abeta <- 1.72e-9          # Plasma to CSF (Abeta) (1/s) [Assumed same as mAb]
k31Abeta <- 4.5e-5           # CSF to plasma (Abeta) (1/s) [Fitted]
# Note: Healthy subjects have k31Abeta = 1.4 * 4.5e-5

k13Aoli  <- 1.72e-9          # Plasma to CSF (Aolig) (1/s) [Assumed same as mAb]
k31Aoli  <- 4.5e-5           # CSF to plasma (Aolig) (1/s) [Assumed same as Abeta]

k13BACE  <- 1.72e-9          # Plasma to CSF (sBACE) (1/s) [Assumed same as mAb]
k31BACE  <- 4.5e-5           # CSF to plasma (sBACE) (1/s) [Assumed same as Abeta]

# --- Plasma <-> Brain ISF ---
k14Abeta <- 1.48e-7          # Plasma to brain ISF (Abeta) (1/s) [Fitted]
k41Abeta <- 1.48e-8          # Brain ISF to plasma (Abeta) (1/s) [Fitted]

k14Aoli  <- 1.48e-7          # Plasma to brain ISF (Aolig) (1/s) [Assumed same as Abeta]
k41Aoli  <- 1.48e-8          # Brain ISF to plasma (Aolig) (1/s) [Assumed same as Abeta]

k14BACE  <- 2.0e-8           # Plasma to brain ISF (sBACE) (1/s) [Fitted]
k41BACE  <- 8.0e-5           # Brain ISF to plasma (sBACE) (1/s) [Fitted to Bernier et al., 2013]

# --- Brain ISF -> CSF (unidirectional) ---
k43Abeta <- 7.5e-5           # Brain ISF to CSF (Abeta) (1/s) [Fitted to SILK data]
k43Aoli  <- 2.25e-6          # Brain ISF to CSF (Aolig) (1/s) [Fitted]
k43BACE  <- 1.5e-5           # Brain ISF to CSF (sBACE) (1/s) [Fitted]

# =============================================================================
# 7. Steady State Concentrations (Reference values from Table S2)
# =============================================================================

# These are target values for model validation
SS_Abeta_plasma <- 6.6e-2    # Abeta monomer in plasma (nM)
SS_Aolig_plasma <- 6.9e-12   # Abeta oligomer in plasma (nM)

SS_Abeta_csf    <- 2.8       # Abeta monomer in CSF (nM)
SS_Aolig_csf    <- 2.0e-2    # Abeta oligomer in CSF (nM)

SS_Abeta_bisf   <- 0.75      # Abeta monomer in brain ISF (nM)
SS_Aolig_bisf   <- 4.4       # Abeta oligomer in brain ISF (nM)
SS_Aplaq_bisf   <- 39        # Abeta plaque in brain ISF (nM)

# =============================================================================
# 8. ADCP Parameters (common values, drug-specific in 04_drug_models.R)
# =============================================================================

konPF  <- 0.001              # Fc:FcgR binding on rate (1/(nM*s))
cFcR_bisf <- 1               # Fc receptor concentration in brain ISF (nM)

# =============================================================================
# 9. Helper Function: Calculate derived parameters
# =============================================================================

# Calculate Kd from kon and koff
calc_Kd <- function(koff, kon) {
  return(koff / kon)
}

# Calculate koff from Kd and kon
calc_koff <- function(Kd, kon) {
  return(Kd * kon)
}

# =============================================================================
# 10. Population-specific adjustments
# =============================================================================

# For healthy volunteers (vs AD patients)
# k31Abeta_healthy = 1.4 * k31Abeta  (faster transport from CSF to circulation)

get_healthy_params <- function() {
  params <- get_AD_params()
  params["k31Abeta"] <- params["k31Abeta"] * 1.4
  return(params)
}

# =============================================================================
# 11. Compile all parameters into a named vector
# =============================================================================

get_AD_params <- function() {
  params <- c(
    # Volumes
    V_plasma     = V_plasma,
    V_peripheral = V_peripheral,
    V_csf        = V_csf,
    V_bisf       = V_bisf,

    # Kinetic rate constants
    konPP     = konPP,
    kcatBACE  = kcatBACE,
    kcatGamma = kcatGamma,
    koffBACE  = koffBACE,
    koffGamma = koffGamma,
    kcleave   = kcleave,

    # Aggregation
    kM2G       = kM2G,
    kG2M_plasma = kG2M_plasma,
    kG2M_per   = kG2M_per,
    kG2M_csf   = kG2M_csf,
    kG2M_bisf  = kG2M_bisf,
    kG2P       = kG2P,
    kP2G_bisf  = kP2G_bisf,

    # Synthesis
    ksynthAPP_per    = ksynthAPP_per,
    ksynthAPP_bisf   = ksynthAPP_bisf,
    ksynthBACE_per   = ksynthBACE_per,
    ksynthBACE_bisf  = ksynthBACE_bisf,
    ksynthGamma_per  = ksynthGamma_per,
    ksynthGamma_bisf = ksynthGamma_bisf,
    ksynthFcR_bisf   = ksynthFcR_bisf,

    # Clearance
    kclearAPP          = kclearAPP,
    kclearBACE         = kclearBACE,
    kclearGamma        = kclearGamma,
    kclearFcR          = kclearFcR,
    kclearAbeta_plasma = kclearAbeta_plasma,
    kclearAbeta_bisf   = kclearAbeta_bisf,
    kclearAolig_plasma = kclearAolig_plasma,
    kclearAolig_bisf   = kclearAolig_bisf,
    kclearAplaq        = kclearAplaq,
    kclearBACEs        = kclearBACEs,
    kclearCTFb         = kclearCTFb,

    # Transport: Plasma <-> Peripheral
    k12Abeta = k12Abeta,
    k21Abeta = k21Abeta,
    k12Aoli  = k12Aoli,
    k21Aoli  = k21Aoli,
    k12BACE  = k12BACE,
    k21BACE  = k21BACE,

    # Transport: Plasma <-> CSF
    k13Abeta = k13Abeta,
    k31Abeta = k31Abeta,
    k13Aoli  = k13Aoli,
    k31Aoli  = k31Aoli,
    k13BACE  = k13BACE,
    k31BACE  = k31BACE,

    # Transport: Plasma <-> Brain ISF
    k14Abeta = k14Abeta,
    k41Abeta = k41Abeta,
    k14Aoli  = k14Aoli,
    k41Aoli  = k41Aoli,
    k14BACE  = k14BACE,
    k41BACE  = k41BACE,

    # Transport: Brain ISF -> CSF
    k43Abeta = k43Abeta,
    k43Aoli  = k43Aoli,
    k43BACE  = k43BACE,

    # ADCP common parameters
    konPF     = konPF,
    cFcR_bisf = cFcR_bisf
  )

  return(params)
}

# =============================================================================
# Print summary of parameters
# =============================================================================

print_params_summary <- function() {
  cat("=== AD QSP Model Parameters Summary ===\n\n")
  cat("Compartment Volumes (L):\n")
  cat(sprintf("  Plasma: %.3f, Peripheral: %.3f, CSF: %.3f, Brain ISF: %.3f\n\n",
              V_plasma, V_peripheral, V_csf, V_bisf))

  cat("Key Rate Constants:\n")
  cat(sprintf("  kcatBACE: %.2e /s, kcatGamma: %.2e /s\n", kcatBACE, kcatGamma))
  cat(sprintf("  kM2G (aggregation): %.2e /s\n", kM2G))
  cat(sprintf("  kG2P (plaque formation): %.2e /s\n", kG2P))
  cat(sprintf("  kclearAplaq: %.2e /s (t1/2 = %.2f years)\n\n",
              kclearAplaq, log(2)/kclearAplaq/3600/24/365))

  cat("Steady State Targets (nM):\n")
  cat(sprintf("  Abeta - Plasma: %.3f, CSF: %.1f, Brain: %.2f\n",
              SS_Abeta_plasma, SS_Abeta_csf, SS_Abeta_bisf))
  cat(sprintf("  Plaque - Brain: %.0f nM\n", SS_Aplaq_bisf))
}
