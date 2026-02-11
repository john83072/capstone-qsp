# =============================================================================
# 02_model_odes.R
# ODE system definition for AD QSP Model - Basic Aβ Biology (No Drugs)
# Based on: Madrasi et al. (2021) Alzheimer's & Dementia
# =============================================================================

library(deSolve)

# =============================================================================
# 1. State Variable Names and Indices
# =============================================================================

# Define state variable names for basic Aβ biology
state_names_biology <- c(
 # Peripheral Compartment
 "APP_per",          # 1: APP in peripheral
 "BACE_per",         # 2: Membrane-bound BACE in peripheral
 "BACEs_per",        # 3: Soluble BACE in peripheral
 "Gamma_per",        # 4: γ-secretase in peripheral
 "APP_BACE_per",     # 5: APP:BACE complex in peripheral
 "CTFb_per",         # 6: CTFβ in peripheral
 "CTFb_Gamma_per",   # 7: CTFβ:γ complex in peripheral
 "Abeta_per",        # 8: Aβ monomer in peripheral
 "Aolig_per",        # 9: Aβ oligomer in peripheral

 # Brain ISF Compartment
 "APP_bisf",         # 10: APP in brain ISF
 "BACE_bisf",        # 11: Membrane-bound BACE in brain ISF
 "BACEs_bisf",       # 12: Soluble BACE in brain ISF
 "Gamma_bisf",       # 13: γ-secretase in brain ISF
 "APP_BACE_bisf",    # 14: APP:BACE complex in brain ISF
 "CTFb_bisf",        # 15: CTFβ in brain ISF
 "CTFb_Gamma_bisf",  # 16: CTFβ:γ complex in brain ISF
 "Abeta_bisf",       # 17: Aβ monomer in brain ISF
 "Aolig_bisf",       # 18: Aβ oligomer in brain ISF
 "Aplaq_bisf",       # 19: Aβ plaque in brain ISF
 "FcR_bisf",         # 20: Fc receptor in brain ISF

 # Plasma Compartment
 "BACEs_plasma",     # 21: Soluble BACE in plasma
 "Abeta_plasma",     # 22: Aβ monomer in plasma
 "Aolig_plasma",     # 23: Aβ oligomer in plasma

 # CSF Compartment
 "BACEs_csf",        # 24: Soluble BACE in CSF
 "Abeta_csf",        # 25: Aβ monomer in CSF
 "Aolig_csf"         # 26: Aβ oligomer in CSF
)

n_states_biology <- length(state_names_biology)

# =============================================================================
# 2. ODE Function for Basic Aβ Biology (No Drugs)
# =============================================================================

ode_abeta_biology <- function(t, state, parameters) {

 with(as.list(c(state, parameters)), {

   # =========================================================================
   # Volume ratios for transport (concentration-based)
   # =========================================================================
   ratio_plasma_per  <- V_plasma / V_peripheral
   ratio_per_plasma  <- V_peripheral / V_plasma
   ratio_plasma_csf  <- V_plasma / V_csf
   ratio_csf_plasma  <- V_csf / V_plasma
   ratio_plasma_bisf <- V_plasma / V_bisf
   ratio_bisf_plasma <- V_bisf / V_plasma
   ratio_bisf_csf    <- V_bisf / V_csf

   # =========================================================================
   # PERIPHERAL COMPARTMENT
   # =========================================================================

   # --- APP ---
   # Synthesis, clearance, binding to BACE, release from complex
   dAPP_per <- ksynthAPP_per / V_peripheral -        # Synthesis (nmol/s -> nM/s)
               kclearAPP * APP_per -                  # Clearance
               konPP * APP_per * BACE_per +           # Binding to BACE
               koffBACE * APP_BACE_per                # Dissociation from complex

   # --- BACE (membrane-bound) ---
   # Synthesis, clearance, cleavage to sBACE, binding/release with APP, catalysis
   dBACE_per <- ksynthBACE_per / V_peripheral -      # Synthesis
                kclearBACE * BACE_per -               # Clearance
                kcleave * BACE_per -                  # Cleavage to soluble BACE
                konPP * APP_per * BACE_per +          # Binding to APP
                koffBACE * APP_BACE_per +             # Dissociation from complex
                kcatBACE * APP_BACE_per               # Release after catalysis

   # --- Soluble BACE (BACEs) ---
   # Formation from BACE, clearance, transport
   dBACEs_per <- kcleave * BACE_per -                # Formation from membrane BACE
                 kclearBACEs * BACEs_per -            # Clearance
                 k21BACE * BACEs_per +                # Transport to plasma
                 k12BACE * BACEs_plasma * ratio_plasma_per  # Transport from plasma

   # --- γ-secretase (Gamma) ---
   # Synthesis, clearance, binding/release with CTFβ, catalysis
   dGamma_per <- ksynthGamma_per / V_peripheral -    # Synthesis
                 kclearGamma * Gamma_per -            # Clearance
                 konPP * CTFb_per * Gamma_per +       # Binding to CTFβ
                 koffGamma * CTFb_Gamma_per +         # Dissociation
                 kcatGamma * CTFb_Gamma_per           # Release after catalysis

   # --- APP:BACE complex ---
   # Formation, dissociation, catalysis
   dAPP_BACE_per <- konPP * APP_per * BACE_per -     # Formation
                    koffBACE * APP_BACE_per -         # Dissociation
                    kcatBACE * APP_BACE_per           # Catalysis (produces CTFβ)

   # --- CTFβ ---
   # Production from APP:BACE, clearance, binding to γ-secretase
   dCTFb_per <- kcatBACE * APP_BACE_per -            # Production
                kclearCTFb * CTFb_per -               # Clearance
                konPP * CTFb_per * Gamma_per +        # Binding to γ-secretase
                koffGamma * CTFb_Gamma_per            # Dissociation

   # --- CTFβ:γ complex ---
   # Formation, dissociation, catalysis
   dCTFb_Gamma_per <- konPP * CTFb_per * Gamma_per - # Formation
                      koffGamma * CTFb_Gamma_per -    # Dissociation
                      kcatGamma * CTFb_Gamma_per      # Catalysis (produces Aβ)

   # --- Aβ monomer (peripheral) ---
   # Production from CTFβ:γ, clearance, aggregation, transport
   dAbeta_per <- kcatGamma * CTFb_Gamma_per -        # Production
                 kclearAbeta_plasma * Abeta_per -     # Clearance (uses plasma rate)
                 kM2G * Abeta_per +                   # Aggregation to oligomer
                 kG2M_per * Aolig_per -               # Dissociation from oligomer
                 k21Abeta * Abeta_per +               # Transport to plasma
                 k12Abeta * Abeta_plasma * ratio_plasma_per  # Transport from plasma

   # --- Aβ oligomer (peripheral) ---
   # Formation from monomer, dissociation, clearance, transport
   dAolig_per <- kM2G * Abeta_per -                  # Formation from monomer
                 kG2M_per * Aolig_per -               # Dissociation to monomer
                 kclearAolig_plasma * Aolig_per -     # Clearance
                 k21Aoli * Aolig_per +                # Transport to plasma
                 k12Aoli * Aolig_plasma * ratio_plasma_per  # Transport from plasma

   # =========================================================================
   # BRAIN ISF COMPARTMENT
   # =========================================================================

   # --- APP ---
   dAPP_bisf <- ksynthAPP_bisf / V_bisf -            # Synthesis
                kclearAPP * APP_bisf -                # Clearance
                konPP * APP_bisf * BACE_bisf +        # Binding to BACE
                koffBACE * APP_BACE_bisf              # Dissociation

   # --- BACE (membrane-bound) ---
   dBACE_bisf <- ksynthBACE_bisf / V_bisf -          # Synthesis
                 kclearBACE * BACE_bisf -             # Clearance
                 kcleave * BACE_bisf -                # Cleavage to sBACE
                 konPP * APP_bisf * BACE_bisf +       # Binding to APP
                 koffBACE * APP_BACE_bisf +           # Dissociation
                 kcatBACE * APP_BACE_bisf             # Release after catalysis

   # --- Soluble BACE ---
   dBACEs_bisf <- kcleave * BACE_bisf -              # Formation
                  kclearBACEs * BACEs_bisf -          # Clearance
                  k41BACE * BACEs_bisf +              # Transport to plasma
                  k14BACE * BACEs_plasma * ratio_plasma_bisf - # From plasma
                  k43BACE * BACEs_bisf                # Transport to CSF

   # --- γ-secretase ---
   dGamma_bisf <- ksynthGamma_bisf / V_bisf -        # Synthesis
                  kclearGamma * Gamma_bisf -          # Clearance
                  konPP * CTFb_bisf * Gamma_bisf +    # Binding
                  koffGamma * CTFb_Gamma_bisf +       # Dissociation
                  kcatGamma * CTFb_Gamma_bisf         # Release after catalysis

   # --- APP:BACE complex ---
   dAPP_BACE_bisf <- konPP * APP_bisf * BACE_bisf -  # Formation
                     koffBACE * APP_BACE_bisf -       # Dissociation
                     kcatBACE * APP_BACE_bisf         # Catalysis

   # --- CTFβ ---
   dCTFb_bisf <- kcatBACE * APP_BACE_bisf -          # Production
                 kclearCTFb * CTFb_bisf -             # Clearance
                 konPP * CTFb_bisf * Gamma_bisf +     # Binding
                 koffGamma * CTFb_Gamma_bisf          # Dissociation

   # --- CTFβ:γ complex ---
   dCTFb_Gamma_bisf <- konPP * CTFb_bisf * Gamma_bisf - # Formation
                       koffGamma * CTFb_Gamma_bisf -     # Dissociation
                       kcatGamma * CTFb_Gamma_bisf       # Catalysis

   # --- Aβ monomer (brain ISF) ---
   dAbeta_bisf <- kcatGamma * CTFb_Gamma_bisf -      # Production
                  kclearAbeta_bisf * Abeta_bisf -     # Clearance
                  kM2G * Abeta_bisf +                 # Aggregation
                  kG2M_bisf * Aolig_bisf -            # Dissociation from oligomer
                  k41Abeta * Abeta_bisf +             # Transport to plasma
                  k14Abeta * Abeta_plasma * ratio_plasma_bisf - # From plasma
                  k43Abeta * Abeta_bisf               # Transport to CSF

   # --- Aβ oligomer (brain ISF) ---
   dAolig_bisf <- kM2G * Abeta_bisf -                # Formation from monomer
                  kG2M_bisf * Aolig_bisf -            # Dissociation to monomer
                  kclearAolig_bisf * Aolig_bisf -     # Clearance
                  kG2P * Aolig_bisf +                 # Formation of plaque
                  kP2G_bisf * Aplaq_bisf -            # Dissociation from plaque
                  k41Aoli * Aolig_bisf +              # Transport to plasma
                  k14Aoli * Aolig_plasma * ratio_plasma_bisf - # From plasma
                  k43Aoli * Aolig_bisf                # Transport to CSF

   # --- Aβ plaque (brain ISF only) ---
   dAplaq_bisf <- kG2P * Aolig_bisf -                # Formation from oligomer
                  kP2G_bisf * Aplaq_bisf -            # Dissociation to oligomer
                  kclearAplaq * Aplaq_bisf            # Clearance

   # --- Fc Receptor (brain ISF) ---
   # For future drug models; included in biology for steady state
   dFcR_bisf <- ksynthFcR_bisf / V_bisf -            # Synthesis
                kclearFcR * FcR_bisf                  # Clearance

   # =========================================================================
   # PLASMA COMPARTMENT
   # =========================================================================

   # --- Soluble BACE ---
   dBACEs_plasma <- -kclearBACEs * BACEs_plasma -    # Clearance
                    k12BACE * BACEs_plasma +          # To peripheral
                    k21BACE * BACEs_per * ratio_per_plasma - # From peripheral
                    k13BACE * BACEs_plasma +          # To CSF
                    k31BACE * BACEs_csf * ratio_csf_plasma - # From CSF
                    k14BACE * BACEs_plasma +          # To brain ISF
                    k41BACE * BACEs_bisf * ratio_bisf_plasma # From brain ISF

   # --- Aβ monomer (plasma) ---
   dAbeta_plasma <- -kclearAbeta_plasma * Abeta_plasma - # Clearance
                    kM2G * Abeta_plasma +             # Aggregation
                    kG2M_plasma * Aolig_plasma -      # Dissociation
                    k12Abeta * Abeta_plasma +         # To peripheral
                    k21Abeta * Abeta_per * ratio_per_plasma - # From peripheral
                    k13Abeta * Abeta_plasma +         # To CSF
                    k31Abeta * Abeta_csf * ratio_csf_plasma - # From CSF
                    k14Abeta * Abeta_plasma +         # To brain ISF
                    k41Abeta * Abeta_bisf * ratio_bisf_plasma # From brain ISF

   # --- Aβ oligomer (plasma) ---
   dAolig_plasma <- kM2G * Abeta_plasma -            # Formation
                    kG2M_plasma * Aolig_plasma -      # Dissociation
                    kclearAolig_plasma * Aolig_plasma - # Clearance
                    k12Aoli * Aolig_plasma +          # To peripheral
                    k21Aoli * Aolig_per * ratio_per_plasma - # From peripheral
                    k13Aoli * Aolig_plasma +          # To CSF
                    k31Aoli * Aolig_csf * ratio_csf_plasma - # From CSF
                    k14Aoli * Aolig_plasma +          # To brain ISF
                    k41Aoli * Aolig_bisf * ratio_bisf_plasma # From brain ISF

   # =========================================================================
   # CSF COMPARTMENT
   # =========================================================================

   # --- Soluble BACE ---
   dBACEs_csf <- k13BACE * BACEs_plasma * ratio_plasma_csf - # From plasma
                 k31BACE * BACEs_csf +                # To plasma
                 k43BACE * BACEs_bisf * ratio_bisf_csf # From brain ISF

   # --- Aβ monomer (CSF) ---
   dAbeta_csf <- -kM2G * Abeta_csf +                 # Aggregation
                 kG2M_csf * Aolig_csf +               # Dissociation
                 k13Abeta * Abeta_plasma * ratio_plasma_csf - # From plasma
                 k31Abeta * Abeta_csf +               # To plasma
                 k43Abeta * Abeta_bisf * ratio_bisf_csf # From brain ISF

   # --- Aβ oligomer (CSF) ---
   dAolig_csf <- kM2G * Abeta_csf -                  # Formation
                 kG2M_csf * Aolig_csf +               # Dissociation
                 k13Aoli * Aolig_plasma * ratio_plasma_csf - # From plasma
                 k31Aoli * Aolig_csf +                # To plasma
                 k43Aoli * Aolig_bisf * ratio_bisf_csf # From brain ISF

   # =========================================================================
   # Return derivatives
   # =========================================================================

   return(list(c(
     # Peripheral
     dAPP_per, dBACE_per, dBACEs_per, dGamma_per,
     dAPP_BACE_per, dCTFb_per, dCTFb_Gamma_per,
     dAbeta_per, dAolig_per,
     # Brain ISF
     dAPP_bisf, dBACE_bisf, dBACEs_bisf, dGamma_bisf,
     dAPP_BACE_bisf, dCTFb_bisf, dCTFb_Gamma_bisf,
     dAbeta_bisf, dAolig_bisf, dAplaq_bisf, dFcR_bisf,
     # Plasma
     dBACEs_plasma, dAbeta_plasma, dAolig_plasma,
     # CSF
     dBACEs_csf, dAbeta_csf, dAolig_csf
   )))
 })
}

# =============================================================================
# 3. Initial Conditions (all zeros for de novo simulation)
# =============================================================================

get_initial_state_biology <- function() {
 state <- rep(0, n_states_biology)
 names(state) <- state_names_biology
 return(state)
}

# =============================================================================
# 4. Helper function to extract results by compartment
# =============================================================================

extract_compartment_results <- function(out) {
 # Convert to data frame
 df <- as.data.frame(out)

 # Calculate total Aβ (monomer + oligomer) in each compartment
 df$total_Abeta_plasma <- df$Abeta_plasma + df$Aolig_plasma
 df$total_Abeta_csf    <- df$Abeta_csf + df$Aolig_csf
 df$total_Abeta_bisf   <- df$Abeta_bisf + df$Aolig_bisf
 df$total_Abeta_per    <- df$Abeta_per + df$Aolig_per

 return(df)
}

# =============================================================================
# 5. Function to check mass balance (for debugging)
# =============================================================================

check_mass_balance <- function(out, params) {
 df <- as.data.frame(out)
 n <- nrow(df)

 # Total Aβ in system (as monomer equivalents)
 # Note: Plaque is also counted as Aβ mass
 total_Abeta <- with(df, {
   (Abeta_plasma + Aolig_plasma) * params["V_plasma"] +
   (Abeta_per + Aolig_per) * params["V_peripheral"] +
   (Abeta_csf + Aolig_csf) * params["V_csf"] +
   (Abeta_bisf + Aolig_bisf + Aplaq_bisf) * params["V_bisf"]
 })

 cat("Total Aβ mass (nmol) at t=0:", total_Abeta[1], "\n")
 cat("Total Aβ mass (nmol) at t=end:", total_Abeta[n], "\n")

 return(total_Abeta)
}

# =============================================================================
# 6. Quick test function
# =============================================================================

test_biology_model <- function(t_end = 3600*24*7, dt = 3600) {
 # Run for 1 week by default

 # Get parameters and initial state
 params <- get_AD_params()
 state <- get_initial_state_biology()

 # Time points
 times <- seq(0, t_end, by = dt)

 # Solve ODE
 out <- ode(y = state, times = times, func = ode_abeta_biology,
            parms = params, method = "lsoda")

 # Convert to data frame
 df <- extract_compartment_results(out)

 # Print final concentrations
 cat("\n=== Test Results (t =", t_end/3600/24, "days) ===\n")
 cat("Abeta_plasma:", tail(df$Abeta_plasma, 1), "nM\n")
 cat("Abeta_csf:", tail(df$Abeta_csf, 1), "nM\n")
 cat("Abeta_bisf:", tail(df$Abeta_bisf, 1), "nM\n")
 cat("Aolig_bisf:", tail(df$Aolig_bisf, 1), "nM\n")
 cat("Aplaq_bisf:", tail(df$Aplaq_bisf, 1), "nM\n")

 return(df)
}
