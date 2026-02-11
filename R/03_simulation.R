# =============================================================================
# 03_simulation.R
# Simulation functions for AD QSP Model
# Based on: Madrasi et al. (2021) Alzheimer's & Dementia
# =============================================================================

library(deSolve)

# =============================================================================
# 1. Time Conversion Constants
# =============================================================================

SECONDS_PER_MINUTE <- 60
SECONDS_PER_HOUR   <- 3600
SECONDS_PER_DAY    <- 86400
SECONDS_PER_WEEK   <- 604800
SECONDS_PER_YEAR   <- 31536000  # 365 days

# Conversion functions
years_to_seconds <- function(years) years * SECONDS_PER_YEAR
days_to_seconds  <- function(days) days * SECONDS_PER_DAY
weeks_to_seconds <- function(weeks) weeks * SECONDS_PER_WEEK
hours_to_seconds <- function(hours) hours * SECONDS_PER_HOUR

seconds_to_years <- function(seconds) seconds / SECONDS_PER_YEAR
seconds_to_days  <- function(seconds) seconds / SECONDS_PER_DAY
seconds_to_weeks <- function(seconds) seconds / SECONDS_PER_WEEK

# =============================================================================
# 2. Steady State Simulation (67 years)
# =============================================================================

#' Run simulation to reach steady state (67 years)
#'
#' @param params Named vector of model parameters (from get_AD_params())
#' @param t_years Duration in years (default: 67)
#' @param n_timepoints Number of output time points (default: 1000)
#' @param verbose Print progress messages (default: TRUE)
#' @return List with simulation output and final state
#'
run_to_steady_state <- function(params = NULL,
                                 t_years = 67,
                                 n_timepoints = 1000,
                                 verbose = TRUE) {

  # Use default parameters if not provided
  if (is.null(params)) {
    params <- get_AD_params()
  }

  # Convert time to seconds
  t_end <- years_to_seconds(t_years)

  # Time points (log-spaced for long simulations)
  # Use log spacing to capture early dynamics and final steady state
  times <- c(0, exp(seq(log(SECONDS_PER_HOUR), log(t_end), length.out = n_timepoints - 1)))

  # Initial conditions (all zeros)
  state <- get_initial_state_biology()

  if (verbose) {
    cat("=== Running Steady State Simulation ===\n")
    cat("Duration:", t_years, "years (", format(t_end, scientific = TRUE), "seconds)\n")
    cat("Time points:", n_timepoints, "\n")
    cat("Running ODE solver...\n")
  }

  # Solve ODE system
  start_time <- Sys.time()

  out <- ode(
    y = state,
    times = times,
    func = ode_abeta_biology,
    parms = params,
    method = "lsoda",        # Adaptive method for stiff systems
    rtol = 1e-8,             # Relative tolerance
    atol = 1e-10             # Absolute tolerance
  )

  end_time <- Sys.time()

  if (verbose) {
    cat("Completed in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")
  }

  # Extract results
  df <- extract_compartment_results(out)

  # Get final state for drug simulations
  final_state <- as.numeric(out[nrow(out), -1])  # Remove time column
  names(final_state) <- state_names_biology

  return(list(
    output = df,
    final_state = final_state,
    params = params,
    t_years = t_years
  ))
}

# =============================================================================
# 3. Validate Steady State Against Table S2
# =============================================================================

#' Expected steady state values from Table S2
#' @return Named list of expected concentrations (nM)
get_expected_steady_state <- function() {
  list(
    # Aβ Monomer
    Abeta_plasma = 0.066,      # 6.6e-2 nM
    Abeta_csf    = 2.8,        # nM
    Abeta_bisf   = 0.75,       # nM

    # Aβ Oligomer
    Aolig_plasma = 6.9e-12,    # nM
    Aolig_csf    = 0.02,       # 2.0e-2 nM
    Aolig_bisf   = 4.4,        # nM

    # Aβ Plaque (Brain ISF only)
    Aplaq_bisf   = 39          # nM
  )
}

#' Validate simulation results against Table S2
#'
#' @param sim_result Result from run_to_steady_state()
#' @param tolerance_fold Acceptable fold difference (default: 2 = within 2-fold)
#' @return Data frame comparing simulated vs expected values
#'
validate_steady_state <- function(sim_result, tolerance_fold = 2) {

  # Get final (steady state) values
  df <- sim_result$output
  final <- tail(df, 1)

  # Expected values
  expected <- get_expected_steady_state()

  # Build comparison table
  comparison <- data.frame(
    Species = c("Abeta_plasma", "Abeta_csf", "Abeta_bisf",
                "Aolig_plasma", "Aolig_csf", "Aolig_bisf",
                "Aplaq_bisf"),
    Compartment = c("Plasma", "CSF", "Brain ISF",
                    "Plasma", "CSF", "Brain ISF",
                    "Brain ISF"),
    Type = c(rep("Monomer", 3), rep("Oligomer", 3), "Plaque"),
    Expected_nM = c(expected$Abeta_plasma, expected$Abeta_csf, expected$Abeta_bisf,
                    expected$Aolig_plasma, expected$Aolig_csf, expected$Aolig_bisf,
                    expected$Aplaq_bisf),
    Simulated_nM = c(final$Abeta_plasma, final$Abeta_csf, final$Abeta_bisf,
                     final$Aolig_plasma, final$Aolig_csf, final$Aolig_bisf,
                     final$Aplaq_bisf),
    stringsAsFactors = FALSE
  )

  # Calculate fold difference
  comparison$Fold_Diff <- comparison$Simulated_nM / comparison$Expected_nM
  comparison$Log10_Ratio <- log10(comparison$Fold_Diff)

  # Check if within tolerance
  comparison$Within_Tolerance <- abs(comparison$Log10_Ratio) <= log10(tolerance_fold)

  # Print results
  cat("\n")
  cat("============================================================\n")
  cat("       STEADY STATE VALIDATION (Table S2 Comparison)        \n")
  cat("============================================================\n\n")

  cat(sprintf("%-15s %-10s %-12s %-12s %-10s %-8s\n",
              "Species", "Compart.", "Expected", "Simulated", "Fold", "Pass?"))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (i in 1:nrow(comparison)) {
    pass_str <- ifelse(comparison$Within_Tolerance[i], "OK", "FAIL")
    cat(sprintf("%-15s %-10s %-12.2e %-12.2e %-10.2f %-8s\n",
                comparison$Species[i],
                comparison$Compartment[i],
                comparison$Expected_nM[i],
                comparison$Simulated_nM[i],
                comparison$Fold_Diff[i],
                pass_str))
  }

  cat(paste(rep("-", 70), collapse = ""), "\n")

  # Summary
  n_pass <- sum(comparison$Within_Tolerance)
  n_total <- nrow(comparison)
  cat(sprintf("\nValidation: %d/%d species within %.1f-fold tolerance\n",
              n_pass, n_total, tolerance_fold))

  if (n_pass == n_total) {
    cat("STATUS: ALL PASSED\n")
  } else {
    cat("STATUS: SOME FAILED - Parameter adjustment may be needed\n")
  }

  return(comparison)
}

# =============================================================================
# 4. Plot Steady State Approach
# =============================================================================

#' Plot time course to steady state
#'
#' @param sim_result Result from run_to_steady_state()
#' @param log_y Use log scale for y-axis (default: TRUE)
#'
plot_steady_state_approach <- function(sim_result, log_y = TRUE) {

  df <- sim_result$output

  # Convert time to years
  df$time_years <- df$time / SECONDS_PER_YEAR

  # Set up plot layout
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  # --- Plot 1: Aβ Monomer ---
  if (log_y) {
    plot(df$time_years, df$Abeta_bisf, type = "l", col = "blue", lwd = 2,
         log = "y", xlab = "Time (years)", ylab = "Concentration (nM)",
         main = "Aβ Monomer", ylim = c(1e-4, 10))
  } else {
    plot(df$time_years, df$Abeta_bisf, type = "l", col = "blue", lwd = 2,
         xlab = "Time (years)", ylab = "Concentration (nM)",
         main = "Aβ Monomer")
  }
  lines(df$time_years, df$Abeta_csf, col = "green", lwd = 2)
  lines(df$time_years, df$Abeta_plasma, col = "red", lwd = 2)
  legend("bottomright", c("Brain ISF", "CSF", "Plasma"),
         col = c("blue", "green", "red"), lwd = 2, cex = 0.8)

  # Add target lines
  abline(h = 0.75, col = "blue", lty = 2)
  abline(h = 2.8, col = "green", lty = 2)
  abline(h = 0.066, col = "red", lty = 2)

  # --- Plot 2: Aβ Oligomer ---
  if (log_y) {
    plot(df$time_years, df$Aolig_bisf, type = "l", col = "blue", lwd = 2,
         log = "y", xlab = "Time (years)", ylab = "Concentration (nM)",
         main = "Aβ Oligomer", ylim = c(1e-14, 10))
  } else {
    plot(df$time_years, df$Aolig_bisf, type = "l", col = "blue", lwd = 2,
         xlab = "Time (years)", ylab = "Concentration (nM)",
         main = "Aβ Oligomer")
  }
  lines(df$time_years, df$Aolig_csf, col = "green", lwd = 2)
  lines(df$time_years, df$Aolig_plasma, col = "red", lwd = 2)
  legend("bottomright", c("Brain ISF", "CSF", "Plasma"),
         col = c("blue", "green", "red"), lwd = 2, cex = 0.8)

  # Add target lines
  abline(h = 4.4, col = "blue", lty = 2)
  abline(h = 0.02, col = "green", lty = 2)
  abline(h = 6.9e-12, col = "red", lty = 2)

  # --- Plot 3: Aβ Plaque ---
  plot(df$time_years, df$Aplaq_bisf, type = "l", col = "purple", lwd = 2,
       xlab = "Time (years)", ylab = "Concentration (nM)",
       main = "Aβ Plaque (Brain ISF)")
  abline(h = 39, col = "purple", lty = 2)
  legend("bottomright", c("Simulated", "Target (39 nM)"),
         col = "purple", lty = c(1, 2), lwd = 2, cex = 0.8)

  # --- Plot 4: Total Aβ in CSF ---
  df$total_csf <- df$Abeta_csf + df$Aolig_csf
  plot(df$time_years, df$total_csf, type = "l", col = "darkgreen", lwd = 2,
       xlab = "Time (years)", ylab = "Concentration (nM)",
       main = "Total Aβ in CSF (for SILK validation)")
  abline(h = 2.8 + 0.02, col = "darkgreen", lty = 2)

  par(mfrow = c(1, 1))
}

# =============================================================================
# 5. Quick Run and Validate Function
# =============================================================================

#' Run steady state simulation and validate against Table S2
#'
#' @param t_years Simulation duration in years (default: 67)
#' @param plot Generate plots (default: TRUE)
#' @return List with simulation result and validation comparison
#'
run_and_validate <- function(t_years = 67, plot = TRUE) {

  cat("\n")
  cat("################################################################\n")
  cat("#           AD QSP Model - Steady State Validation             #\n")
  cat("################################################################\n\n")

  # Run simulation
  sim_result <- run_to_steady_state(t_years = t_years)

  # Validate against Table S2
  comparison <- validate_steady_state(sim_result)

  # Plot if requested
  if (plot) {
    plot_steady_state_approach(sim_result)
  }

  return(list(
    simulation = sim_result,
    validation = comparison
  ))
}

# =============================================================================
# 6. Check Rate of Plaque Change at Steady State
# =============================================================================

#' Check if plaque is at steady state (< 1% change per year)
#'
#' @param sim_result Result from run_to_steady_state()
#' @return Annual percent change in plaque
#'
check_plaque_steady_state <- function(sim_result) {

  df <- sim_result$output
  n <- nrow(df)

  # Get plaque values at last two time points
  t1 <- df$time[n-1]
  t2 <- df$time[n]
  p1 <- df$Aplaq_bisf[n-1]
  p2 <- df$Aplaq_bisf[n]

  # Calculate annual rate of change
  dt_years <- (t2 - t1) / SECONDS_PER_YEAR
  annual_pct_change <- ((p2 - p1) / p1) / dt_years * 100

  cat("\n=== Plaque Steady State Check ===\n")
  cat(sprintf("Plaque at t=%0.1f years: %.4f nM\n",
              seconds_to_years(t1), p1))
  cat(sprintf("Plaque at t=%0.1f years: %.4f nM\n",
              seconds_to_years(t2), p2))
  cat(sprintf("Annual rate of change: %.4f %%/year\n", annual_pct_change))

  if (abs(annual_pct_change) < 1) {
    cat("Status: STEADY STATE REACHED (< 1%/year change)\n")
  } else {
    cat("Status: NOT YET AT STEADY STATE\n")
  }

  return(annual_pct_change)
}

# =============================================================================
# 7. Save/Load Steady State
# =============================================================================

#' Save steady state for later use
#'
#' @param sim_result Result from run_to_steady_state()
#' @param filename File path to save
#'
save_steady_state <- function(sim_result, filename = "data/steady_state.rds") {
  saveRDS(sim_result, filename)
  cat("Steady state saved to:", filename, "\n")
}

#' Load previously saved steady state
#'
#' @param filename File path to load
#' @return Simulation result
#'
load_steady_state <- function(filename = "data/steady_state.rds") {
  if (!file.exists(filename)) {
    stop("Steady state file not found: ", filename)
  }
  sim_result <- readRDS(filename)
  cat("Steady state loaded from:", filename, "\n")
  return(sim_result)
}
