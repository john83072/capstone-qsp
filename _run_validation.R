source("R/01_parameters.R")
source("R/02_model_odes.R")
source("R/03_simulation.R")

result <- run_and_validate(t_years = 67, plot = FALSE)

check_plaque_steady_state(result$simulation)
