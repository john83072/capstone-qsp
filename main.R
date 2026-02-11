# =============================================================================
# main.R
# Main execution script for AD QSP Model
# Based on: Madrasi et al. (2021) Alzheimer's & Dementia
# =============================================================================

# Clear environment
rm(list = ls())

# Set working directory (adjust as needed)
# setwd("C:/Users/baeje/AD")

# Load required packages ------------------------------------------------------
# install.packages(c("deSolve", "ggplot2", "tidyverse"))
library(deSolve)
library(ggplot2)
library(tidyverse)

# Source model files ----------------------------------------------------------
source("R/01_parameters.R")
source("R/02_model_odes.R")
source("R/03_simulation.R")
source("R/04_drug_models.R")
source("R/05_validation.R")
source("R/06_analysis.R")

# Run Model -------------------------------------------------------------------

# 1. Calculate steady state (67 years without drug)

# 2. Validate against SILK data

# 3. Simulate drug treatments

# 4. Generate figures

# 5. Export results
