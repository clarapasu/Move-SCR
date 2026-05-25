# ==============================================================================
# MARTEN POPULATION DENSITY ESTIMATION (SURVEY 2) USING OUMOVE-SCR MODEL
# ==============================================================================
# Author: Clara Panchaud
# Description: Analysis of the second survey of marten camera trap data using Markov-Modulated 
#              Poisson Process (MMPP) to estimate population density. OUMove-SCR model used here to model
#              the activity centre reverting movement. 
# ==============================================================================

# ------------------------------------------------------------------------------
#  SETUP AND DEPENDENCIES
# ------------------------------------------------------------------------------

# Load required libraries
library(ggplot2)
library(secr)
library(TMB)
library(expm)
library(tidyverse)
library(fields)
library(dplyr)
library(purrr)

# Source custom functions
Rcpp::sourceCpp("Functions/like_CTSCR.cpp")
source("Functions/Extra_CTSCR.r")
source("Functions/Extra_Func.r")

# Compile and load TMB model
compile("Functions/like_MMPP_ac.cpp")
dyn.load(dynlib("Functions/like_MMPP_ac"))

# Set random seed for reproducibility
set.seed(1)



# ------------------------------------------------------------------------------
#  LOAD AND PREPARE DATA
# ------------------------------------------------------------------------------

# Read the camera trap data and the trap locations
data <- read.csv("Data/marten_data2.csv")
traps <- read.csv("Data/marten_traps2.csv")

# Create a trap polygon
trap <- make.poly(x = traps$x, y = traps$y)

# ------------------------------------------------------------------------------
# CREATE DUAL MESH SYSTEM
# ------------------------------------------------------------------------------

# === DUAL MESH APPROACH ===
# Create two different meshes with different purposes:

# 1. FINE MESH: For state space (animal movement between locations)
#    This represents the spatial resolution of the biological process
state_spacing <- 0.5 # Fine spacing for detailed spatial process 

mask_states <- make.mask(trap, buffer = 2, spacing = state_spacing, type = "trapbuffer")



# Cell spacing (km)
h <- attr(mask_states, "spacing")

# Area per cell (km²)
a <- h * h

# Total study area (km²)
A <- nrow(mask_states) * a

# 2. COARSE MESH: For activity center integration (numerical integration only)
#    This reduces computational load while maintaining accuracy
#ac_spacing <- 1.5  # Coarser spacing to reduce computational burden 
ac_spacing <- 1.5 # doable with 1 but long
mask_ac <- make.mask(trap, buffer = 2, spacing = ac_spacing, type = "trapbuffer")

n_ac <- dim(mask_ac)[1]

traps <- as.matrix(traps)

# ------------------------------------------------------------------------------
# MAP TRAPS TO STATE SPACE
# ------------------------------------------------------------------------------

# === TRAP LOCATION MAPPING ===
# Map traps to the state space mesh (where animals are detected)
distance_traps <- proxy::dist(traps, mask_states)
traps_loc <- apply(distance_traps, 1, which.min)

# Get mesh dimensions
S <- dim(mask_states)[1]  # Number of states = state space mesh size
S_AC <- dim(mask_ac)[1]   # Number of AC integration points

# Remap trap identifiers in data to mask locations
data$y <- traps_loc[data$y]

# Set observation time period
Time <- 11

# Build neighbour matrix 
neighbour <- neighbour_matrix(mask_states)

# ------------------------------------------------------------------------------
# PREPARE CAPTURE HISTORY DATA
# ------------------------------------------------------------------------------

# Group the data by individual ID and sort by time
df_grouped <- data %>% 
  group_by(id) %>% 
  arrange(Time, .by_group = TRUE)

# Extract the capture histories (spatial locations and times) for each individual
mmpp <- list(
  ss = df_grouped %>% group_split() %>% map(~ .x$y),
  tt = df_grouped %>% group_split() %>% map(~ .x$Time)
)

# Add an NA at the end as a convention
mmpp$ss[[length(mmpp$ss) + 1]] <- NA
mmpp$tt[[length(mmpp$tt) + 1]] <- NA

# Number of observed individuals
observed_ind <- length(mmpp$ss) - 1

# ------------------------------------------------------------------------------
# CONSTRUCT MODEL MATRICES
# ------------------------------------------------------------------------------

# Lambda (observation rate) structure:
# - Fixed at 1 for non-trap locations (no detections possible)
# - Free parameters at trap locations
lambda_fixed <- rep(1, S)
lambda_fixed[traps_loc] <- 0


# Count cameras at each location (handle multiple traps in same cell)
camera_count <- tabulate(traps_loc, nbins = S)


# Q matrix structure:
# - Fixed at 1 for non-neighbours (structural zeros)
# - Template has 2 for neighbours (to be filled with transition rates)
Q_fixed <- 1 - neighbour

# ------------------------------------------------------------------------------
# BUILD TMB DATA OBJECT WITH DUAL MESH INFORMATION
# ------------------------------------------------------------------------------


tmbdata_ou <- list(
  Us = lengths(mmpp$tt),
  t = unlist(mmpp$tt),
  s = as.integer(unlist(mmpp$ss) - 1L),  # 0-based indexing for C++
  n_obs = observed_ind,
  mesh_spacing = state_spacing,
  t0 = 0,
  Time = Time,
  n_states = S,
  n_ac = S_AC,
  n_indiv = length(mmpp$tt),
  Q_fixed = Q_fixed,
  lambda_fixed = lambda_fixed,
  camera_count = camera_count,
  state_coords = as.matrix(mask_states),      
  ac_coords = as.matrix(mask_ac)  
)


# ------------------------------------------------------------------------------
#  MODEL FITTING
# ------------------------------------------------------------------------------

# Initial parameter values: (log_lambda, theta1, theta2)
theta_init <- c(0, 0, 0)

# Record start time
start <- Sys.time()

# Create the TMB autodiff function
obj <- MakeADFun(
  data = tmbdata_ou, 
  parameters = list(theta = theta_init), 
  DLL = "like_MMPP_ac", 
  silent = FALSE
)

# Optimise the model
out <- nlminb(theta_init, obj$fn, obj$gr)

# Record end time and calculate duration
end <- Sys.time()
optimization_time <- difftime(end, start, units = "hours")
# ------------------------------------------------------------------------------
# COMPUTE STANDARD ERRORS
# ------------------------------------------------------------------------------

# Get standard errors using sdreport
tryCatch({
  rep <- sdreport(obj, par.fixed = out$par)
}, error = function(e) {
  cat("Standard error computation failed:", e$message, "\n")
})
optimization_time

#---------------------------------------------------------------------------
# COMPUTE CONFIDENCE INTERVALS
# ------------------------------------------------------------------------------

Pop <- confint_pop_mmmpp_ac(out$par, Time, mask_states, 
                            state_spacing,
                            as.matrix(mask_ac),
                            observed_ind, 
                            S,
                            neighbour,
                            traps_loc,
                            rep$cov,
                            n_ac)



# ------------------------------------------------------------------------------
# DISPLAY RESULTS
# ------------------------------------------------------------------------------

# Extract parameter estimates with confidence intervals
param_estimates <- summary(rep)


param_1 <- param_estimates[4, "Estimate"]
param_1_SE <- param_estimates[4, "Std. Error"]
param_1_CI_lower <- param_1 - 1.96 * param_1_SE
param_1_CI_upper <- param_1 + 1.96 * param_1_SE

param_2 <- param_estimates[5, "Estimate"]
param_2_SE <- param_estimates[5, "Std. Error"]
param_2_CI_lower <- param_2 - 1.96 * param_2_SE
param_2_CI_upper <- param_2 + 1.96 * param_2_SE

param_3 <- param_estimates[6, "Estimate"]
param_3_SE <- param_estimates[6, "Std. Error"]
param_3_CI_lower <- param_3 - 1.96 * param_3_SE
param_3_CI_upper <- param_3 + 1.96 * param_3_SE



# Calculate confidence intervals for N 
CI_lower_N <- Pop[1] - 1.96 * Pop[2]      # Lower 95% CI
CI_upper_N <- Pop[1] + 1.96 * Pop[2]     # Upper 95% CI


# Create summary table
results_table <- data.frame(
  Parameter = c("Density (animals/km²)", "sigma^2", "alpha", "lambda"),
  Estimate = c(round(Pop[1], 4), round(param_1, 4), round(param_2, 4), round(param_3, 4)),
  SE = c(round(Pop[2], 4), round(param_1_SE, 4), round(param_2_SE, 4), round(param_3_SE, 4)),
  CI_Lower = c(round(CI_lower_N, 4), round(param_1_CI_lower, 4), 
               round(param_2_CI_lower, 4), round(param_3_CI_lower, 4)),
  CI_Upper = c(round(CI_upper_N, 4), round(param_1_CI_upper, 4), 
               round(param_2_CI_upper, 4), round(param_3_CI_upper, 4)),
  stringsAsFactors = FALSE
)



# Print the results
cat("=== FINAL RESULTS ===\n")
print(results_table, row.names = FALSE)

# Display optimisation results
cat("\n=== OPTIMISATION TIME ===\n")
cat("Time:", round(optimization_time, 2), "hours\n")


# Calculate AIC
nll <- obj$fn(out$par)
AIC <- 2 * 3 + nll * 2  

# Display AIC
cat("AIC:", AIC, "\n\n")
