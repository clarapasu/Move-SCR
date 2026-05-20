# ==============================================================================
# MARTEN POPULATION DENSITY ESTIMATION (SURVEY 2) USING ACMOVE-SCR MODEL
# ==============================================================================
# Author: Clara Panchaud
# Description: Analysis of the second survey of marten camera trap data using Markov-Modulated 
#              Poisson Process (MMPP) to estimate population density. AcMove-SCR model used here to model
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
source("MSCR/Fit_Func.r")
source("MSCR/Sim_Func.r")
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
n_ac
# ------------------------------------------------------------------------------
# CALCULATE SPATIAL DISTANCES
# ------------------------------------------------------------------------------

# === DISTANCE CALCULATIONS ===
#.          Use the same distance scaling for both meshes to maintain 
#           parameter interpretability across the dual mesh system

# Convert traps to matrix format
traps <- as.matrix(traps)

# Calculate raw distances for both meshes
#mesh_dist_raw_states <- rdist(mask_states, mask_states)
#mesh_dist_raw_ac <- rdist(mask_ac, mask_ac)

# Use the state space mesh to define the distance scale (the "biological" scale)
#mean_dist <- mean(mesh_dist_raw_states[mesh_dist_raw_states > 0])
#cat("Distance scale (from state mesh):", mean_dist, "\n")

# Scale both distance matrices using the same scaling factor
#mesh_dist_states <- mesh_dist_raw_states / mean_dist  # For spatial transitions
#mesh_dist_ac <- mesh_dist_raw_ac / mean_dist         # For AC integration

# Distance from AC mesh to state mesh (for mapping between meshes)
#ac_to_states_dist <- rdist(mask_ac, mask_states) / mean_dist

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




# Calculate density estimates (convert to per 100 ha)
D <- Pop[1] / A          # Density (animals/100 ha)
SE.D <- Pop[2] / A       # Standard error of density
CI_lower_D <- D - 1.96 * SE.D       # Lower 95% CI
CI_upper_D <- D + 1.96 * SE.D       # Upper 95% CI

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

# Create summary table
results_table <- data.frame(
  Parameter = c("Density (animals/km²)", "sigma^2", "alpha", "lambda"),
  Estimate = c(round(D, 4), round(param_1, 4), round(param_2, 4), round(param_3, 4)),
  SE = c(round(SE.D, 4), round(param_1_SE, 4), round(param_2_SE, 4), round(param_3_SE, 4)),
  CI_Lower = c(round(CI_lower_D, 4), round(param_1_CI_lower, 4), 
               round(param_2_CI_lower, 4), round(param_3_CI_lower, 4)),
  CI_Upper = c(round(CI_upper_D, 4), round(param_1_CI_upper, 4), 
               round(param_2_CI_upper, 4), round(param_3_CI_upper, 4)),
  stringsAsFactors = FALSE
)



# Print the results
cat("=== FINAL RESULTS ===\n")
print(results_table, row.names = FALSE)

# Display density-specific results
cat("\n=== DENSITY ESTIMATES ===\n")
cat("Density:", round(D, 4), "animals/100 ha\n")
cat("Standard Error:", round(SE.D, 4), "\n")
cat("95% CI: [", round(CI_lower_D, 4), ",", round(CI_upper_D, 4), "]\n")

# Display optimisation results
cat("\n=== OPTIMISATION TIME ===\n")
cat("Time:", round(optimization_time, 2), "hours\n")


# Calculate AIC
nll <- obj$fn(out$par)
AIC <- 2 * 3 + nll * 2  

# Display AIC
cat("AIC:", AIC, "\n\n")


# time: 2.16544 hours

#     Parameter Estimate     SE CI_Lower CI_Upper
#Density (animals/km²)   0.1621 0.0385   0.0867   0.2374
#sigma^2   0.4018 0.0902   0.2250   0.5786
#alpha   0.3671 0.4708  -0.5556   1.2899
#lambda   3.8550 0.6211   2.6377   5.0723

# AIC: 294.0837 
