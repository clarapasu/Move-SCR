# ==============================================================================
# MARTEN POPULATION DENSITY ESTIMATION (SURVEY 2) USING UMOVE-SCR MODEL
# ==============================================================================
# Author: Clara Panchaud
# Description: Analysis of the second survey of marten camera trap data using Markov-Modulated 
#              Poisson Process (MMPP) to estimate population density. AcMove-SCR model used here to model uniform movement. 
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

compile("Functions/like_MMPP.cpp")
dyn.load(dynlib("Functions/like_MMPP"))


# Set random seed for reproducibility
set.seed(1)


# ------------------------------------------------------------------------------
#  LOAD AND PREPARE DATA
# ------------------------------------------------------------------------------

# Read the camera trap data and the trap locations
data <- read.csv("Data/marten_data2.csv")
traps <- read.csv("Data/marten_traps2.csv")


mesh_spacing <- 0.5
# Create a trap polygon and a state space mask
trap <- make.poly(x = traps$x, y = traps$y)
mask <- make.mask(trap, buffer = 2, spacing = mesh_spacing, type = "trapbuffer")

dim(mask)

# Convert the traps to the matrix format
traps <- as.matrix(traps)

# Get the number of states (spatial cells in mask)
S <- dim(mask)[1]
print(paste("Number of states (S):", S))


# ------------------------------------------------------------------------------
#  MAP TRAPS TO MASK LOCATIONS
# ------------------------------------------------------------------------------

# Calculate the distances between traps and mask points
distance_traps <- proxy::dist(traps, mask)

# Find the closest mask point for each trap
traps_loc <- apply(distance_traps, 1, which.min)

# Count the number of traps
n_traps <- 30

# Initialise the camera counts (one per mask cell)
camera_count <- rep(1, S)

# Check that at most one trap is located in each grid cell
trap_counts <- table(traps_loc)
print("Trap counts by location:")
print(trap_counts)


# Create neighbour matrix for spatial structure
neighbour <- neighbour_matrix(mask)

# Remap trap identifiers in data to mask locations
data$y <- traps_loc[data$y]

# Set observation time period
Time <- 11 

# ------------------------------------------------------------------------------
#  CALCULATE SPATIAL METRICS
# ------------------------------------------------------------------------------

# Cell spacing (km)
h <- attr(mask, "spacing")

# Area per cell (km²)
a <- h * h

# Total study area (km²)
A <- nrow(mask) * a

# ------------------------------------------------------------------------------
# PREPARE DATA FOR TMB MODEL
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
observed_ind <- length(mmpp$tt) - 1
total_detections <- sum(lengths(mmpp$tt))
# ------------------------------------------------------------------------------
# CONSTRUCT MODEL MATRICES
# ------------------------------------------------------------------------------

# Q matrix structure:
# - Fixed at 1 for non-neighbours (structural zeros)
# - Template has 1 for neighbours (to be filled with transition rates)
Q_fixed <- 1 - neighbour
Q_template <- neighbour

# Lambda (observation rate) structure:
# - Fixed at 1 for non-trap locations (no detections possible)
# - Free parameters (0 in fixed, 2 in template) at trap locations
lambda_fixed <- rep(1, S)
lambda_fixed[traps_loc] <- 0

lambda_template <- rep(0, S)
lambda_template[traps_loc] <- 2

# Initial state distribution (uniform across all states)
f <- rep(1/S, S)



# ------------------------------------------------------------------------------
# BUILD TMB DATA OBJECT
# ------------------------------------------------------------------------------

tmbdata <- list(
  Us = lengths(mmpp$tt),              # Number of detections per individual
  t = unlist(mmpp$tt),                # Detection times (flattened)
  s = as.integer(unlist(mmpp$ss) - 1L), # Detection locations (0-indexed for C++)
  f = f,                              # Initial state distribution
  n_obs = observed_ind,               # Number of observed individuals
  mesh_spacing = mesh_spacing,
  t0 = 0,                             # Start time
  Time = Time,                        # Survey length (end time)
  n_states = S,                       # Number of spatial states
  n_indiv = length(mmpp$tt),          # Total individuals (including the last empty one)
  Q_template = Q_template,            # Transition rate matrix template
  Q_fixed = Q_fixed,                  # Fixed elements of Q
  lambda_template = lambda_template,  # Observation rate template
  lambda_fixed = lambda_fixed,        # Fixed elements of lambda
  camera_count = camera_count        # In this case 1 trap max per cell
)

# Set length to zero for individuals with no capture history
tmbdata$Us[unlist(lapply(mmpp$tt, \(x) all(is.na(x))))] <- 0


tmbdata
# ------------------------------------------------------------------------------
#  MODEL FITTING
# ------------------------------------------------------------------------------

# Initial parameter values 
theta_init <- c(-1.2, -2.5)  # (lambda, alpha). lambda is on the log-scale

# Record the start time
start <- Sys.time()



# Create the TMB autodiff function
obj <- MakeADFun(
  data = tmbdata, 
  parameters = list(theta = theta_init), 
  DLL = "like_MMPP"
)



# Optimise the model (suppress convergence output)
invisible(capture.output({
  fit <- nlminb(
    start = theta_init, 
    objective = obj$fn, 
    gradient = obj$gr
  )
}))

# Record the end time and calculate duration of the optimisation
end <- Sys.time()
optimization_time <- difftime(end, start, units = "secs")

# Get the parameter estimates with standard errors
report <- sdreport(obj, par.fixed = fit$par)
param <- summary(report)


nll <- obj$fn(fit$par)


AIC <- 2*2 + nll*2 
print(AIC)

# ------------------------------------------------------------------------------
#  POPULATION ESTIMATION 
# ------------------------------------------------------------------------------

fit$par # alpha, lambda 

# Calculate the population size and confidence intervals
Pop <- confint_pop_mmmpp(fit$par, Time, mesh_spacing, observed_ind, S, neighbour, 
                              traps_loc, f, report$cov)


# Calculate the density estimates
D <- Pop[1] / A           # Density (animals/km²)
SE.D <- Pop[2] / A        # Standard error of density
CI_lower_D <- D - 1.96 * SE.D  # Lower 95% CI
CI_upper_D <- D + 1.96 * SE.D  # Upper 95% CI

CI_upper_D

# ------------------------------------------------------------------------------
# DISPLAY RESULTS
# ------------------------------------------------------------------------------

# Extract parameters 
param_2 <- param[3, "Estimate"]
param_2_SE <- param[3, "Std. Error"]
param_2_CI_lower <- param_2 - 1.96 * param_2_SE
param_2_CI_upper <- param_2 + 1.96 * param_2_SE

param_3 <- param[4, "Estimate"]
param_3_SE <- param[4, "Std. Error"]
param_3_CI_lower <- param_3 - 1.96 * param_3_SE
param_3_CI_upper <- param_3 + 1.96 * param_3_SE

# Create summary table (Table 4.1 of the thesis)
param_table <- data.frame(
  Parameter = c("Density (animals/km²)", "sigma^2", "lambda"),
  Estimate = c(round(D, 2), round(param_2, 2), round(param_3, 2)),
  SE = c(round(SE.D, 3), round(param_2_SE, 2), round(param_3_SE, 2)),
  CI_Lower = c(round(CI_lower_D, 2), round(param_2_CI_lower, 2), round(param_3_CI_lower, 2)),
  CI_Upper = c(round(CI_upper_D, 2), round(param_2_CI_upper, 2), round(param_3_CI_upper, 2)),
  stringsAsFactors = FALSE
)

#Print the results 
print(param_table, row.names = FALSE)

#   Parameter Estimate    SE CI_Lower CI_Upper
# Density (animals/km²)     0.17 0.038     0.09     0.24
# sigma^2     0.38 0.080     0.22     0.54
# lambda     3.81 0.610     2.61     5.00

optimization_time # 218.0094 secs

AIC # 293.2175
