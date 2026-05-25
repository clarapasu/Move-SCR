# ==============================================================================
# MARTEN POPULATION DENSITY ESTIMATION USING UMOVE-SCR MODEL
# ==============================================================================
# Author: Clara Panchaud
# Description: Analysis of the second survey of marten camera trap data using Markov-Modulated 
#              Poisson Process (MMPP) to estimate population density. RWMove-SCR model used here to model uniform movement. 
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

compile("Functions/like_MMPP.cpp")
dyn.load(dynlib("Functions/like_MMPP"))

# Set random seed for reproducibility
set.seed(1)


# ------------------------------------------------------------------------------
#  LOAD AND PREPARE DATA
# ------------------------------------------------------------------------------

# Read the camera trap data and the trap locations
data <- read.csv("Data/marten_data.csv")
traps <- read.csv("Data/marten_traps.csv")


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

# ------------------------------------------------------------------------------
#  POPULATION ESTIMATION 
# ------------------------------------------------------------------------------

# Calculate the population size and confidence intervals
Pop <- confint_pop_mmmpp(fit$par, Time, mesh_spacing, observed_ind, S, neighbour, 
                              traps_loc, f, report$cov)


# Calculate the population confidence intervals
CI_lower_N <- Pop[1] - 1.96 * Pop[2]  # Lower 95% CI
CI_upper_N <- Pop[1] + 1.96 * Pop[2] # Upper 95% CI

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


# Summary table
param_table <- data.frame(
  Parameter = c("Population size", "sigma^2", "lambda"),
  Estimate = c(round(Pop[1] , 2), round(param_2, 2), round(param_3, 2)),
  SE = c(round(Pop[2], 3), round(param_2_SE, 2), round(param_3_SE, 2)),
  CI_Lower = c(round(CI_lower_N, 2), round(param_2_CI_lower, 2), round(param_3_CI_lower, 2)),
  CI_Upper = c(round(CI_upper_N, 2), round(param_2_CI_upper, 2), round(param_3_CI_upper, 2)),
  stringsAsFactors = FALSE
)


#Print the results 
print(param_table, row.names = FALSE)

optimization_time 

AIC 
