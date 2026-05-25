################################################################################
# Simulation study from the OUMove-SCR model with parameter alpha = 0.5
# Author: Clara Panchaud
# 
################################################################################

# Load required libraries ------------------------------------------------------
library(ggplot2)
library(secr)
library(TMB)
library(expm)
library(tidyverse)


# Source custom functions and load TMB model ------------------
source("Functions/sim_MMPP.r")
source("Functions/Extra_Func.r")
source("Functions/Extra_CTSCR.r")


compile("Functions/like_MMPP.cpp")
dyn.load(dynlib("Functions/like_MMPP"))

compile("Functions/like_MMPP_ac.cpp")
dyn.load(dynlib("Functions/like_MMPP_ac"))

Rcpp::sourceCpp("Functions/like_CTSCR.cpp")

# Set random seed for reproducibility 
set.seed(12)


# Load and prepare trap data ---------------------------------------------------

# Define the number of traps 
n_traps<- 30 

# Define the number of states
S <- 100 

mesh_spacing <- 1

# Create the mesh as a 10x10 grid of 100 square cells at a distance 1km from each other 
traps = make.grid(sqrt(S), sqrt(S) , spacex = 1, detector = "multi")
mesh = make.mask(traps,buffer=0.5,type="traprect",spacing=mesh_spacing) 
a <-  attr(mesh, "a")

# Get the bounding box of the state space mesh
x_range <- range(mesh[,1])
y_range <- range(mesh[,2])

# Create an inner area in which cameras can be located
buffer <- 2
inner_states <- which(
  mesh[,1] >= x_range[1] + buffer & 
    mesh[,1] <= x_range[2] - buffer & 
    mesh[,2] >= y_range[1] + buffer & 
    mesh[,2] <= y_range[2] - buffer 
)



# Create a fine grid covering the same area
fine_spacing <- 0.1  # resolution of the mesh to draw activity centres from
x_seq <- seq(x_range[1], x_range[2], by = fine_spacing)
y_seq <- seq(y_range[1], y_range[2], by = fine_spacing)

mesh_ac <- expand.grid(x = x_seq, y = y_seq)
class(mesh_ac) <- c("mask", "data.frame")
attr(mesh_ac, "spacing") <- fine_spacing
attr(mesh_ac, "a") <- fine_spacing^2



# For the coarser AC integration mesh
coarse_spacing <- 1.0  # resolution of integration mesh
x_seq_coarse <- seq(x_range[1], x_range[2], by = coarse_spacing)
y_seq_coarse <- seq(y_range[1], y_range[2], by = coarse_spacing)

mesh_ac_fit <- expand.grid(x = x_seq_coarse, y = y_seq_coarse)
class(mesh_ac_fit) <- c("mask", "data.frame")
attr(mesh_ac_fit, "spacing") <- coarse_spacing
attr(mesh_ac_fit, "a") <- coarse_spacing^2


# Now set areas consistently
a_ac <- attr(mesh_ac_fit, "a") * 10000
n_ac <- nrow(mesh_ac_fit)
A <- a_ac * n_ac  # Total area matches the integration grid
area_ratio <- a_ac / A 


# Set simulation parameters ----------------------------------------------------

# Movement parameter 
sigma_sq <- 1


# Attraction to activity centre
alpha <- 0.5

theta <- c(log(sigma_sq), log(alpha)) 

# Detection rate parameter 
l <- 0.5 # natural scale


# Number of simulations
n_sim <- 100

# Study duration
Time <- 11 # days

# True population size for the simulation study
N <- 20

# Initialise results data frame
col_names <- c(
  "sim_id", "alpha_true", "beta_true", "lambda_true", "N_true", "Time", "S", "n_ac",
  "n_obs", "n_events",
  "lambda_hat", "alpha_hat", "beta_hat",
  "se_lambda", "se_alpha", "se_beta", "N_hat", "SE_N",
  "lambda_hat_ac", "alpha_hat_ac", "beta_hat_ac",
  "se_lambda_ac", "se_alpha_ac", "se_beta_ac", "N_hat_ac", "SE_N_ac",
  "h0", "sigma_sq", "se_h0", "se_sigma_sq", "N_SCR", "SE_N_SCR",
  "elapsed_sec", "elapsed_sec_ac", "elapsed_sec_SCR")
results <- as.data.frame(matrix(NA, nrow = n_sim, ncol = length(col_names)))
colnames(results) <- col_names


# Run simulations --------------------------------------------------------------


for (j in 1:n_sim) {
  
  cat("\n========================================\n")
  cat("Starting simulation", j, "of", n_sim, "\n")
  cat("========================================\n")
  
  set.seed(j+200)
  
  #Sample the trap locations
  traps_on<-sort(sample( inner_states, size=n_traps))
  trap <- mesh[traps_on,]
  
  # Convert traps to matrix format
  traps <- as.matrix(trap)
  
  # Count cameras at each location (handle multiple traps at same location)
  camera_count <- rep(1, S)
  
  # Set up detection rates at each location
  lambda <- rep(0, S)
  lambda[traps_on] <- l * camera_count[traps_on]
  
  # Sample N activity centres
  ac<-mesh_ac[sample(1:nrow(mesh_ac), N, replace = TRUE), ]
  
  pre_data<-make_Q_OU(theta, ac, mask, N )
  
  Q_list<-pre_data[[1]]
  pi_list<-pre_data[[2]]
  s0<-pre_data[[3]]
  
  
  # Simulate continuous-time Markov chain (movement paths)
  obj <- sim_CTMC_ac(Q_list,Time,s0,S,N)
  
  # Simulate detections using Markov-modulated Poisson process
  mmpp <- sim_MMPP(obj, lambda)
  
  # plot_trajectory(obj, 4, mesh, trap)
  
  
  # Filter to observed individuals only (those with at least one detection)
  observed_idx <- which(!vapply(mmpp$tt, is.null, logical(1)))
  
  # Subset detection times and locations
  mmpp$tt <- mmpp$tt[observed_idx]
  mmpp$ss <- mmpp$ss[observed_idx]
  
  # Count observed individuals
  observed_ind <- length(mmpp$tt) - 1
  
  # Total number of detections
  total_detections <- sum(lengths(mmpp$tt))
  
  cat("  Observed individuals:", observed_ind, "\n")
  cat("  Total detections:", total_detections, "\n")
  
  
  
  ############## FIT WITH UMOVE FIRST ############
  ################################################
  
  cat("\n--- Fitting UMOVE model ---\n")
  
  # compute distances between mask and potential activity centres
  ac_fit_dist <- fields::rdist(mesh_ac_fit, mesh)
  
  neighbour<-neighbour_matrix(mask)
  
  # Construct Q and lambda templates
  Q_fixed <- 1 - neighbour         # Positions of structural zeros
  Q_template <- neighbour          # Positions to fill with transition rates
  
  
  lambda_fixed <- rep(1, S)
  lambda_fixed[traps_on] <- 0     # Free parameters at active trap locations only
  lambda_template <- rep(0, S)
  lambda_template[traps_on] <- 2
  
  
  # Initial distribution (uniform across all states)
  f <- rep(1/S, S)
  
  # Build TMB data list
  tmbdata <- list(
    Us = lengths(mmpp$tt),
    t = unlist(mmpp$tt),
    s = as.integer(unlist(mmpp$ss) - 1L),  # 0-based indexing for C++
    f = f,
    n_obs = observed_ind,
    mesh_spacing = mesh_spacing,
    t0 = 0,
    Time = Time,
    n_states = S,
    n_indiv = length(mmpp$tt),
    Q_template = Q_template,
    Q_fixed = Q_fixed,
    lambda_template = lambda_template,
    lambda_fixed = lambda_fixed,
    camera_count = camera_count
  )
  
  # Handle individuals with no detection history (set length to zero)
  tmbdata$Us[unlist(lapply(mmpp$tt, \(x) all(is.na(x))))] <- 0
  
  
  # Fit MMPP model -------------------------------------------------------------
  
  # Initial parameter values
  theta_init <- c( -4.5, -2.5) 
  
  # Record start time
  start_umove <- Sys.time()
  
  # Create autodiff function
  obj_umove <- MakeADFun(data = tmbdata, parameters = list(theta = theta_init), 
                         DLL = "like_MMPP")
  
  # Perform optimization (suppress output)
  invisible(capture.output({
    fit_umove <- nlminb(start = theta_init, objective = obj_umove$fn, gradient = obj_umove$gr)
  }))
  
  # Get parameter estimates and standard errors
  report_umove <- sdreport(obj_umove, par.fixed = fit_umove$par)
  param_umove <- summary(report_umove)
  
  
  
  # Record end time  
  end_umove <- Sys.time()
  elapsed_umove <- as.numeric(difftime(end_umove, start_umove, units = "secs"))
  
  
  # Calculate population size estimate with confidence interval
  Pop_size_umove <- confint_pop_mmmpp(fit_umove$par, Time, mesh_spacing, observed_ind, S, neighbour, 
                                      traps_on, f, report_umove$cov)
  
  
  cat("  UMOVE N estimate:", Pop_size_umove[1], "± SE:", Pop_size_umove[2], "\n")
  cat("  UMOVE elapsed time:", round(elapsed_umove, 2), "seconds\n")
  
  ##################################### FIT WITH OU MODEL #####################################
  #############################################################################################
  
  cat("\n--- Fitting OU (activity center) model ---\n")
  
  
  tmbdata_ou <- list(
    Us = lengths(mmpp$tt),
    t = unlist(mmpp$tt),
    s = as.integer(unlist(mmpp$ss) - 1L),  # 0-based indexing for C++
    n_obs = observed_ind,
    mesh_spacing = mesh_spacing,
    t0 = 0,
    Time = Time,
    n_states = S,
    n_ac = n_ac,
    n_indiv = length(mmpp$tt),
    Q_fixed = Q_fixed,
    lambda_fixed = lambda_fixed,
    camera_count = camera_count,
    state_coords = as.matrix(mesh),      
    ac_coords = as.matrix(mesh_ac_fit)  
  )
  
  
  tmbdata_ou$Us[unlist(lapply(mmpp$tt, \(x) all(is.na(x))))] <- 0
  
  
  # Initial parameters: log(sigma^2), log(theta), log(lambda)
  theta_init_ou <- c( 0.5, -2, -3.5)
  
  # Record start time 
  start_ou <- Sys.time()
  
  # Fit
  obj_ou <- MakeADFun(data = tmbdata_ou, parameters = list(theta = theta_init_ou), 
                      DLL = "like_MMPP_ac")
  
  
  fit_ou <- nlminb(start = theta_init_ou, objective = obj_ou$fn, gradient = obj_ou$gr)
  
  
  report_ou <- sdreport(obj_ou, par.fixed = fit_ou$par)
  param_ou <- summary(report_ou)
  
  # Record end time
  end_ou <- Sys.time()
  elapsed_ou <- as.numeric(difftime(end_ou, start_ou, units = "secs"))
  
  
  Pop_size_ac <- confint_pop_mmmpp_ac(fit_ou$par, Time, mesh, mesh_spacing, as.matrix(mesh_ac_fit), observed_ind, S, neighbour,
                                      traps_on, report_ou$cov, n_ac)
  
  
  cat("  OU N estimate:", Pop_size_ac[1], "± SE:", Pop_size_ac[2], "\n")
  cat("  OU elapsed time:", round(elapsed_ou, 2), "seconds\n")
  
  
  
  
  ##################################### FIT WITH CT SCR MODEL #####################################
  #############################################################################################
  
  r2<-50
  cams<-mesh[traps_on,]
  mesh_fit = make.mask(cams,buffer=2,type="traprect",spacing=0.5)
  
  
  new_format <- mmpp_to_df(mmpp,traps_on,Time,r2, mesh)
  ddfmat_sim <- new_format[[1]]
  dfrows_sim <- new_format[[2]]
  
  # fit CT SCR
  theta_init<-c(1,0.01)
  start_time_SCR <- Sys.time()
  fit_nomem <- optim(theta_init[1:2], LikelihoodCnoMem, trap = traps ,
                     df = ddfmat_sim, dfrows =dfrows_sim, mesh = as.matrix(mesh_fit), endt = Time, hessian=TRUE)
  end_time_SCR <- Sys.time()
  
  param_SCR <- fit_nomem$par
  
  N_est_SCR <- confint_pop(fit_nomem,Time,cams,mesh_ac_fit,observed_ind)
  SE_SCR <- sqrt(diag(solve(fit_nomem$hessian)))
  
  elapsed_SCR <- as.numeric(difftime(end_time_SCR, start_time_SCR, units="secs"))
  
  cat("  SCR N estimate:", as.numeric(N_est_SCR[1]), "± SE:", as.numeric(N_est_SCR[2]), "\n")
  cat("  SCR elapsed time:", round(as.numeric(elapsed_SCR), 2), "seconds\n")
  
  
  # Store results --------------------------------------------------------------
  
  results[j, "sim_id"] <- j
  results[j, "alpha_true"] <- sigma_sq
  results[j, "beta_true"] <- alpha
  results[j, "lambda_true"] <- l
  results[j, "N_true"] <- N
  results[j, "Time"] <- Time
  results[j, "S"] <- S
  results[j, "n_ac"] <- n_ac
  results[j, "n_obs"] <- observed_ind
  results[j, "n_events"] <- total_detections
  
  
  
  # UMOVE results
  results[j, "lambda_hat"] <- param_umove[3, 1]
  results[j, "alpha_hat"] <- param_umove[4, 1]
  results[j, "se_lambda"] <- param_umove[3, 2]
  results[j, "se_alpha"] <- param_umove[4, 2]
  results[j, "N_hat"] <- Pop_size_umove[1]
  results[j, "SE_N"] <- Pop_size_umove[2]
  results[j, "elapsed_sec"] <- elapsed_umove
  
  # OU results
  results[j, "alpha_hat_ac"] <- param_ou[4, 1]
  results[j, "beta_hat_ac"] <- param_ou[5, 1]
  results[j, "lambda_hat_ac"] <- param_ou[6, 1]
  results[j, "se_alpha_ac"] <- param_ou[4, 2]
  results[j, "se_beta_ac"] <- param_ou[5, 2]
  results[j, "se_lambda_ac"] <- param_ou[6, 2]
  results[j, "N_hat_ac"] <- Pop_size_ac[1]
  results[j, "SE_N_ac"] <- Pop_size_ac[2]
  results[j, "elapsed_sec_ac"] <- elapsed_ou
  
  # SCR results
  
  
  results[j, "h0"] <-  exp(param_SCR[1])
  results[j, "sigma_sq"] <-  exp(param_SCR[2])
  results[j, "se_h0"] <- exp(param_SCR[1]) * SE_SCR[1]
  results[j, "se_sigma_sq"] <- exp(param_SCR[2]) * SE_SCR[2]
  results[j, "N_SCR"] <- N_est_SCR[1]
  results[j, "SE_N_SCR"] <- N_est_SCR[2]
  results[j, "elapsed_sec_SCR"] <- elapsed_SCR
  results
  cat("\nSimulation", j, "complete!\n")
}


# Print summary of results -----------------------------------------------------
cat("\n========================================\n")
cat("All simulations complete!\n")
cat("========================================\n\n")


# Save results to CSV
write.csv(results, "simulation_results_1.csv", row.names = FALSE)


# Print summary statistics
cat("\n--- Summary Statistics ---\n")
cat("True N:", N, "\n")
cat("UMOVE mean N estimate:", round(mean(results$N_hat, na.rm = TRUE), 2), "\n")
cat("UMOVE mean SE:", round(mean(results$SE_N, na.rm = TRUE), 2), "\n")
cat("OU mean N estimate:", round(mean(results$N_hat_ac, na.rm = TRUE), 2), "\n")
cat("OU mean SE:", round(mean(results$SE_N_ac, na.rm = TRUE), 2), "\n")
cat("CT mean N estimate:", round(mean(results$N_SCR, na.rm = TRUE), 2), "\n")
cat("CT mean SE:", round(mean(results$SE_N_SCR, na.rm = TRUE), 2), "\n")

boxplot(results$N_hat,
        results$N_hat_ac,
        results$N_SCR,
        names = c("Category 1", "Category 2", "Category 3"),
        col = c("skyblue", "salmon", "lightgreen"),
        main = "Boxplot of Three Categories",
        ylab = "Values")


sqrt(mean((results$N_true - results$N_hat)^2))


sqrt(mean((results$N_true - results$N_SCR)^2))



mean(results$N_true - results$N_hat)
mean(results$N_true - results$N_SCR)



