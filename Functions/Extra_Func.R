# ==============================================================================
# POPULATION SIZE ESTIMATION UNDER THE Move-SCR MODEL
# By Clara Panchaud
# ==============================================================================

# A helper function to start: 

#' Create Neighbour Matrix for Mesh
#' Identifies neighbouring points in a spatial mesh based on minimum distances
#' 
#' @param mesh Matrix of spatial coordinates
#' @param tol Tolerance for considering points as neighbours (default: 1e-6)
#' @return Binary matrix where 1 indicates neighbouring relationship
neighbour_matrix <- function(mesh, tol = 1e-6) {
  
  # Calculate the pairwise distances between all mesh points
  dist_mat <- as.matrix(dist(mesh))
  
  # Set the diagonal to infinity to exclude self-distances
  diag(dist_mat) <- Inf
  
  # Find minimum distance for each point
  min_dist <- apply(dist_mat, 1, min)
  
  # Initialise neighbour matrix with zeros
  neighbour <- matrix(0, nrow = nrow(dist_mat), ncol = ncol(dist_mat))
  
  # Mark neighbours: points within tolerance of minimum distance
  for (i in 1:nrow(dist_mat)) {
    neighbour[i, dist_mat[i, ] - min_dist[i] < tol] <- 1
  }
  
  return(neighbour)
}



#' Form the transition matrix Q
#'
#' @param theta the two parameters included in Q
#' @param ac_dist a matrix of the distances between the activity centre and the grid cells,
#'                such that ac_dist[i,k] is the distance between the activity centre i
#'                and the centre of grid cell k
#'                
#' @return  three lists of length equal to the number of simulated individuals, where Q_list contains
#'          the Q matrices that are now individual specific, same for pi_list containing the initial state probabilities
#'          and s0 the initial state locations. 


make_Q<-function(theta,ac_dist,mesh,N,S){
  
  Q_list <- list()
  pi_list <- list()
  s0 <- numeric(N)
  
  neighbour<-neighbour_matrix(mesh)
  
  for (i in 1:N){
    Q<-neighbour
    for (j in 1:S){
      for (k in 1:S){
        if(Q[j,k]==1){
          Q[j,k]= exp( theta[1] - theta[2] * ac_dist[i,k]) 
        }
      }
    }
    diag(Q)<- - rowSums(Q)
    
    eig <- eigen(t(Q))
    pi <- Re(eig$vectors[, which.min(abs(eig$values))])
    pi <- abs(pi) / sum(abs(pi))  # Ensure positive and normalized
    
    # Safety check
    if(any(pi < 0) || any(!is.finite(pi))) {
      pi <- rep(1/S, S)  # Fallback to uniform distribution
    }
    
    initial_state <- sample(1:S, size = 1, prob = pi)
    
    Q_list[[i]] <- Q
    pi_list[[i]] <- pi
    s0[i] <- initial_state
  }
  
  return(list(Q_list,pi_list,s0))
}





#' Form the transition matrix Q for OU version 
#'
#' @param theta the two parameters included in Q
#' @param ac_dist a matrix of the distances between the activity centre and the grid cells,
#'                such that ac_dist[i,k] is the distance between the activity centre i
#'                and the centre of grid cell k
#'                
#' @return  three lists of length equal to the number of simulated individuals, where Q_list contains
#'          the Q matrices that are now individual specific, same for pi_list containing the initial state probabilities
#'          and s0 the initial state locations. 


make_Q_OU <- function(theta , ac , mesh , N){
    
    Q_list <- list()
    pi_list <- list()
    s0 <- numeric(N)
    
    S <- nrow(mesh)
    
    a <- attr(mesh, "area")* 10000
    mesh_spacing <- sqrt(a)
    
    # Build neighbour matrix
    neighbour<-neighbour_matrix(mesh)
    
    # Parameters
    sigma_sq <- exp(theta[1])  # Diffusion coefficient
    alpha <- exp(theta[2])      # Mean reversion strength
    base_rate <- sigma_sq / (2 * mesh_spacing^2)
    
    
    for (i in 1:N){
      
      mu_x <- ac[i,1]
      mu_y <- ac[i,2]
      
      Q <- matrix(0, nrow = S, ncol = S)
      
      for (j in 1:S){
        
        # Coordinates of current state j 
        x_j <- mesh[j,1]
        y_j <- mesh[j,2]
        
        for (k in 1:S){
          
          if (j!=k && neighbour[j,k] == 1) {
            
            # Get coordinates of neighbour (state k) from mesh
            x_k <- mesh[k, 1]
            y_k <- mesh[k, 2]
            
            # Jump direction vector
            jump_x <- x_k - x_j
            jump_y <- y_k - y_j
            
            # Displacement from activity center
            disp_x <- x_j - mu_x
            disp_y <- y_j - mu_y
            
            # Drift dot product: (r_j - mu) · (r_k - r_j)
            drift_dot <- disp_x * jump_x + disp_y * jump_y
            
            # OU process rate formula
            drift_component <- alpha * drift_dot / (2 * mesh_spacing)
            Q[j, k] <- base_rate - drift_component
            
            # Handle negative rates
            if (Q[j, k] < 0) {
              warning(paste("Negative rate at (", j, ",", k, ") for individual", i, 
                            ": Q =", Q[j, k]))
              Q[j, k] <- 1e-10
              
            }
            
          }
          
        }
      }
      
      diag(Q)<- - rowSums(Q)
      
      # Compute stationary distribution
      
      if (alpha < 1e-4) {  # Essentially zero
        pi <- rep(1/S, S)  # Uniform distribution
      } else {
        eig <- eigen(t(Q))
        pi <- Re(eig$vectors[, which.min(abs(eig$values))])
        pi <- abs(pi) / sum(abs(pi))
      }
    
      
      # Safety check
      if(any(pi < 0) || any(!is.finite(pi))) {
        pi <- rep(1/S, S)  # Fallback to uniform distribution
      }
      
      # Sample initial state
      initial_state <- sample(1:S, size = 1, prob = pi)
      
      Q_list[[i]] <- Q
      pi_list[[i]] <- pi
      s0[i] <- initial_state
    }
    
    return(list(Q_list,pi_list,s0))
  }







#' Estimate population size and detection probability under the UMove-SCR model
#'
#' Computes the estimated population size (N̂) and the overall detection 
#' probability (1 - P_unobserved) based on the model parameters and the survey setup.
#'
#' @param theta numeric(2). Model parameters on the log-scale: 
#'        c(q, log(lambda)).
#'         where exp(q)     : baseline transition rate between states
#                lambda     : detection rate in states with traps on
#' @param Time numeric(1). Total duration of the survey.
#' @param observed_ind integer(1). Number of individuals observed at least once.
#' @param S integer(1). Number of spatial states in the mesh.
#' @param neighbour S x S matrix. Adjacency matrix defining
#'        neighbouring relationships between spatial states.
#' @param traps_on integer vector. Indices of states where traps are active.
#' @param f numeric vector of length S. Initial probability distribution over states.
#'
#' @return Numeric vector of length 2:
#'         \itemize{
#'           \item N_est: estimated total population size
#'           \item obs_prob: overall probability of being observed at least once
#'         }
#'         
#'         
Pop_size_mmmpp <- function(theta, Time, mesh_spacing, observed_ind, S, neighbour, traps_on, f) {
  q_hat <- exp(theta[1]) # transition rate
  l_hat <- exp(theta[2]) # detection rate
  
  # Construct generator matrix for movement
  Q_est <- neighbour * q_hat / (2 * mesh_spacing^2)
  diag(Q_est) <- -rowSums(Q_est)
  
  
  # Construct detection rate matrix (lambda_i = 0 outside trap states)
  camera_count = tabulate(traps_on, nbins = S)
  lambda_est <- rep(0, S)
  lambda_est[traps_on] <- l_hat * camera_count[traps_on]
  Lambda_mat <- diag(lambda_est)
  
  # Compute the probability of remaining unobserved over the entire survey
  # via matrix exponential of (Q - Lambda)
  unobs_prob <- sum(f %*% expm(Time * (Q_est - Lambda_mat)))
  
  # Detection probability
  obs_prob <- 1 - unobs_prob
  
  # Population size estimated via Horvitz–Thompson-like estimator
  N_est <- observed_ind / obs_prob
  
  return(c(N_est, obs_prob))
}



#' Compute population size estimate and its standard error
#'
#' @param theta_est numeric(2). Estimated parameters.
#' @param Time numeric(1). Total survey duration.
#' @param observed_ind integer(1). Number of individuals observed in the data.
#' @param S integer(1). Number of spatial states.
#' @param neighbour S x S matrix. Adjacency matrix.
#' @param traps_on integer vector. Indices of states with traps.
#' @param f numeric vector of length S. Initial state distribution.
#' @param cov 2x2 covariance matrix of parameter estimates.
#' @param distribution character(1). "poisson" or "binomial", determines the 
#'        sampling variance component.
#'     
#' @return Numeric vector of length 2:
#'         \itemize{
#'           \item N_est: estimated total population size
#'           \item SE: standard error of N_est
#'         }
confint_pop_mmmpp <- function(theta_est, Time, mesh_spacing, observed_ind, S, neighbour, traps_on, f, cov,distribution = "binomial") {
  
  # Compute point estimates of population size and detection probability using the previous function
  Pop_size <- Pop_size_mmmpp(theta_est, Time, mesh_spacing, observed_ind, S, neighbour, traps_on, f)
  N_est <- Pop_size[1]
  obs_prob <- Pop_size[2]
  
  # Numerical gradient of N_est with respect to theta (for delta method)
  grad_vec <- numDeriv::grad(
    func = function(theta) Pop_size_mmmpp(theta, Time, mesh_spacing, observed_ind, S, neighbour, traps_on, f)[1],
    x = theta_est
  )
  
  # Sampling variance component: depends on assumed distribution
  s2 <- switch (tolower(distribution),
                poisson  = observed_ind/obs_prob^2,
                binomial = observed_ind*((1-obs_prob)/(obs_prob^2)))
  
  # Total variance
  Var <- as.numeric(t(grad_vec) %*% cov %*% grad_vec) + s2
  SE <- sqrt(Var)
  
  return(c(N_est=N_est, SE=SE))
}



# ==============================================================================
# Repeat with the ACMove-SCR model 
# ==============================================================================


stationary_dist <- function(Q) {
  S <- nrow(Q)
  A <- rbind(t(Q), rep(1, S))
  b <- c(rep(0, S), 1)
  pi <- qr.solve(A, b)
  pi <- pmax(pi, 0)
  pi / sum(pi)
}


#' Estimate population size and detection probability under the ACMove-SCR model
#'
#' Computes the estimated population size (N̂) and the overall detection 
#' probability (1 - P_unobserved) based on the model parameters and the survey setup.
#'
#' @param theta numeric(3). Model parameters on the log-scale: 
#'        c(sigma^2, alpha, lambda) .
#' @param Time numeric(1). Total duration of the survey.
#' @param observed_ind integer(1). Number of individuals observed at least once.
#' @param S integer(1). Number of spatial states in the mesh.
#' @param neighbour S x S matrix. Adjacency matrix defining
#'        neighbouring relationships between spatial states.
#' @param traps_on integer vector. Indices of states where traps are active.
#' @param n_ac numeric(1). Number of cells in the activity centre mesh used for integration. 
#' @param ac_to_states_dist N_ac x S matrix.  Matrix defining the distances between each potential activity centre locations and the state space. 
#' 
#' 
#' @return Numeric vector of length 2:
#'         \itemize{
#'           \item N_est: estimated total population size
#'           \item obs_prob: overall probability of being observed at least once
#'         }

Pop_size_full_ac <- function(theta, Time, mesh,  mesh_spacing, ac, observed_ind, S, neighbour, traps_on, n_ac ) {
  # Parameters
  sigma_sq <- exp(theta[1])  # Diffusion coefficient
  alpha <- exp(theta[2])      # Mean reversion strength
  base_rate <- sigma_sq / (2 * mesh_spacing^2)
  l <- exp(theta[3]) # Detection rate lambda
  
  
  # Construct detection rate matrix with the camera count multiplying lambda
  camera_count = tabulate(traps_on, nbins = S)
  lambda_est <- rep(0, S)
  lambda_est[traps_on] <- l * camera_count[traps_on]

  Lambda_mat <- diag(lambda_est)
  
  unobs_prob <- 0
  
  # Construct generator matrix for movement
  for (i in 1:n_ac) {
    
    
    mu_x <- ac[i,1]
    mu_y <- ac[i,2]
    
    Q <- matrix(0, nrow = S, ncol = S)
    
    for (j in 1:S){
      
      # Coordinates of current state j 
      x_j <- mesh[j,1]
      y_j <- mesh[j,2]
      
      for (k in 1:S){
        
        if (j!=k && neighbour[j,k] == 1) {
          
          # Get coordinates of neighbour (state k) from mesh
          x_k <- mesh[k, 1]
          y_k <- mesh[k, 2]
          
          # Jump direction vector
          jump_x <- x_k - x_j
          jump_y <- y_k - y_j
          
          # Displacement from activity center
          disp_x <- x_j - mu_x
          disp_y <- y_j - mu_y
          
          # Drift dot product: (r_j - mu) · (r_k - r_j)
          drift_dot <- disp_x * jump_x + disp_y * jump_y
          
          # OU process rate formula
          drift_component <- alpha * drift_dot / (2 * mesh_spacing)
          Q[j, k] <- base_rate - drift_component
          
          # Handle negative rates
          if (Q[j, k] < 0) {
            warning(paste("Negative rate at (", j, ",", k, ") for individual", i, 
                          ": Q =", Q[j, k]))
            Q[j, k] <- 1e-10
            
          }
          
        }
        
      }
    }
    
    diag(Q)<- - rowSums(Q)
    
    # Use the stationary distribution as the initial distribution
    #eig <- eigen(t(Q))
    #pi <- Re(eig$vectors[, which.min(abs(eig$values))])
    #pi <- pi / sum(pi)
    pi <- stationary_dist(Q)
    
    # Compute the probability of remaining unobserved over the entire survey
    # via matrix exponential of (Q - Lambda)
    unobs_prob <- unobs_prob + sum(pi %*% expm(Time * (Q - Lambda_mat)))
  }
  
  p_unobs = unobs_prob / n_ac
  p_obs   <- max(1 - p_unobs, 1e-12)  
  
  # Population size estimated via Horvitz–Thompson-like estimator
  N_est <- observed_ind / p_obs
  
  return(c(N_est, p_obs))  
}


#' This function returns ONLY the population estimate from the previous function, as a scalar
Pop_size_mmmpp_ac <- function(theta, Time,mesh,  mesh_spacing, ac, observed_ind, S, neighbour, traps_on, n_ac) {
  Pop_full <- Pop_size_full_ac(theta, Time, mesh,  mesh_spacing, ac, observed_ind, S, neighbour, traps_on, n_ac)
  return(Pop_full[1])  
}


#' Compute population size estimate and its standard error
#'
#' @param theta_est numeric(3). Estimated parameters.
#' @param Time numeric(1). Total survey duration.
#' @param observed_ind integer(1). Number of individuals observed in the data.
#' @param S integer(1). Number of spatial states.
#' @param neighbour S x S matrix. Adjacency matrix.
#' @param traps_on integer vector. Indices of states with traps.
#' @param cov 2x2 covariance matrix of parameter estimates.
#' @param n_ac numeric(1). Number of cells in the activity centre mesh used for integration. 
#' @param distribution character(1). "poisson" or "binomial", determines the 
#'        sampling variance component.
#'     
#' @return Numeric vector of length 2:
#'         \itemize{
#'           \item N_est: estimated total population size
#'           \item SE: standard error of N_est
#'         }
# Fixed confidence interval function
confint_pop_mmmpp_ac <- function(theta_est, Time, mesh,  mesh_spacing, ac, observed_ind, S, neighbour, traps_on, cov, n_ac, distribution = "binomial") {
  
  # Compute point estimates of population size and detection probability using the previous function
  Pop_full <- Pop_size_full_ac(theta_est, Time, mesh,  mesh_spacing, ac, observed_ind, S, neighbour, traps_on,   n_ac)
  N_est <- Pop_full[1]
  obs_prob <- Pop_full[2]
  
  # Calculate gradient using the scalar function
  grad_vec <- numDeriv::grad(
    func = function(theta) Pop_size_mmmpp_ac(theta, Time, mesh,  mesh_spacing, ac, observed_ind, S, neighbour, traps_on, n_ac),
    x = theta_est
  )
  
  # Sampling variance component: depends on assumed distribution
  s2 <- switch(tolower(distribution),
               poisson  = observed_ind/obs_prob^2,
               binomial = observed_ind*((1-obs_prob)/(obs_prob^2)))
  
  # Total variance
  Var <- as.numeric(t(grad_vec) %*% cov %*% grad_vec) + s2
  SE <- sqrt(Var)
  
  return(c(N_est=N_est, SE=SE))
}




# Function to extract and plot individual trajectories (black & white, filled traps)
plot_individual_trajectory <- function(sim_obj, individual_id, mask, traps_loc = NULL, 
                                          title = NULL) {
  
  # Extract trajectory for specific individual
  traj <- sim_obj$ss[[individual_id]]
  coords <- mask[traj, ]
  trajectory_data <- data.frame(
    x = coords$x,
    y = coords$y,
    state = traj,
    individual = individual_id
  )
  
  # Base plot
  p <- ggplot() +
    # Mask points as light background
    geom_point(data = mask, aes(x = x, y = y), color = "gray90", size = 0.6)+
  
    
    # Trajectory path in black
    geom_path(data = trajectory_data, aes(x = x, y = y),
              color = "black", size = 0.6) +
    
    # Traps (filled black diamonds)
    {
      if (!is.null(traps_loc)) {
        trap_coords <- mask[unique(traps_loc), ]
        geom_point(data = trap_coords, aes(x = x, y = y),
                   shape = 23, size = 2.5, fill = "black", color = "black")
      }
    } +
    
    labs(
      x = "Easting", y = "Northing"
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 11)
    ) +
    coord_equal()
  
  return(p)
}



create_coarse_grid <- function(fine_mesh, coarsening_factor) {
  # Extract fine grid properties
  fine_coords <- as.matrix(fine_mesh)
  fine_spacing <- attr(fine_mesh, "spacing")
  
  # New spacing for coarse grid
  coarse_spacing <- fine_spacing * coarsening_factor
  
  # Find the BOUNDING BOX of the fine grid (extent of grid cells, not just centers)
  x_min <- min(fine_coords[, 1]) - fine_spacing/2
  x_max <- max(fine_coords[, 1]) + fine_spacing/2
  y_min <- min(fine_coords[, 2]) - fine_spacing/2
  y_max <- max(fine_coords[, 2]) + fine_spacing/2
  
  # Calculate how many coarse cells fit in this bounding box
  n_coarse_x <- floor((x_max - x_min) / coarse_spacing)
  n_coarse_y <- floor((y_max - y_min) / coarse_spacing)
  
  # Center the coarse grid within the same bounding box
  # Start from the center of the first coarse cell
  x_start <- x_min + coarse_spacing/2
  y_start <- y_min + coarse_spacing/2
  
  # Create coarse grid centers
  x_seq <- x_start + (0:(n_coarse_x - 1)) * coarse_spacing
  y_seq <- y_start + (0:(n_coarse_y - 1)) * coarse_spacing
  
  coarse_coords <- expand.grid(x = x_seq, y = y_seq)
  
  # Convert to mask object
  coarse_mesh <- as.data.frame(coarse_coords)
  class(coarse_mesh) <- c("mask", "data.frame")
  attr(coarse_mesh, "spacing") <- coarse_spacing
  attr(coarse_mesh, "area") <- coarse_spacing^2 * 1e-4  # hectares
  attr(coarse_mesh, "boundingbox") <- c(x_min, x_max, y_min, y_max)
  
  return(coarse_mesh)
}


