# ==============================================================================
# Functions to prepare the data and help fit the CT-SCR model
# By Clara Panchaud
# ==============================================================================

#' Discretise one individual's capture history on a regular time grid
#' 
#' Splits [0,T] into equal segments and then further splits at each observed capture time
#' so that event times are on the grid.
#' 
#' @param df data frame with columns Time (days) and y (trap index). 
#' Only the rows corresponding to the capture history of a single indivdual should be passed here
#' @param T_days numeric(1). Total survey duration (days)
#' @param r integer(1). Number of segments in the base discretisation
#' 
#' @return data frame with columns:
#'         t - sorted times including 0, grid breaks, and event times
#'         y - trap index at time t, 0 indicates no capture
discretize<-function(df,T_days,r){
  t<-round(seq(0, T_days, T_days/r),digits=3)
  y<-rep(0,r+1)
  for (j in 1:length(df$Time)){
    time<-round(df$Time[j],digits=3)
    i<-findInterval(time, t)
    # Handle computational issues as specific cases
    if (round(time,digits=2)==0.3){ 
      y[i+1]<-df$y[j]
    }
    else if (round(time,digits=2)==0.6){ 
      y[i+1]<-df$y[j]
    }
    else if (round(time,digits=2)==0.7){  
      y[i+1]<-df$y[j]
    }
    else if(round(time,digits=3)==round(t[i],digits=3)){
      y[i]<-df$y[j]
    }
    else{
      # insert the event time as its own grid point between t[i] and t[i+1]
      t<-append(t,time,after=i)
      y<-append(y,df$y[j],after=i)
    }
  }
  data<-data.frame(t,y)
  return(data)
}



#' Convert MMPP simulation output to discretised detection data used for CT SCR
#'
#' Turns the MMMPP data into a discretized detection history format suitable 
#' for analysis with CT SCR. The function maps 
#' detections to trap locations, removes duplicate detections within a time 
#' threshold, and discretises continuous time data into intervals.
#'
#' @param mmpp List. MMPP simulation output containing:
#'        \itemize{
#'          \item ss: spatial location index of detections
#'          \item tt: detection times
#'        }
#' @param traps_on Integer vector. Indices of active trap locations.
#' @param T Numeric(1). Total survey duration.
#' @param r2 Numeric(1). Parameter passed to discretise function.
#' @param mesh Data frame. Mesh defining the state space of the MMMPP, here defining the trap locations with columns x and y.
#'
#' @return List of length 2:
#'         \itemize{
#'           \item ddfmat_sim: matrix of discretised detection data
#'           \item dfrows_sim: numeric vector with number of detections per individual
#' 
mmpp_to_df<-function(mmpp,traps_on,T,r2,mesh){
  # Convert MMPP simulation output to observation data frame
  # Map over spatial locations (ss) and times (tt) to create tibbles
  obs_df <- map2(
    .x = mmpp$ss,
    .y = mmpp$tt,
    .f = ~ tibble(y = .x, Time = .y)
  ) %>%
    # Combine all observations into single data frame with ID column for each individual
    bind_rows(.id = "id") %>%
    mutate(id = as.integer(id))
  
  # Convert to base R data frame
  obs_df<-as.data.frame(obs_df)
  
  # Remove last row as it is an artifact 
  obs_df<-obs_df[-nrow(obs_df),]
  
  # Add trap coordinates from mesh based on trap locations
  obs_df$trap_x <- mesh[obs_df$y,]$x
  obs_df$trap_y <- mesh[obs_df$y,]$y
  
  # Remap trap IDs to match the active traps
  obs_df$y<- match(obs_df$y, traps_on)
  
  # Remove duplicate detections that occur within 1 hour (0.04166667 days)
  # This filters out rapid re-detections at the same trap
  df_sim <- obs_df %>%
    arrange(id, Time) %>% 
    filter(!(id == lag(id) & (Time - lag(Time) < 0.04166667)))
  
  # Pre-allocate data frame for discretized observations
  n_rows <- sum(table(df_sim$id))
  ddf_sim <- data.frame(
    t = numeric(n_rows),
    y = integer(n_rows),
    id = integer(n_rows)
  )
  
  row_idx <- 1
  # Loop through each individual and discretise their detection history
  for(i in unique(df_sim$id)){
    # Discretise observations for individual i into time intervals using r2
    data <- discretize(df_sim[df_sim$id == i, ],T, r2)
    data$id <- i
    
    # Determine how many rows this individual contributes
    n_new <- nrow(data)
    
    # Insert discretised data into pre-allocated data frame
    ddf_sim[row_idx:(row_idx + n_new - 1), ] <- data
    
    # Update row index for next individual
    row_idx <- row_idx + n_new
    
    # Convert to matrix format
    ddfmat_sim = as.matrix(ddf_sim)
    
    # Count detections per individual
    dfrows_sim = as.numeric(table(ddf_sim$id))
  }
  
  return(list(ddfmat_sim, dfrows_sim))
}

#' Half-normal hazard at a trap for a given activity centre
#' 
#' @param k numeric(2). Coordinates of a camera trap 
#' @param theta numeric(>=2). Model parameters (h0, sigma) on the log-scale
#' @param s numeric(2). Activity centre coordinates 
#' @param logscale logical. Returns the function value on the log-scale if TRUE. 
#' 
#' @return Numeric scalar: the (log-)hazard at trap k
halfnormal<-function(k, theta, s, logscale = FALSE){
  kx <- k[1]
  ky <- k[2]
  h0 <- exp(theta[1])
  sigma2 <- exp(theta[2]) #represents a variance sigma^2
  h <- h0*exp(-( (kx - s[1])^2 + (ky - s[2])^2 ) / (2*sigma2))
  if(logscale){h <- log(h)}
  return(as.numeric(h))
}




#' Probability of being observed at least once, conditional on the activity centre location
#'
#' @param T_days numeric(1). Total survey duration (days)
#' @param trap Matrix with 2 columns. A row contains a camera trap's coordinates
#' @param s numeric(2). Activity centre coordinates 
#' @param theta numeric(>=2). Model parameters (h0, sigma, beta) on the log-scale. Ignore beta if using CT SCR model
#' 
#' @return Numeric scalar in [0,1] representing the probability that an individual with activity centre s is observed at least once in the survey
seen<-function(T_days,trap,s,theta){
  h<-sum(apply(trap,1,halfnormal,theta=theta,s=s))
  U<-1-exp(-T_days*h)
  return(U)
}

#' Integrate the previous function over activity centres
#' 
#' @param T_days numeric(1). Total survey duration (days) 
#' @param trap Matrix with 2 columns. A row contains a camera trap's coordinates
#' @param mesh mesh object (see secr package) with 2 columns, representing the potential activity centre locations for integration
#' @param theta numeric(>=2). Model parameters (h0, sigma, beta) on the log-scale. Ignore beta if using CT SCR model
#' 
#' @return Numeric scalar: probability of being observed at least once in the survey
Seen_int<-function(T_days,trap,mesh,theta){
  a = attr(mesh,"a")
  D<-dim(mesh)[1]
  A = a*D
  mesh<-matrix(c(mesh[,1],mesh[,2]),ncol=2,nrow=D)
  S<-0
  for (i in 1:D){
    S<-S+seen(T_days,trap,as.double(mesh[i,]),theta)
  }
  return(log(S*a/A))
}





#' Wald confidence intervals for model parameters from an optimiser fit
#'
#' @param fit Result of an optimiser (e.g. from using optim) with elements par and hessian
#' @return data.frame with columns value, upper and lower
#' 
#' Computes the confidence intervals on the original (log) parameter scale.
confint_param<-function(fit){
  theta_est<-fit$par
  Hessian<-fit$hessian
  cov<-solve(Hessian)
  se<-sqrt(diag(cov)) 
  upper<-fit$par+1.96*se
  lower<-fit$par-1.96*se
  interval<-data.frame(value=fit$par, upper=upper, lower=lower)
  return(interval)
}


#' Add (possibly asymmetric) confidence limits
#'
#' @param df data.frame with columns estimate and SE.estimate
#' @param alpha numeric(1). Significance level (default 0.05)
#' @param loginterval logical. If TRUE, compute asymmetric log-normal CIs. Appropriate when the variance is on the log scale and back-transformed. 
#' @param lowerbound numeric(1). Lower bound (default 0)
#' 
#' @return The same data frame df with added columns lcl and ucl.
add.cl <- function (df, alpha, loginterval, lowerbound = 0){
  z <- abs(qnorm(1 - alpha/2))
  if (loginterval) {
    delta <- df$estimate - lowerbound
    df$lcl <- delta/exp(z * sqrt(log(1 + (df$SE.estimate/delta)^2))) + lowerbound
    df$ucl <- delta * exp(z * sqrt(log(1 + (df$SE.estimate/delta)^2))) + lowerbound
  } else {
    df$lcl <- pmax(lowerbound, df$estimate - z * df$SE.estimate)
    df$ucl <- df$estimate + z * df$SE.estimate
  }
  df
}


#' Population size estimate and confidence interval
#'
#' @param fit Result of an optimiser (e.g. from using optim) with elements par and hessian
#' @param T_days numeric(1). Total survey duration (days)
#' @param trap Matrix with 2 columns. A row contains a camera trap's coordinates
#' @param mesh mesh object (see secr package) with 2 columns, representing the potential activity centre locations for integration
#' @param n integer(1). Number of observed individuals 
#' @param distribution character. "binomial" (default) or "poisson" for the count model of n. 
#' @param loginterval logical. if TRUE (default) return asymmetric CIs on N (log-normal style)
#' @param alpha numeric(1). Significance level (default 0.05)
#' 
#' @return data frame with the estimate of N, the standard error, lower confidence interval bound and upper confidence interval bound. 
confint_pop<-function(fit,T_days,trap,mesh,n,distribution = "binomial", loginterval = TRUE, alpha = 0.05){
  theta_est<-fit$par
  Hessian<-fit$hessian
  cov<-solve(Hessian)
  
  #P(seen at least once) under estimated parameters
  Seen_est<-exp(Seen_int(T_days,trap,mesh,theta_est)) 
  
  # Define N as a function of theta so that we can take the gradients
  Pop_size<-function(theta,T_days,trap,mesh,n){
    Seen_est_t<-exp(Seen_int(T_days,trap,mesh,theta)) 
    N_est<-n/Seen_est_t
    return(N_est)
  }
  
  N_est<-Pop_size(theta_est,T_days,trap,mesh,n)
  d <- numDeriv::grad(Pop_size,theta_est, T_days=T_days,trap=trap,mesh=mesh,n=n)
  
  s2 <- switch (tolower(distribution),
                poisson  = n/Seen_est^2,
                binomial = n*((1-Seen_est)/(Seen_est^2)))
  Var <- d%*% cov %*%d + s2
  temp <- data.frame(row.names = c('N'), estimate = N_est,  SE.estimate = sqrt(Var))
  temp <- add.cl(temp, alpha=alpha, loginterval)
  return(temp)
}



