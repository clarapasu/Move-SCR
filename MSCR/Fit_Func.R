# ==============================================================================
# FITTING THE MSCR AND CT SCR MODELS
# By Clara Panchaud
# ==============================================================================
# Conventions:
# - theta = c(log(h0), log(sigma^2), log(beta)) for MSCR; CT-SCR ignores beta.
# - Time units are days; trap, mesh, s, z are 2D coordinates (matrix with 2 columns).
# - memory indicator: 1 = MSCR, 0 = CT-SCR.
# - This file depends on halfnormal(), hazard(), Survival() defined in Sim_Func.R.
# - a data set consists of capture histories where each row contains; 
#     id: the individual's unique ID 
#     Time: the time of the capture (assuming that the survey started at time 0)
#     y: the index of the trap where the capture occurred 
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


#' Individual likelihood conditional on activity centre
#' 
#' Computes the likleihood contribution of one individual, given
#' their discretised capture history and their activity centre s.
#' 
#' @param theta numeric(>=2). Model parameters (h0, sigma, beta) on the log-scale. Ignore beta if using CT SCR model
#' @param trap Matrix with 2 columns. A row contains a camera trap's coordinates
#' @param data data frame obtained from discretize() for this individual
#' @param s numeric(2). Activity centre coordinates 
#' @param memory integer(1) Memory indicator (1 = MSCR, 0 = CT SCR)
#' 
#' @return Numeric scalar: likelihood value for the individual's capture history, conditional on its activity centre location
Likelihood_ind<-function(theta,trap,data,s,memory){
  m<-dim(data)[1]
  L<-1
  t<-data$t
  y<-data$y
  capt_time<-t[1]
  
  # Initial hazard (no prior capture): use the halfnormal distribution
  start_hazards<-sum(apply(trap,1,halfnormal,theta=theta,s=unlist(s)))
  
  # First capture time and location
  first_capt<-which(data$y!=0)[1]
  capt_time<-t[first_capt]
  z<-trap[y[first_capt],]
  
  # Likelihood = survival to first capture * hazard at first capture
  L<-exp(-capt_time*start_hazards)*halfnormal(z,theta,s)
  j<-first_capt
  
  # Subsequent intervals: survival between captures and hazard at captures. Now using the OU hazard
  for (i in (first_capt+1):m){
    L<-L*Survival(trap,theta,t,i,j,z,s,memory)
    if (y[i]!=0){
      k<-trap[y[i],]
      L<-L*hazard(k,theta,(t[i]-capt_time),z,s,memory)
      z<-k
      capt_time<-t[i]
      j<-i
    }
  }
  return(L)
}


#' Integrate individual likelihood over activity-centre locations
#'
#' @param theta numeric(>=2). Model parameters (h0, sigma, beta) on the log-scale. Ignore beta if using CT SCR model
#' @param trap Matrix with 2 columns. A row contains a camera trap's coordinates
#' @param data data frame obtained from discretize() for this individual
#' @param mesh mesh object (see secr package) with 2 columns, representing the potential activity centre locations for integration
#' @param memory integer(1) Memory indicator (1 = MSCR, 0 = CT SCR)
#' 
#' @return Numeric scalar: likelihood value for the individual's capture history, with activity centre integrated out
Likelihood_integrate<-function(theta,trap,data,mesh,memory){
  a = attr(mesh,"a") # grid cell area
  D<-dim(mesh)[1] # number of grid cells in the mesh
  A = a*D # total mesh area
  mesh<-matrix(c(mesh[,1],mesh[,2]),ncol=2,nrow=D)
  L<-0
  for (i in 1:D){
    L<-L+Likelihood_ind(theta,trap,data,mesh[i,],memory)
  }
  return(log(L*a/A))
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


#' Negative log-likelihood for the full dataset with AC integrated out
#'
#' @param theta numeric(>=2). Model parameters (h0, sigma, beta) on the log-scale. Ignore beta if using CT SCR model
#' @param trap Matrix with 2 columns. A row contains a camera trap's coordinates
#' @param df data frame of all detections with columns id, Time and y
#' @param mesh mesh object (see secr package) with 2 columns, representing the potential activity centre locations for integration
#' @param T_days numeric(1). Total survey duration (days)
#' @param r integer(1). Number of segments in the base discretisation
#' @param memory integer(1) Memory indicator (1 = MSCR, 0 = CT SCR)
#' 
#' @return Numeric scalar: negative log-likelihood.
Likelihood<-function(theta,trap,df,mesh,T_days,r,memory){
  n<-length(unique(df$id))
  L<-0
  for (i in 1:n){ 
    data<-discretize(df[df$id==i,],T_days,r)
    L<-L+Likelihood_integrate(theta,trap,data,mesh,memory)
  }
  L<-L-n*Seen_int(T_days,trap,mesh,theta)
  return(-L)
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


