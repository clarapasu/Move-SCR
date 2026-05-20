# ==============================================================================
# Simulation of a MMMPP
# By Clara Panchaud (adapted from code by Paul Blackwell)
# ==============================================================================

#' Simulate a continuous-time Markov chain, given its generator Q
#'
#' @param Q  matrix of transition rates
#' @param Time  duration of simulation
#' @param s0  initial state(s)
#' @param N  number of realisations
#' @param t0  initial time
#'
#' @return  object specifying one or more Markov chain trajectories
#' @export
#'
#' @examples
#' 
sim_CTMC <- function(Q,Time,s0,N=1,t0=0)
{
  #
  # Housekeeping
  #
  S <- dim(Q)[1]
  stopifnot(s0%in%(1:S))
  stopifnot(t0<Time)
  if (length(s0)<N) s0 <- rep(s0,length=N)
  #
  # Derived parameters
  #
  r <- -diag(Q)
  W <- sweep(Q,1,r,"/")
  diag(W) <- 0
  #
  # Loop over individuals
  #
  n <- rep(NA,N) # empty vector
  ss <- tt <- vector("list",N+1) # empty list
  for (i in 1:N)
  {
    j <- 1
    current_t <- t0
    current_s <- s0[i]
    s <- c(current_s)
    t <- c(current_t)
    #
    # Loop over time
    #
    while(current_t<Time)
    {
      rate <- r[current_s]
      delta_t <- if (rate>0) rexp(1,rate) else Inf
      if(current_t + delta_t < Time)
      {
        # Switch
        j <- j+1
        current_t <- current_t + delta_t
        t[j] <- current_t
        current_s <- sample(1:S,prob=W[current_s,],size=1)
        s[j] <- current_s
      } else current_t <- Time
    }
    n[i] <- j
    ss[[i]] <- s
    tt[[i]] <- t
  }
  ss[[N+1]] <- tt[[N+1]] <- NA
  return(list(N=N,s0=s0,t0=t0,Time=Time,S=S, # always known
              Q=Q, # known because simulated
              n=n,ss=ss,tt=tt)) # stochastic
}




#' Simulate a Markov modulated Poisson process based on one realisation of a CTMC
#'
#' @param traj  realisation of CTMC
#' @param lambda  vector of observation rates
#'
#' @return  an object representing the MMPP
#' @export
#'
#' @examples
#' 
sim_MMPP_one <- function(traj,lambda)
  # traj is list(n,s0,s,t0,t,Time)
{
  times <- c(traj$t0,traj$t,traj$Time)
  durations <- diff(times)
  states <- c(traj$s0,traj$s)
  n_intervals <- traj$n+1
  s <- t <- c()
  for (j in 1:n_intervals)
  {
    state_here <- states[j]
    n_here <- rpois(1,durations[j]*lambda[state_here])
    if (n_here>0)
    {
      #cat("Observations in interval ",j,"when state is",state_here,"\n")
      times_here <- sort(runif(n_here,times[j],times[j+1]))
      t <- c(t,times_here)
      s <- c(s,rep(state_here,n_here))
    }
  }
  return(list(s=s,t=t))
}






#' Simulate a Markov modulated Poisson process based on a CTMC
#'
#' @param obj  object representing one of more realisations of a CTMC
#' @param lambda  vector of observation rates
#'
#' @return  an object representing the MMPP
#' @export  
#'
#' @examples
#' 
sim_MMPP <- function(obj,lambda)
{
   S <- obj$S
   N <- obj$N
   Time <- obj$Time
   t0 <- obj$t0
   n <- obj$n
   s0 <- obj$s0
   ss <- obj$ss
   tt <- obj$tt
   #
   # Loop over individuals
   #
   ss_out <- tt_out <- vector("list",N+1) # empty list
   for (i in 1:N)
   {
      #can(i)
      traj_i <- list(n=n[i],s0=s0[i],s=ss[[i]],t0=t0,t=tt[[i]],Time=Time)
      #print(traj_i)
      proc_i <- sim_MMPP_one(traj_i,lambda)
      #print(proc_i)
      ss_out[[i]] <- proc_i$s
      tt_out[[i]] <- proc_i$t
   }
   ss_out[[N+1]] <- tt_out[[N+1]] <- NA
   return(list(S=S,N=N,ss=ss_out,tt=tt_out,t0=t0,Time=Time))
}



#' Simulate a continuous-time Markov chain with ACMove and individual-specific Q matrices
#'
#' @param Q_list list of generator matrices for each individual
#' @param Time  duration of simulation
#' @param s0  initial state(s)
#' @param S  State space size
#' @param N  number of realisations (Individuals)
#' @param t0  initial time
#'
#' @return  object specifying one or more Markov chain trajectories
#' @export
#'
#' @examples
#' 
sim_CTMC_ac <- function(Q_list,Time,s0,S,N,t0=0)
{
  #
  # Housekeeping
  #
  stopifnot(s0%in%(1:S))
  stopifnot(t0<Time)
  if (length(s0)<N) s0 <- rep(s0,length=N)
  #
  # Loop over individuals
  #
  n <- rep(NA,N) # empty vector
  ss <- tt <- vector("list",N+1) # empty list
  for (i in 1:N)
  {
    
    Q <- Q_list[[i]]
    r <- -diag(Q)
    W <- sweep(Q,1,r,"/")
    diag(W) <- 0
    
    s <- c()
    t <- c()
    j <- 0
    current_t <- t0
    current_s <- s0[i]
    #
    # Loop over time
    #
    while(current_t<Time)
    {
      rate <- r[current_s]
      delta_t <- if (rate>0) rexp(1,rate) else Inf
      if(current_t + delta_t < Time)
      {
        # Switch
        j <- j+1
        current_t <- current_t + delta_t
        t[j] <- current_t
        current_s <- sample(1:S,prob=W[current_s,],size=1)
        s[j] <- current_s
      } else current_t <- Time
    }
    n[i] <- j
    ss[[i]] <- s
    tt[[i]] <- t
  }
  ss[[N+1]] <- tt[[N+1]] <- NA
  return(list(N=N,s0=s0,t0=t0,Time=Time,S=S, # always known
              n=n,ss=ss,tt=tt)) # stochastic
}




