# ==============================================================================
# SIMULATING MSCR AND CT SCR DATA
# By Clara Panchaud
# ==============================================================================

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

# Hazard function with an OU shape, for an individual 
#' @param k numeric(2). Coordinates of a location in the landscape (in the likelihood only evaluated at camera trap coordinates)
#' @param theta numeric(>=2). Model parameters (h0, sigma, beta) on the log-scale. Ignore beta if using CT SCR model
#' @param t numeric(1). Time since previous capture of the individual (days)
#' @param z numeric(2). Coordinates of previous location
#' @param s numeric(2). Activity centre coordinates 
#' @param m integer(1) Memory indicator (1 = MSCR, 0 = CT SCR)
#' @param logscale logical. Returns the function value on the log-scale if TRUE. 
#' 
#' @return Numeric scalar: the (log-)hazard function at location k
hazard<-function(k, theta, t, z, s, m, logscale = FALSE){
  kx <- k[1]
  ky <- k[2]
  h0<-exp(theta[1])   
  sigma2 <-exp(theta[2]) # represents sigma^2
  if (m==0){
    # the model reduces to CT SCR when beta tends to infinity, here we use -beta directly so setting 
    # beta to -exp(100) reduces to the CT SCR model
    beta<--exp(100)}
  if (m==1){
    beta <- -exp(theta[3])} # the negative sign is included in beta here and the exponential ensures the correct sign 
  B<-exp(beta*t) 
  denom <- sigma2 * (1 - B^2) #  = sigma2 if beta = -exp(100)
  denom <- pmax(denom, .Machine$double.eps)
  mu<-B*z+(1-B)*s # = s if beta = -exp(100)
  if(logscale){
    h <- log(h0) + ((-(1/2) * ( (kx - mu[1])^2 + (ky - mu[2])^2 )  / denom))
  } else {
    h <- h0 * exp(-(1/2) * ( (kx - mu[1])^2 + (ky - mu[2])^2 )  / denom)
  }
  return(h)
}

# Cumulative hazard function (sum of the hazard function over all traps) for an individual
#' @param trap Matrix with 2 columns. A row contains a camera trap's coordinates
#' @param theta numeric(>=2). Model parameters (h0, sigma, beta) on the log-scale. Ignore beta if using CT SCR model
#' @param t numeric(1). Time since previous capture of the individual (days)
#' @param z numeric(2). Coordinates of previous location
#' @param s numeric(2). Activity centre coordinates 
#' @param m integer(1) Memory indicator (1 = MSCR, 0 = CT SCR)
#' 
#' @return Numeric scalar: sum of OU hazard over all traps
total_hazard<-function(trap,theta,t,z,s,m){
  h.<-sum(apply(trap,1,hazard,theta=theta,t=t,z=unlist(z),s=unlist(s),m=m))
  return(h.)
}

# Survival function approximation between two times via midpoint hazard
#' @param trap Matrix with 2 columns. A row contains a camera trap's coordinates
#' @param theta numeric(>=2). Model parameters (h0, sigma, beta) on the log-scale. Ignore beta if using CT SCR model
#' @param t numeric vector of times (days), increasing from 0 to the end of the survey T
#' @param i integer. Current time index (i>j)
#' @param j integer. Index of most recent capture time
#' @param z numeric(2). Coordinates of previous location
#' @param s numeric(2). Activity centre coordinates 
#' @param m integer(1) Memory indicator (1 = MSCR, 0 = CT SCR)
#' 
#' @return Numeric scalar in (0,1]: approximation of the survival function
Survival<-function(trap,theta,t,i,j,z,s,m){
  # Midpoint rule for integrated hazard over (t[i-1], t[i])
  Surv<-exp(-(t[i]-t[i-1])*total_hazard(trap,theta, t[i-1]-t[j]+((t[i]-t[i-1])/2),z,s,m))
  return(Surv)
}


#' Simulates capture histories under MSCR or CT SCR
#' 
#' Simulates N individuals over time [0,T]. Only the observed individuals are returned
#'
#' @param N integer. True total population size
#' @param T_days numeric. Total survey time (days)
#' @param r integer. Number of intervals used in the time discretisation
#' @param theta numeric(>=2). Model parameters (h0, sigma, beta) on the log-scale. Ignore beta if using CT SCR model
#' @param trap matrix with 2 columns. A row contains a camera trap's coordinates
#' @param ac matrix of dimension N x 2. Contains the activity center coordinates for the N individuals
#' @param m integer(1) Memory indicator (1 = MSCR, 0 = CT SCR)
#' @param mesh matrix with 2 columns. Contains the coordinates of mesh points forming a discretisation of the landscape, used as potential starting locations
#' 
#' @return A data frame with columns:
#' id: The unique id of the individual corresponding to that observation
#' Time: The time of the observation (days)
#' y: The index of the trap leading to that observation
sim_data<-function(N,T_days,r,theta,trap,ac,m,mesh){
  data<-data.frame(matrix(ncol = 3, nrow = 0))
  names<-c("id","t","y")
  colnames(data)<-names
  time<-seq(0,T_days,T_days/r)
  for (n in 1:N){
    s<-ac[n,]
    y<-c()
    capt_time<-time[1] 
    start_hazards<-apply(mesh,1,halfnormal,theta=theta,s=s)
    z<-mesh[sample(1:nrow(mesh), size = 1, prob = start_hazards),]
    j<-1
    for (i in 2:length(time)){
      ti<-time[i]
      proba<-Survival(trap,theta,time,i,j,z,s,m)
      u<-runif(1)
      if (u<=proba){
        x<-0
      }
      else{
        traps_proba<-apply(trap,1,hazard,theta=theta,t=(ti-capt_time),z=z,s=s,m=m)
        p<-traps_proba/sum(traps_proba)
        x<-sample(length(traps_proba),1,prob=p)
        z<-trap[x,]
        capt_time<-ti
        j<-i
      }
      y<-c(y,x)
    }
    id<-rep(n,length(time)-1)
    Time<-time[-1]
    dat<-data.frame(id,Time,y)
    data<-rbind(data,dat)
  }
  df<-data
  for (j in 1:N){
    if (identical(df[df$id==j,]$y,rep(0,length(time)-1))){
      df<-df[!(df$id==j),]
    }
  }
  # keep only the actual detections, "0" means there are no captures at this step
  df<-df[!df$y==0,]
  rownames(df) <- NULL
  return(df)
}


##  Renumber IDs after removing unobserved individuals
#' @param df data frame from sim_data, containing only observed individuals but with original numbering of individual IDs
#' 
#' @return data frame with new IDs from 1 to n
re_id<-function(df){
  n<-length(unique(df$id)) 
  ids<-as.data.frame(table(df$id))
  ids$Var1<-c(1:n)
  ids <- ids[rep(ids$Var1,ids$Freq),1:(ncol(ids)-1)]
  df$id<-ids
  return(df)
}
