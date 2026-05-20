// ==============================================================================
// FITTING THE MSCR MODEL
// By Clara Panchaud & Ian Durbach
// ==============================================================================
// C++ implementation of the likelihood functions similar as in Fit_Func.R, but through TMB to increase computational speed. Contains functions for:
// 1. Memory Spatial Capture-Recapture (MSCR) - incorporates memory effects
// 2. Continuous Time SCR without memory (CT SCR) - no memory component
//
// These functions are called from R using Rcpp and are used for maximum
// likelihood estimation of model parameters.
// ==============================================================================

#include <Rcpp.h>
using namespace Rcpp;

// ==============================================================================
// HAZARD FUNCTION FROM MSCR
// ==============================================================================
/**
 * Calculate the OU-inspired hazard rate from MSCR
 * 
 * The hazard incorporates memory by using an exponentially weighted average
 * of the last capture location and the activity centre, inspired by the distribution
 * of an OU process. The more recent is the previous capture, the more it impacts the hazard.
 * 
 * @param k0, k1: Location coordinates (x, y), corresponding to a trap location in the likelihood
 * @param theta: Parameter vector [h0, log(sigma^2), log(beta)]
 *               - h0: baseline hazard rate (log scale)
 *               - sigma^2: spatial scale parameter (log scale)  
 *               - beta: memory  parameter (log scale)
 * @param t: Time since previous capture
 * @param z0, z1: Coordinates of previous capture location
 * @param s0, s1: Activity centre coordinates
 * 
 * @return: Hazard rate at the location (k0,k1)
 */
double hazardC (
    double k0,         // Location x-coordinate
    double k1,         // Location y-coordinate
    NumericVector theta,       // Parameter vector [log(h0), log(sigma^2), log(beta)]
    double t,          // Time since last capture
    double z0,         // Last capture x-coordinate
    double z1,         // Last capture y-coordinate
    double s0,         // Activity center x-coordinate
    double s1          // Activity center y-coordinate
)
{
  // Memory decay factor: B = exp(-beta * t)
  // As t increases, B approaches 0 (memory effect fades)
  // As beta increases, memory decays faster
  double B = exp(-exp(theta[2]) * t);
  
  // The expected location is the weighted average of the last capture and the activity center
  // mu = B * last_capture + (1-B) * activity_center
  // When B=1 (t=0): expected location = last capture location
  // When B=0 (t=∞): expected location = activity centre
  double mu0 = B * z0 + (1-B) * s0;
  double mu1 = B * z1 + (1-B) * s1;
  
  double effective_variance = exp(theta[1]) - B * B * exp(theta[1]);
  double distance_squared = (k0 - mu0) * (k0 - mu0) + (k1 - mu1) * (k1 - mu1);
  double hz = exp(theta[0]) * exp(-0.5 * distance_squared / effective_variance);
  return hz;
}

// ==============================================================================
// STANDARD HALF-NORMAL DETECTION FUNCTION
// ==============================================================================
/**
 * Calculate half-normal hazard function
 * 
 * Standard continuous-time spatial capture-recapture hazard function based on distance
 * between trap and activity center. No memory effects.
 * 
 * @param k0, k1: Location coordinates (only evaluated at traps in the likelihood)
 * @param theta: Parameter vector [log(h0), log(sigma^2)]
 * @param s0, s1: Activity centre coordinates
 * 
 * @return: Hazard function at location (k0,k1)
 */
double halfnormalC (
    double k0,         // Location x-coordinate
    double k1,         // Location y-coordinate
    NumericVector theta,        // Parameter vector
    double s0,         // Activity center x-coordinate
    double s1          // Activity center y-coordinate
)
{
  double hn = exp(theta[0]) * exp(-0.5 * ((k0 - s0) * (k0 - s0)  + (k1 - s1) * (k1 - s1))  / exp(theta[1]));
  return hn;
}


// ==============================================================================
// MAIN LIKELIHOOD FUNCTION WITH MEMORY (MSCR)
// ==============================================================================
/**
 * Calculate negative log-likelihood for MSCR
 * 
 * This function implements the full likelihood for the MSCR model
 * 
 * @param theta: Parameter vector [log(h0), log(sigma^2), log(beta)]
 * @param trap: Matrix of trap coordinates (nTraps × 2)
 * @param df: Detection data matrix (nDetections × 3) [time, trap_id, individual_id]
 * @param dfrows: Vector indicating number of detections per individual
 * @param mesh: Spatial mesh for integration (nMeshPoints × 2)
 * @param endt: End time of survey period
 * @return: Negative log-likelihood value
 */
// [[Rcpp::export]]
double LikelihoodC(
    NumericVector theta,    // Parameter vector
    NumericMatrix trap,     // Trap coordinates
    NumericMatrix df,       // Data of capture histories
    IntegerVector dfrows,   // Detections per individual
    NumericMatrix mesh,     // Spatial mesh for integration
    double endt) {          // Survey end time
  
  // Get number of observed individuals and mesh dimensions
  int n = unique(df(_, 2)).length();  // Number of observed individuals
  int D = mesh.nrow();                // Number of mesh points
  int nT = trap.nrow();               // Number of traps
  
  
  // Variables for likelihood calculation
  double S;        // Probability of being detected at least once in the study (any trap, any time)
  double U;        // Detection probability for single mesh point
  double LS;       // Log probability of being detected
  double LL;       // Log-likelihood accumulator
  
  
  // ===========================================================================
  // PART 1: CALCULATE PROBABILITY OF BEING SEEN
  // ===========================================================================
  // For each possible activity center, calculate probability of being detected
  // at least once during the survey period. 
  
  
  NumericVector start_hazards(D);
  
  S = 0;
  for (int i = 0; i < D; i++) {
    start_hazards[i] = 0;
    for (int j = 0; j < nT; j++) {
      // Sum hazard rates across all traps (using half-normal, no memory initially)
      start_hazards[i] += halfnormalC(trap(j, 0), trap(j, 1), theta, mesh(i, 0), mesh(i, 1));
    }
    // Probability of detection = 1 - exp(-total_hazard * survey_duration)
    U = 1 - exp(-endt * start_hazards[i]);
    S += U;
  }
  LS = log(S);  // Log probability of detection (up to a constant)
  
  // ===========================================================================
  // PART 2: CALCULATE LIKELIHOOD FOR EACH INDIVIDUAL'S CAPTURE HISTORY
  // ===========================================================================
  
  LL = 0; // Initialise log-likelihood
  
  // Variables for processing individual capture histories
  int data_startrow;
  int data_endrow;
  int m;                // Number of detections for current individual
  NumericMatrix data;   // Detection data for current individual
  NumericVector t;      // Detection times
  NumericVector y;      // Trap IDs (y[i] > 0 means detection, y[i] = 0 means no detection)
  double capt_time;     // Time of most recent capture
  
  // Initialise data indices
  data_startrow = 0;
  data_endrow = dfrows[0] - 1;
  
  // Loop through each captured individual
  for (int hh = 0; hh < n; hh++) {
    
    // Extract detection data for current individual
    data = df(Range(data_startrow, data_endrow) , Range(0, 1));
    m = data.nrow();
    t = data(_, 0);
    y = data(_, 1);
    
    // Variables for likelihood calculation for this individual
    int seen_ind;        // Index of most recent capture
    double z0, z1;       // Coordinates of last capture location
    double k0, k1;       // Current trap coordinates
    double s0, s1;       // Activity center coordinates (unknown, integrated over)
    double total_hazard; // Total hazard rate across all traps
    double survival;     // Survival probability (no detection)
    double L;            // Likelihood for this individual
    double Lj;           // Likelihood for specific activity center j
    
    
    // Find first actual detection (skip any leading zeros)
    int row_ind = 0;
    int fcy = y[row_ind];
    while (fcy == 0) {
      row_ind += 1;
      fcy = y[row_ind];
    }
    
    // Initialise likelihood for this individual
    L = 0;
    
    // =========================================================================
    // INTEGRATE OVER ALL POSSIBLE ACTIVITY CENTERS
    // =========================================================================
    
    for (int j = 0; j < D; j++) {
      
      // Current activity center coordinates
      s0 = mesh(j, 0);
      s1 = mesh(j, 1);
      
      // Time and location of first capture
      capt_time = t[row_ind]; 
      k0 = trap(fcy - 1, 0); // Trap coordinates (convert to 0-based indexing)
      k1 = trap(fcy - 1, 1); 
      
      // Likelihood starts with:
      // 1. Survival from start to first capture: exp(-capt_time * total_hazard)
      // 2. Detection probability at first capture location
      Lj = exp(-(capt_time * start_hazards[j])) * halfnormalC(k0, k1, theta, s0, s1);
      seen_ind = row_ind;
      
      // Update memory: most recent capture location becomes reference point
      z0 = k0;
      z1 = k1;
      
      // =====================================================================
      // PROCESS REMAINING DETECTION HISTORY
      // =====================================================================
      
      
      for (int i = (row_ind + 1); i < m; i++) {
        
        // Calculate total hazard rate across all traps
        // Uses memory-based hazard function with time since last capture
        total_hazard = 0;
        
        for (int jj = 0; jj < nT; jj++) {
          // Time argument: time since last capture + half of current interval (for midpoint rule)
          // This accounts for continuous-time process within discrete intervals
          total_hazard += hazardC(trap(jj, 0), trap(jj, 1), theta, t[i-1] - t[seen_ind] + 0.5 * (t[i] - t[i-1]), z0, z1, s0, s1);
        }
        
        // Survival probability for this time interval
        survival = exp(-(t[i]-t[i-1]) * total_hazard);
        Lj *= survival;
        
        // Check if there was a detection at time t[i]
        if (y[i] > 0.5){ // Detection occurred
          // Get trap coordinates for this detection
          k0 = trap(y[i]-1, 0); 
          k1 = trap(y[i]-1, 1); 
          
          // Add detection to likelihood
          Lj *= hazardC(k0, k1, theta, (t[i] - capt_time), z0, z1, s0, s1);
          
          // Update memory: new capture becomes reference point
          z0 = k0;
          z1 = k1;
          capt_time = t[i];
          seen_ind = i;
        }
        // If y[i] = 0, no detection occurred (already accounted for in survival)
      }
      
      // Add contribution from this activity center to total likelihood
      L += Lj; 
    }
    
    // Add log-likelihood for this individual
    LL += log(L);
    
    // Update data indices for next individual
    data_startrow += dfrows[hh];
    data_endrow += dfrows[hh+1];
    
  }
  
  // ===========================================================================
  // FINAL LIKELIHOOD CALCULATION
  // ===========================================================================
  
  // Complete likelihood = product of individual likelihoods / (detection probability)^n
  // In log scale: log(L) = sum(log(L_i)) - n * log(S)
  LL = LL - n * LS;
  
  // Return negative log-likelihood (for minimisation)
  return -LL;
}





// ==============================================================================
// MAIN LIKELIHOOD FUNCTION WITHOUT MEMORY (CT SCR)
// ==============================================================================
/**
 * Calculate negative log-likelihood for CT SCR)
 * 
 * This function implements the full likelihood for the CT SCR model
 * 
 * @param theta: Parameter vector [log(h0), log(sigma^2)], beta is omitted 
 * @param trap: Matrix of trap coordinates (nTraps × 2)
 * @param df: Detection data matrix (nDetections × 3) [time, trap_id, individual_id]
 * @param dfrows: Vector indicating number of detections per individual
 * @param mesh: Spatial mesh for integration (nMeshPoints × 2)
 * @param endt: End time of survey period
 * @return: Negative log-likelihood value
 */
// [[Rcpp::export]]
double LikelihoodCnoMem(
    NumericVector theta,    // Parameter vector
    NumericMatrix trap,     // Trap coordinates
    NumericMatrix df,       // Data of capture histories
    IntegerVector dfrows,   // Detections per individual
    NumericMatrix mesh,     // Spatial mesh for integration
    double endt) {          // Survey end time
  
  // Get number of observed individuals and mesh dimensions 
  int n = unique(df(_, 2)).length(); // Number of observed inidividuals
  int D = mesh.nrow(); // Number of mesh points
  int nT = trap.nrow(); // Number of traps
  
  // Variables for likelihood calculation
  double S;     // Probability of being detected at least once in the survey
  double U;     // Detection probability for a single mesh point
  double LS;    // Log probability of being detected
  double LL;    // Log-likelihood accumulator
  
  // ===========================================================================
  // PART 1: CALCULATE PROBABILITY OF BEING SEEN
  // ===========================================================================
  // For each possible activity center, calculate probability of being detected
  // at least once during the survey period. 
  
  NumericVector start_hazards(D);
  
  S = 0;
  for (int i = 0; i < D; i++) {
    start_hazards[i] = 0;
    for (int j = 0; j < nT; j++) {
      // Sum hazard rates across all traps (CT SCR halfnormal hazard)
      start_hazards[i] += halfnormalC(trap(j, 0), trap(j, 1), theta, mesh(i, 0), mesh(i, 1));
    }
    // Probability of detection = 1 - exp(-total_hazard * survey_duration)
    U = 1 - exp(-endt * start_hazards[i]);
    S += U;
  }
  LS = log(S); // Log probability of detection (up to a constant)
  
  // ===========================================================================
  // PART 2: CALCULATE LIKELIHOOD FOR EACH INDIVIDUAL'S CAPTURE HISTORY
  // ===========================================================================
  
  
  LL = 0; // Initialise log-likelihood
  
  // Variables for processing individual capture histories
  int data_startrow;
  int data_endrow;
  int m;                // Number of detections for current individual
  NumericMatrix data;   // Detection data for current individual
  NumericVector t;      // Detection times
  NumericVector y;      // Trap IDs (y[i] > 0 means detection, y[i] = 0 means no detection)
  double capt_time;     // Time of most recent capture
  
  // Initialise data indices
  data_startrow = 0;
  data_endrow = dfrows[0] - 1;
  
  // Loop through each captured individuals
  for (int hh = 0; hh < n; hh++) {
    
    // Extract detection data for current individual
    data = df(Range(data_startrow, data_endrow) , Range(0, 1));
    m = data.nrow();
    t = data(_, 0);
    y = data(_, 1);
    
    // Variables for likelihood calculation for this individual
    int seen_ind;        // Index of most recent capture
    double z0, z1;       // Coordinates of last capture location
    double k0, k1;       // Current trap coordinates
    double s0, s1;       // Activity center coordinates (unknown, integrated over)
    double total_hazard; // Total hazard rate across all traps
    double survival;     // Survival probability (no detection)
    double L;            // Likelihood for this individual
    double Lj;           // Likelihood for specific activity center j
    
    // Find first actual detection (skip any leading zeros)
    int row_ind = 0;
    int fcy = y[row_ind];
    while (fcy == 0) {
      row_ind += 1;
      fcy = y[row_ind];
    }
    
    // Initialise likelihood for this individual
    L = 0;
    
    // =========================================================================
    // INTEGRATE OVER ALL POSSIBLE ACTIVITY CENTERS
    // =========================================================================
    
    for (int j = 0; j < D; j++) {
      
      // Current activity center coordinates
      s0 = mesh(j, 0);
      s1 = mesh(j, 1);
      
      // Time and location of first capture
      capt_time = t[row_ind]; 
      k0 = trap(fcy - 1, 0); 
      k1 = trap(fcy - 1, 1); 
      
      // Likelihood starts with:
      // 1. Survival from start to first capture: exp(-capt_time * total_hazard)
      // 2. Detection probability at first capture location
      Lj = exp(-(capt_time * start_hazards[j])) * halfnormalC(k0, k1, theta, s0, s1);
      seen_ind = row_ind;
      
      // =====================================================================
      // PROCESS REMAINING DETECTION HISTORY
      // =====================================================================
      
      for (int i = (row_ind + 1); i < m; i++) {
        
        // Calculate total hazard rate across all traps (using halfnormal hazard)
        total_hazard = 0;
        
        for (int jj = 0; jj < nT; jj++) {
          total_hazard += halfnormalC(trap(jj, 0), trap(jj, 1), theta, s0, s1);
        }
        
        // Survival probability for this time interval
        survival = exp(-(t[i]-t[i-1]) * total_hazard);
        Lj *= survival;
        
        // Check if there was a detection at time t[i]
        if (y[i] > 0.5){ // Detection occurred
          // Get trap coordinates for this detection
          k0 = trap(y[i]-1, 0); 
          k1 = trap(y[i]-1, 1); 
          
          // Add detection to likelihood
          Lj *= halfnormalC(k0, k1, theta, s0, s1);
          seen_ind = i;
        }
        // If y[i] = 0, no detection occurred (already accounted for in survival)
      }
      // Add contribution from this activity center to total likelihood
      L += Lj; 
      
    }
    
    // Add log-likelihood for this individual
    LL += log(L);
    
    // Update data indices for next individual
    data_startrow += dfrows[hh];
    data_endrow += dfrows[hh+1];
    
  }
  
  // ===========================================================================
  // FINAL LIKELIHOOD CALCULATION
  // ===========================================================================
  
  // Complete likelihood = product of individual likelihoods / (detection probability)^n
  // In log scale: log(L) = sum(log(L_i)) - n * log(S)
  LL = LL - n * LS;
  
  // Return negative log-likelihood (for minimisation)
  return -LL;
}
