// ==============================================================================
// LIKELIHOOD FOR UMove-SCR model
// By Clara Panchaud (adapted from code by Paul Blackwell and C++ version by David L. Miller)
// ==============================================================================
// Description:
// This file implements the continuous-time UMove-SCR likelihood for capture–recapture
// data, using the Template Model Builder (TMB) framework for efficient
// parameter estimation.
//
// The MMPP represents individual movement as a continuous-time Markov process
// with transition rate matrix Q, and detection events as a Poisson process with
// state-dependent intensity λ (lambda). The likelihood integrates over all
// unobserved state transitions between detections.
// 
// The code loops over individuals, computing each individual’s likelihood
// based on their capture times and corresponding states. 
// ==============================================================================


#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() () {
  
  // --------------------------------------------------------------------------
  // PARAMETER AND DATA INPUTS
  // --------------------------------------------------------------------------
  PARAMETER_VECTOR(theta); // log-transformed model parameters
  
  DATA_IVECTOR(Us);  // vector giving the number of observations per individual
  DATA_VECTOR(t);    // concatenated vector of capture times
  DATA_IVECTOR(s);   // corresponding vector of observed states (indices)
  DATA_VECTOR(f);    // initial state distribution (probabilities)
  
  
  DATA_SCALAR(t0);   // survey start time
  DATA_SCALAR(Time); // survey end time
  

  DATA_INTEGER(n_states); // number of spatial/latent states
  DATA_INTEGER(n_indiv);  // number of individuals (includes an extra empty capture history at the end)
  DATA_INTEGER(n_obs);    // number of individuals observed at least once
  DATA_SCALAR(mesh_spacing);  // spacing between cells of the mesh 
  
  DATA_IVECTOR(Q_template);         // template for transition parameters
  DATA_VECTOR(Q_fixed);             // indicator: 1 = fixed (no parameter), 0 = estimated
  DATA_IVECTOR(lambda_template);    // template for λ (detection) parameters
  DATA_VECTOR(lambda_fixed);        // indicator: 1 = fixed, 0 = estimated
  DATA_VECTOR(camera_count);        // vector of length S, at cells containing at least a trap
                                    // it contains the number of traps in that cell. Ones elsewhere. 
  
  // --------------------------------------------------------------------------
  // INITIALISE OBJECTS
  // --------------------------------------------------------------------------
  Type ll = 0.0;                           // negative log-likelihood accumulator
  vector<Type> etheta = exp(theta);        // parameters on natural scale
  
  // Model components
  vector<Type> lambda(n_states);           // detection rates
  matrix<Type> Lambda(n_states, n_states); // diagonal matrix of λ
  matrix<Type> Q(n_states, n_states);      // transition rate matrix
  matrix<Type> A(n_states, n_states);      // temporary matrix (for exponentials)
  
  // --------------------------------------------------------------------------
  // BUILD TRANSITION RATE MATRIX Q
  // --------------------------------------------------------------------------
  int Qit = 0;
  for(int j=0; j<n_states; j++){
    for(int k=0; k<n_states; k++){
      if(Q_fixed(Qit) == 1){
        Q(j,k) = 0.0; // fixed element (not estimated)
      }else{
        Q(j,k) = etheta(Q_template(Qit)-1) / (2.0 * mesh_spacing * mesh_spacing); // estimated transition rate, uniform movement. 
      }
      Qit++;
    }
  }
  // Diagonal elements ensure rows sum to zero
  Q.diagonal() << -1.0 * Q.rowwise().sum();
  
  // --------------------------------------------------------------------------
  // BUILD DETECTION RATE MATRIX Λ
  // --------------------------------------------------------------------------
  for(int j=0; j<n_states; j++){
    if(lambda_fixed(j) == 1){
      lambda(j) = 0.0;  // inactive state (no detection)
    }else{
      lambda(j) = etheta(lambda_template(j)-1) *camera_count(j); // estimated detection rate
    }
  }
  Lambda.setZero();
  Lambda.diagonal() << lambda;
  
  // Combined matrix representing movement when no observations occur
  matrix<Type>G = Q-Lambda;
  
  // --------------------------------------------------------------------------
  // INDIVIDUAL LIKELIHOODS
  // --------------------------------------------------------------------------
  int istart = 0; // running index to extract data for each individual
  Type ill;       // individual likelihood contribution
  
  for(int i=0; i<n_indiv; i++){
    
    // Number of capture events for individual i
    int U = Us(i);
    
    // extract capture history for this individual
    vector<Type> ti = t.segment(istart, U);
    vector<int> si = s.segment(istart, U);
    istart = istart + U;
    
    // --------------------------------------------------
    // Case 1: Individual not observed (U = 0)
    // This case will most often not be used, but we leave the option. 
    // --------------------------------------------------
    if (U==0){
      // Compute probability of remaining unobserved over full survey duration
      A = atomic::expm(matrix<Type>((Time-t0)*G));
      ill = (f.matrix().transpose()*A).sum();
    }
    // --------------------------------------------------
    // Case 2: Individual observed (U > 0)
    // --------------------------------------------------
    else{
      // Include start (t0) and end (Time) times
      vector<Type> ts(ti.size()+2);
      ts << t0, ti, Time;
      
      // Durations between successive capture or boundary times
      vector<Type> tau = diff(ts);
      vector<Type> w(U+1);
      
      // First interval (state before first detection unknown)
      A = atomic::expm(matrix<Type>(tau(0)*G));
      
      vector<Type> omega1 = A.col(si(0));
      w(0) = (omega1.array()*f.array()).sum();
      w(0) = w(0) * lambda(si(0)); // multiply by detection rate at the first capture 
      
      // Intermediate intervals (between detections)
      if(U>1){
        for(int u=1; u<U; u++){
          A = atomic::expm(matrix<Type>(tau(u)*G));
          w(u) = A(si(u-1), si(u)) * lambda(si(u));
        }
      }
      
      // Final interval (after last detection until end)
      A = atomic::expm(matrix<Type>(tau(U)*G));
      vector<Type> omega_Uplus1 = A.row(si(U-1));
      w(U) = sum(omega_Uplus1);
      
      // Product of all interval contributions
      ill = w.prod();
    }
    
    // Add log-likelihood contribution
    ll -= log(ill);
  } 
  
  // --------------------------------------------------------------------------
  // UNOBSERVED INDIVIDUALS CONTRIBUTION
  // --------------------------------------------------------------------------
  // Compute overall probability of being unobserved
  vector<Type> ones(n_states);
  ones.setOnes();
  A = atomic::expm(matrix<Type>(Time * G));
  Type p_unobs = (f.matrix().transpose() * A * ones.matrix())(0,0);
  Type p_obs = 1.0 - p_unobs;
  
  // Add term for the n_obs observed individuals
  ll += n_obs * log(p_obs);
  
  // --------------------------------------------------------------------------
  // REPORT PARAMETERS ON NATURAL SCALE and their standard errors
  // --------------------------------------------------------------------------
  ADREPORT(etheta); // avoids explicit delta-method transformations
  
  // Return total negative log-likelihood
  return(ll);
} 
