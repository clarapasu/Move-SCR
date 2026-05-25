// ==============================================================================
// LIKELIHOOD FOR Move-SCR model
// By Clara Panchaud (adapted from code by Paul Blackwell and C++ version by David L. Miller)
// ==============================================================================


// Optimizations added for performance
#include <TMB.hpp>
#include <iostream>
#include <vector>

// OPTIMIZED: More efficient stationary distribution using improved power iteration
template<class Type>
vector<Type> StationaryDist(matrix<Type>& Q) {
  
  int S = Q.rows();
  
  // Initialize with uniform distribution
  vector<Type> pi(S);
  for (int i = 0; i < S; i++) {
    pi[i] = Type(1.0) / Type(S);
  }
  
  // Find maximum diagonal rate more efficiently
  Type max_rate = Type(0);
  for (int i = 0; i < S; i++) {
    Type rate = -Q(i, i);
    if (rate > max_rate) max_rate = rate;
  }
  max_rate = max_rate + Type(1e-10);
  
  // Pre-compute P = I + Q/max_rate (avoid recomputing in loop)
  matrix<Type> P(S, S);
  for (int i = 0; i < S; i++) {
    for (int j = 0; j < S; j++) {
      P(i, j) = (i == j) ? Type(1) + Q(i, j) / max_rate : Q(i, j) / max_rate;
    }
  }
  
  // Power iteration with convergence check
  Type tolerance = Type(1e-12);
  for (int iter = 0; iter < 50; iter++) {  // Reduced max iterations
    vector<Type> new_pi(S);
    
    // Manual matrix-vector multiplication to avoid Eigen issues
    for (int j = 0; j < S; j++) {
      new_pi[j] = Type(0);
      for (int i = 0; i < S; i++) {
        new_pi[j] += pi[i] * P(i, j);
      }
    }
    
    // Normalize
    Type sum_pi = Type(0);
    for (int i = 0; i < S; i++) {
      sum_pi += new_pi[i];
    }
    if (sum_pi < Type(1e-10)) sum_pi = Type(1);
    
    for (int i = 0; i < S; i++) {
      new_pi[i] /= sum_pi;
    }
    
    // Check convergence
    Type diff = Type(0);
    for (int i = 0; i < S; i++) {
      Type delta = new_pi[i] - pi[i];
      diff += delta * delta;
    }
    diff = sqrt(diff);
    
    pi = new_pi;
    if (diff < tolerance) break;
  }
  
  return pi;
}



template<class Type> 
Type ind_likelihood(int i, int n_states, const vector<Type>& lambda,
                    const matrix<Type>& Q,
                    const vector<Type>& pi,
                    int U, const vector<Type>& t, 
                    const vector<int>& s, int& istart, Type t0, Type Time,
                    const vector<Type>& theta, int state_idx) {
  
  matrix<Type> G = Q;
  // Create G = Q - Lambda by modifying diagonal elements
  for(int k = 0; k < n_states; k++) {
    G(k, k) -= lambda(k);
  }
  
  Type ill;
  
  if (istart + U > t.size() || istart + U > s.size()) {
    error("istart + U out of bounds");
  }
  
  // get data for this individual
  vector<Type> ti = t.segment(istart, U);
  vector<int> si = s.segment(istart, U);
  istart = istart + U;
  
  // Bounds checking
  for(int idx = 0; idx < si.size(); idx++) {
    if (si(idx) < 0 || si(idx) >= n_states) {
      error("State index out of bounds");
    }
  }
  
  if (U == 0) {
    matrix<Type> A = atomic::expm( matrix<Type>(G * (Time-t0) ) );

    ill = (pi.matrix().transpose() * A).sum();
  } else { 
    // Pre-allocate vectors
    vector<Type> ts(ti.size() + 2);
    ts(0) = t0;
    ts.segment(1, ti.size()) = ti;
    ts(ti.size() + 1) = Time;
    
    vector<Type> tau(U + 1);
    for(int k = 0; k < U + 1; k++) {
      tau(k) = ts(k + 1) - ts(k);
    }
    
    // Pre-allocate w vector
    vector<Type> w(U + 1);
    
    // First interval
    matrix<Type> A = atomic::expm( matrix<Type>(G * tau(0)) );
    vector<Type> omega1 = A.col(si(0));
    w(0) = (omega1.array() * pi.array()).sum() * lambda(si(0));
    
    // Middle intervals
    for(int u = 1; u < U; u++) {
      matrix<Type> A_u = atomic::expm( matrix<Type>(G * tau(u)) );
      w(u) = A_u(si(u-1), si(u)) * lambda(si(u));
    }
    
    // Final interval
    matrix<Type> A_final = atomic::expm( matrix<Type>(G * tau(U)) );
    Type omega_sum = Type(0);
    for (int k = 0; k < n_states; k++) {
      omega_sum += A_final(si(U-1), k);
    }
    w(U) = omega_sum;
    
    // Compute product more efficiently
    ill = Type(1);
    for(int k = 0; k < w.size(); k++) {
      ill *= w(k);
    }
  }
  
  return ill;
}

template<class Type>
Type objective_function<Type>::operator() () {
  
  PARAMETER_VECTOR(theta);
  
  // actual data
  DATA_IVECTOR(Us);
  DATA_VECTOR(t);
  DATA_IVECTOR(s);
  
  DATA_SCALAR(t0);
  DATA_SCALAR(Time);
  DATA_INTEGER(n_states);         // State space mesh size
  DATA_INTEGER(n_ac);             // Activity center mesh size  
  DATA_INTEGER(n_indiv);
  
  DATA_VECTOR(Q_fixed);
  DATA_VECTOR(lambda_fixed);
  DATA_VECTOR(camera_count);
  DATA_SCALAR(mesh_spacing);
  
  //Coordinates for state space and AC points
  DATA_MATRIX(state_coords);    // n_states x 2
  DATA_MATRIX(ac_coords);        // n_ac x 2
  
  
  Type ll = Type(0);
  vector<Type> etheta = exp(theta);
  
  
  // Parameters
  Type sigma_sq = etheta(0);           // Diffusion coefficient
  Type theta_reversion = etheta(1);    // Mean reversion strength
  
  // Pre-allocate matrices outside loops
  vector<Type> lambda(n_states);
  matrix<Type> Q(n_states, n_states);
  
  // OPTIMIZED: Build lambda once (based on state space mesh)
  for(int j = 0; j < n_states; j++) {
    lambda(j) = (lambda_fixed(j) == 1) ? Type(0) : etheta(2)*camera_count(j);
  }
  
  
  // OPTIMIZED: Pre-allocate storage for Q matrices and stationary distributions
  std::vector<matrix<Type>> Q_list;
  std::vector<vector<Type>> pi_list;
  Q_list.reserve(n_ac);  // Reserve space for AC mesh size, not state mesh size
  pi_list.reserve(n_ac);

  
  
  // Build Q matrices for each AC point
  for(int ac_idx = 0; ac_idx < n_ac; ac_idx++) {
    
    // Get activity center coordinates
    Type mu_x = ac_coords(ac_idx, 0);
    Type mu_y = ac_coords(ac_idx, 1);
    
    Type base_rate = sigma_sq / (Type(2) * mesh_spacing * mesh_spacing);
    
    int Qit = 0;
    for(int j = 0; j < n_states; j++) {
      
      // Current state coordinates
      Type x_j = state_coords(j, 0);
      Type y_j = state_coords(j, 1);
      
      for(int k = 0; k < n_states; k++) {
        
        if(Q_fixed(Qit) == 1) {
          Q(j, k) = Type(0);
        } 
        else if (j == k) {
          Q(j, k) = Type(0);  // Diagonal set later
        } 
        else {
          // Neighbor state coordinates
          Type x_k = state_coords(k, 0);
          Type y_k = state_coords(k, 1);
          
          // Jump direction vector
          Type jump_x = x_k - x_j;
          Type jump_y = y_k - y_j;
          
          // Displacement from activity center
          Type disp_x = x_j - mu_x;
          Type disp_y = y_j - mu_y;
          
          // Drift dot product: (r_j - mu) · (r_k - r_j)
          Type drift_dot = disp_x * jump_x + disp_y * jump_y;
          
          // OU process rate formula
          Type drift_component = theta_reversion * drift_dot / (Type(2) * mesh_spacing);
          Q(j, k) = base_rate - drift_component;
          
          // Handle potential negative rates
          if(Q(j, k) < Type(0)) {
            Q(j, k) = Type(1e-10);
          } 
        }
        Qit++;
      } 
    }
    
    // Set diagonal elements
    for(int j = 0; j < n_states; j++) {
      Type row_sum = Type(0);
      for(int k = 0; k < n_states; k++) {
        if(j != k) row_sum += Q(j, k);
      } 
      Q(j, j) = -row_sum;
    } 
    
    vector<Type> pi = StationaryDist(Q);
    
    Q_list.push_back(Q);
    pi_list.push_back(pi);
  }
  
  // OPTIMIZED: More efficient individual loop
  for(int i = 0; i < n_indiv - 1; i++) {
    
    // Calculate starting index more efficiently
    int original_istart = 0;
    for(int prev = 0; prev < i; prev++) {
      original_istart += Us(prev);
    }
    
    Type ill = Type(0);
    
    // Loop over AC integration points (not all states!)
    for(int j = 0; j < n_ac; j++) {
      int istart = original_istart;
      
      Type contrib = ind_likelihood<Type>(i, n_states, lambda, Q_list[j], pi_list[j], 
                                          Us(i), t, s, istart, t0, Time, theta, j);
      
      ill += contrib ;
    }
    
    ll -= log(ill / n_ac );
  }
  
  // Add unobserved-individual likelihood component
  Type p_obs = Type(0);
  for(int j = 0; j < n_ac; j++) {
    matrix<Type> G = Q_list[j];
    for(int i = 0; i < n_states; i++) {
      G(i, i) -= lambda(i);
    }
    matrix<Type> A = atomic::expm( matrix<Type>(G * Time) );
    
    vector<Type> ones(n_states);
    ones.setOnes();
    p_obs += (1-(pi_list[j].matrix().transpose() * A * ones.matrix())(0,0));
  }
  
  if (p_obs <= Type(1e-12)) p_obs = Type(1e-12);
  
  ll += (n_indiv - 1) * log(p_obs / n_ac);
  
  
  ADREPORT(etheta);
  
  return ll;
}
