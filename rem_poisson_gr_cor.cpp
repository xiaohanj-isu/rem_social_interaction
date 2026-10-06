#include <TMB.hpp>

//============================================================
// TMB REM Poisson model
//
// Model:
//
//   y_m ~ Poisson(mu_m)
//
//   log(mu_m) =
//       X_m beta
//     + g[giver_m]
//     + r[receiver_m]
//     + gpe[giver_m]
//     + rpe[receiver_m]
//     + sg[group_m]
//     + offset_m
//
// Random effects:
//
//   [g, r] ~ N(0, G ⊗ Sigma_g)
//
//   Sigma_g =
//      [ var_g    cov_gr ]
//      [ cov_gr   var_r  ]
//
//   [gpe_i, rpe_i] ~ N(0, Sigma_pe)
//
//   Sigma_pe =
//      [ var_gpe       cov_gpe_rpe ]
//      [ cov_gpe_rpe   var_rpe     ]
//
//   sg_k ~ N(0, var_sg)
//============================================================

template<class Type>
Type objective_function<Type>::operator() ()
{
  //============================================================
  // 1. Data
  //============================================================
  
  DATA_VECTOR(y);
  
  DATA_MATRIX(X);
  
  DATA_VECTOR(offset);
  
  DATA_IVECTOR(giver_id);
  DATA_IVECTOR(receiver_id);
  
  DATA_IVECTOR(group_id);
  
  DATA_SPARSE_MATRIX(Ginv);
  
  DATA_SCALAR(logdet_G);
  
  DATA_INTEGER(n_animal);
  DATA_INTEGER(n_group);
  
  //============================================================
  // 2. Parameters
  //============================================================
  
  PARAMETER_VECTOR(beta);
  
  // genetic effects
  PARAMETER_VECTOR(g);   // giver genetic effect
  PARAMETER_VECTOR(r);   // receiver genetic effect
  
  // permanent environmental effects
  PARAMETER_VECTOR(gpe);
  PARAMETER_VECTOR(rpe);
  
  // social group effect
  PARAMETER_VECTOR(sg);
  
  // variance components on log scale
  PARAMETER(log_var_g);
  PARAMETER(log_var_r);
  
  PARAMETER(log_var_gpe);
  PARAMETER(log_var_rpe);
  
  PARAMETER(log_var_sg);
  
  // correlation parameters on unconstrained scale
  PARAMETER(theta_gr);
  PARAMETER(theta_gpe_rpe);
  
  //============================================================
  // 3. Transform variance components
  //============================================================
  
  Type var_g   = exp(log_var_g);
  Type var_r   = exp(log_var_r);
  
  Type var_gpe = exp(log_var_gpe);
  Type var_rpe = exp(log_var_rpe);
  
  Type var_sg  = exp(log_var_sg);
  
  // transform theta to correlation in (-1, 1)
  Type rho_gr =
    Type(2.0) / (Type(1.0) + exp(-theta_gr)) - Type(1.0);
  
  Type rho_gpe_rpe =
    Type(2.0) / (Type(1.0) + exp(-theta_gpe_rpe)) - Type(1.0);
  
  // covariance
  Type cov_gr = rho_gr * sqrt(var_g * var_r);
  Type cov_gpe_rpe = rho_gpe_rpe * sqrt(var_gpe * var_rpe);
  
  // standard deviation for social group
  Type sd_sg = sqrt(var_sg);
  
  // constants
  Type pi = Type(3.14159265358979323846);
  Type log2pi = log(Type(2.0) * pi);
  
  //============================================================
  // 4. Negative log likelihood
  //============================================================
  
  Type nll = Type(0.0);
  
  //============================================================
  // 5. Correlated genetic prior for giver and receiver effects
  //
  // Let U = [g, r], an n_animal x 2 matrix.
  //
  // Assumption:
  //
  //   vec(U) ~ N(0, G ⊗ Sigma_g)
  //
  // where:
  //
  //   Sigma_g =
  //      [ var_g    cov_gr ]
  //      [ cov_gr   var_r  ]
  //
  // The negative log density is:
  //
  //   0.5 * n_animal * log|Sigma_g|
  // + 0.5 * 2 * log|G|
  // + 0.5 * 2*n_animal * log(2*pi)
  // + 0.5 * trace(Sigma_g^{-1} U' G^{-1} U)
  //
  //============================================================
  
  vector<Type> Ginv_g = Ginv * g;
  vector<Type> Ginv_r = Ginv * r;
  
  Type g_Ginv_g = Type(0.0);
  Type r_Ginv_r = Type(0.0);
  Type g_Ginv_r = Type(0.0);
  
  for(int i = 0; i < n_animal; i++){
    g_Ginv_g += g(i) * Ginv_g(i);
    r_Ginv_r += r(i) * Ginv_r(i);
    g_Ginv_r += g(i) * Ginv_r(i);
  }
  
  Type det_Sigma_g = var_g * var_r - cov_gr * cov_gr;
  
  Type invSigma_g_11 = var_r / det_Sigma_g;
  Type invSigma_g_22 = var_g / det_Sigma_g;
  Type invSigma_g_12 = -cov_gr / det_Sigma_g;
  
  nll += Type(0.5) * (
    Type(n_animal) * log(det_Sigma_g)
    + Type(2.0) * logdet_G
  + Type(2.0) * Type(n_animal) * log2pi
  );
  
  nll += Type(0.5) * (
    invSigma_g_11 * g_Ginv_g
  + Type(2.0) * invSigma_g_12 * g_Ginv_r
  + invSigma_g_22 * r_Ginv_r
  );
  
  //============================================================
  // 6. Correlated PE prior for giver PE and receiver PE
  //
  // For each animal i:
  //
  //   [gpe_i, rpe_i]' ~ N(0, Sigma_pe)
  //
  // where:
  //
  //   Sigma_pe =
  //      [ var_gpe       cov_gpe_rpe ]
  //      [ cov_gpe_rpe   var_rpe     ]
  //
  //============================================================
  
  Type det_Sigma_pe = var_gpe * var_rpe - cov_gpe_rpe * cov_gpe_rpe;
  
  Type invSigma_pe_11 = var_rpe / det_Sigma_pe;
  Type invSigma_pe_22 = var_gpe / det_Sigma_pe;
  Type invSigma_pe_12 = -cov_gpe_rpe / det_Sigma_pe;
  
  Type gpe2 = Type(0.0);
  Type rpe2 = Type(0.0);
  Type gpe_rpe = Type(0.0);
  
  for(int i = 0; i < n_animal; i++){
    gpe2 += gpe(i) * gpe(i);
    rpe2 += rpe(i) * rpe(i);
    gpe_rpe += gpe(i) * rpe(i);
  }
  
  nll += Type(0.5) * (
    Type(n_animal) * log(det_Sigma_pe)
    + Type(2.0) * Type(n_animal) * log2pi
  );
  
  nll += Type(0.5) * (
    invSigma_pe_11 * gpe2
  + Type(2.0) * invSigma_pe_12 * gpe_rpe
  + invSigma_pe_22 * rpe2
  );
  
  //============================================================
  // 7. Social group prior
  //
  //   sg_k ~ N(0, var_sg)
  //
  //============================================================
  
  for(int k = 0; k < n_group; k++){
    nll -= dnorm(sg(k), Type(0.0), sd_sg, true);
  }
  
  //============================================================
  // 8. Poisson likelihood
  //
  // For risk-set observation m:
  //
  //   y_m ~ Poisson(mu_m)
  //
  //   log(mu_m) =
  //       X_m beta
  //     + g[giver_m]
  //     + r[receiver_m]
  //     + gpe[giver_m]
  //     + rpe[receiver_m]
  //     + sg[group_m]
  //     + offset_m
  //
  //============================================================
  
  int n = y.size();
  
  for(int m = 0; m < n; m++){
    
    Type xb = Type(0.0);
    
    for(int p = 0; p < X.cols(); p++){
      xb += X(m, p) * beta(p);
    }
    
    Type eta =
      xb
      + g(giver_id(m))
      + r(receiver_id(m))
      + gpe(giver_id(m))
      + rpe(receiver_id(m))
      + sg(group_id(m))
      + offset(m);
      
      Type mu = exp(eta);
      
      nll -= dpois(y(m), mu, true);
  }
  
  //============================================================
  // 9. Report transformed parameters
  //============================================================
  
  ADREPORT(var_g);
  ADREPORT(var_r);
  ADREPORT(cov_gr);
  ADREPORT(rho_gr);
  
  ADREPORT(var_gpe);
  ADREPORT(var_rpe);
  ADREPORT(cov_gpe_rpe);
  ADREPORT(rho_gpe_rpe);
  
  ADREPORT(var_sg);
  
  //============================================================
  // 10. Report random effects
  //============================================================
  
  REPORT(g);
  REPORT(r);
  
  REPORT(gpe);
  REPORT(rpe);
  
  REPORT(sg);
  
  //============================================================
  // 11. Return total negative log likelihood
  //============================================================
  
  return nll;
}