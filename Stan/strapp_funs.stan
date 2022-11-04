
//' PWE log likelihood
//' @param y               failure / censoring time
//' @param eta             linear predictor
//' @param lambda          vector of baseline hazards
//' @param log_lambda      vector of log(lambda)
//' @param breaks          (J+1)-dim vector giving intervals
//' @param j               index of the interval for which the individual failed / was censored: j \in {1, ..., J}
//' @param death_ind       integer giving whether individual died (death_ind == 1) or was censored (death_ind == 0)
//' @param lambda_d_breaks (J-1)-dim vector giving lambda[j] * (s[j] - s[j-1]), j = 1, ..., J-1
real pwe_lpdf(real y, real eta, vector lambda, vector log_lambda, vector breaks, int j, int death_ind, vector cumblhaz ) {
  real loglik;
  
  // Initialize cumhaz to lambda[j] ( y - s[j] )
  real cumhaz = lambda[j] * ( y - breaks[j] );
  
  // add on (sum_{g=1}^{j-1} lambda[j] * ( s[j] - s[j-1] )
  if ( j > 1 )
    cumhaz += cumblhaz[j-1];
  
  // Initialize loglik = - cumhaz * exp(x'beta) = log(Survival)
  loglik = -cumhaz * exp(eta);
  
  // If individual experienced event, add log hazard: log(lambda[j]) + x'beta
  if ( death_ind == 1 )
    loglik += log_lambda[j] + eta;
    
  return(loglik);
}


//' Promotion time cure rate likelihood
//' L(lambda, beta | y) = [theta[i] * f(y[i] | psi[i])]^(nu[i]) * exp(-theta[i] * F(y[i] | psi[i]))
//' @param y               failure / censoring time
//' @param eta             linear predictor
//' @param lambda          vector of baseline hazards
//' @param log_lambda      vector of log(lambda)
//' @param breaks          (J+1)-dim vector giving intervals
//' @param j               index of the interval for which the individual failed / was censored: j \in {1, ..., J}
//' @param death_ind       integer giving whether individual died (death_ind == 1) or was censored (death_ind == 0)
//' @param cumblhaz       (J-1)-dim vector giving cumulative sum of lambda[j-1] * (s[j] - s[j-1]), j = 2, ..., J
real promotion_time_cure_rate_lpdf(real y, real eta, vector lambda, vector log_lambda, vector breaks, int j, int death_ind, vector cumblhaz ) {
  real cdf;
  real loglik;
  real cumhaz;
  // compute cumulative hazard
  cumhaz = lambda[j] * (y - breaks[j]);
  if ( j > 1)
    cumhaz += cumblhaz[j-1];
  
  // Compute 1 - survival time = cdf for individual
  cdf = -expm1(-cumhaz);       // -expm1(x) = -(exp(x) - 1) = 1 - exp(x)
  
  // Initialize log likelihood = -exp(eta) * (1 - S(y))
  loglik = -exp(eta) * cdf;
  
  // If individual has event, add on linear predictor and log of density:
  //   log(f(y)) = log(h(y)*S(y)) = log(h(y)) - H(y) since S(y) = exp(-H(y))
  if ( death_ind )
    loglik += eta + log_lambda[j] - cumhaz;
    
  return(loglik);
}



// function to return square root of matrix via spectral decomposition
matrix matrix_sqrt(matrix A) {
  int p                 = rows(A);
  vector[p]   A_values  = eigenvalues_sym(A);
  matrix[p,p] A_vectors = eigenvectors_sym(A);
  return A_vectors * diag_pre_multiply(sqrt(A_values), A_vectors');
}


// function to return square root fisher information matrix dropping intercepts
// dropping intercept
matrix sqrt_fisher_logistic(matrix X, vector beta, real intercept) {
  int n = rows(X);
  vector[n] W = inv_logit(intercept + X*beta);
  W = W .* (1 - W); // p * (1 - p)
  return matrix_sqrt( X' * diag_pre_multiply(W, X) );
}


matrix sqrt_fisher_promotiontime(matrix X, vector beta, real intercept) {
  int n = rows(X);
  vector[n] W = exp(intercept + X * beta);
  return matrix_sqrt( X' * diag_pre_multiply(W, X) );
}

vector system_strapp(vector beta, vector beta0_int0_int, real[] x_r, int[] x_i) {
  int n = x_i[1];
  int p = x_i[2];
  matrix[n,p] X0 = to_matrix(x_r, n, p);
  vector[p] beta_hist = beta0_int0_int[1:p];
  real intercept_hist = beta0_int0_int[p+1];
  real intercept = beta0_int0_int[p+2];
  return 
      sqrt_fisher_promotiontime(X0, beta, intercept) * beta
    - sqrt_fisher_logistic(X0, beta_hist, intercept_hist) * beta_hist
  ;
}


vector system_genstrapp(vector beta, vector beta0_int0_c0_int, real[] x_r, int[] x_i) {
  int n = x_i[1];
  int p = x_i[2];
  matrix[n,p] X0 = to_matrix(x_r, n, p);
  vector[p] beta_hist = beta0_int0_c0_int[1:p];
  real intercept_hist = beta0_int0_c0_int[p+1];
  vector[p] c0 = beta0_int0_c0_int[(p+2):(2*p + 1)];
  real intercept = beta0_int0_c0_int[2*p + 2];
  return (
      sqrt_fisher_promotiontime(X0, beta, intercept) * beta
    - sqrt_fisher_logistic(X0, beta_hist, intercept_hist) * beta_hist
    - c0
  );
}




matrix A_kron_I(matrix A, int p) {
  int m = rows(A);
  matrix[m*p,m*p] C = rep_matrix(0, m*p, m*p);
  for ( i in 1:m ) {
    for ( j in 1:m ) {
      int row_start = (i - 1) * p + 1;
      int row_end = (i - 1) * p + p;
      int col_start = (j - 1) * p + 1;
      int col_end = (j - 1) * p + p;
      // diagonal(C[row_start:row_end, col_start:col_end]) = rep_vector(A[i,j], p);
      C[row_start:row_end, col_start:col_end] = add_diag(C[row_start:row_end, col_start:col_end], A[i,j]);
    }
  }
  return C;
}
matrix I_kron_A(matrix A, int p) {
  int m = rows(A);
  matrix[p*m, p*m] C = rep_matrix(0, m*p, m*p);
  for ( k in 1:p ) {
    int start = (k-1)*m + 1;
    int end   = k*m;
    C[start:end,start:end] = A;
  }
  return C;
}

real log_abs_det_jacob_promotion(vector beta1, real intercept1, real sigma0, matrix X, real log_det_invsqrt_XtX) {
  int n                   = rows(X);
  int p                   = rows(beta1);
  vector[n]   exp_eta1    = exp( intercept1 + X*beta1 );
  matrix[p,p] fisher      = X' * diag_pre_multiply(exp(intercept1 + X*beta1), X);
  matrix[p,p] sqrt_fisher = matrix_sqrt(fisher);
  vector[n] d_W_j_vec;
  matrix[p,p] d_fisher_j;
  vector[p*p] vec_d_sqrt_fisher_j;
  matrix[p,p] d_sqrt_fisher_j;
  matrix[p*p,p*p] inv_kron_sum = inverse_spd( I_kron_A(sqrt_fisher, p) + A_kron_I(sqrt_fisher, p) );
  matrix[p,p] jacob;
  // initialize result of jacobian to LHS of transformation
  // compute derivative of square root of fisher info
  for ( j in 1:p ) {
    d_W_j_vec  = exp_eta1 .* X[, j];                   // derivative of W matrix wrt beta1[j]
    d_fisher_j = X' * diag_pre_multiply(d_W_j_vec, X); // derivative of fisher info wrt beta1[j]
      
      // Compute vectorized derivative of sqrt of fisher --> convert to pxp matrix
      vec_d_sqrt_fisher_j = inv_kron_sum * to_vector(d_fisher_j);
      d_sqrt_fisher_j     = to_matrix(vec_d_sqrt_fisher_j, p, p);
      
      // Compute jacobian
      jacob[, j] = d_sqrt_fisher_j * beta1 + sqrt_fisher[, j];
  }
  return p * log(sigma0) + log_det_invsqrt_XtX + log_determinant(jacob);
}
matrix jacob_promotion(vector beta1, real intercept1, real sigma0, matrix X, matrix invsqrt_XtX) {
  int n                   = rows(X);
  int p                   = rows(beta1);
  vector[n]   exp_eta1    = exp( intercept1 + X*beta1 );
  matrix[p,p] fisher      = X' * diag_pre_multiply(exp(intercept1 + X*beta1), X);
  matrix[p,p] sqrt_fisher = matrix_sqrt(fisher);
  vector[n] d_W_j_vec;
  matrix[p,p] d_fisher_j;
  vector[p*p] vec_d_sqrt_fisher_j;
  matrix[p,p] d_sqrt_fisher_j;
  matrix[p*p,p*p] inv_kron_sum = inverse_spd( I_kron_A(sqrt_fisher, p) + A_kron_I(sqrt_fisher, p) );
  matrix[p,p] jacob;
  // initialize result of jacobian to LHS of transformation
  // compute derivative of square root of fisher info
  for ( j in 1:p ) {
    d_W_j_vec  = exp_eta1 .* X[, j];                   // derivative of W matrix wrt beta1[j]
    d_fisher_j = X' * diag_pre_multiply(d_W_j_vec, X); // derivative of fisher info wrt beta1[j]
      
      // Compute vectorized derivative of sqrt of fisher --> convert to pxp matrix
      vec_d_sqrt_fisher_j = inv_kron_sum * to_vector(d_fisher_j);
      d_sqrt_fisher_j     = to_matrix(vec_d_sqrt_fisher_j, p, p);
      
      // Compute jacobian
      jacob[, j] = d_sqrt_fisher_j * beta1 + sqrt_fisher[, j];
  }
  return sigma0 * invsqrt_XtX * jacob;
}



