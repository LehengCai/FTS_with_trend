library(MASS)
StaionarySCB = function(Y, n, N, alpha,
                        J_T, J_N, p_T, p_N,
                        N0 = 50,
                        B_boot = 1000){
  
  ############################################################
  # 1. B-spline bases
  ############################################################
  
  lattice = seq(1/N, 1, 1/N)
  
  knotsN = seq(0, 1, length = J_N + 2)
  
  B_N = bs(
    lattice,
    knots = knotsN[-c(1, J_N + 2)],
    degree = p_N - 1,
    intercept = TRUE,
    Boundary.knots = c(0, 1)
  )
  
  
  lattice_time = seq(1/n, 1, 1/n)
  
  knotsT = seq(0, 1, length = J_T + 2)
  
  # IMPORTANT:
  # unrestricted temporal spline estimator
  # Do NOT use B_T = 1 for the time-invariance test
  B_T = bs(
    lattice_time,
    knots = knotsT[-c(1, J_T + 2)],
    degree = p_T - 1,
    intercept = TRUE,
    Boundary.knots = c(0, 1)
  )
  
  
  ############################################################
  # 2. Unrestricted estimator m_hat(tau,x)
  ############################################################
  
  B = kronecker(B_T, B_N)
  
  coef_hat =
    solve(t(B) %*% B) %*%
    t(B) %*%
    c(t(Y))
  
  m_hat =
    matrix(
      B %*% coef_hat,
      nrow = n,
      ncol = N,
      byrow = TRUE
    )
  
  
  ############################################################
  # 3. Estimate m_0(x)
  #
  # Appendix F:
  # m_hat_0(x) = T^{-1} sum_t eta_hat_t(x)
  ############################################################
  
  trajectory_coef =
    solve(t(B_N) %*% B_N) %*%
    t(B_N) %*%
    t(Y)
  
  m0_vec =
    as.vector(
      B_N %*%
        solve(t(B_N) %*% B_N) %*%
        t(B_N) %*%
        colMeans(Y)
    )
  
  m0_surface =
    matrix(
      rep(m0_vec, n),
      nrow = n,
      ncol = N,
      byrow = TRUE
    )
  
  
  ############################################################
  # Optional numerical check:
  # Appendix F implies
  # m_hat_0(x) = T^{-1} sum_t m_hat(t/T,x)
  ############################################################
  
  check_m0 =
    max(abs(m0_vec - colMeans(m_hat)))
  
  if(check_m0 > 1e-6){
    warning(
      paste(
        "Two representations of m0 differ by",
        check_m0
      )
    )
  }
  
  
  ############################################################
  # 4. Estimate long-run covariance Gamma
  #
  # Residuals must be based on the unrestricted estimator
  ############################################################
  
  pre_smoothing =
    t(
      sapply(
        1:n,
        function(tt){
          
          B_t = matrix(B_T[tt, ], nrow = 1)
          
          as.vector(
            B_N %*% trajectory_coef[, tt] -
              kronecker(B_t, B_N) %*% coef_hat
          )
        }
      )
    )
  
  
  calculate_ans <- function(BB){
    
    ans = matrix(0, N, N)
    
    for(l in 1:(n - BB + 1)){
      
      start_index = l
      end_index = l + BB - 1
      
      block_mean =
        colMeans(
          pre_smoothing[start_index:end_index, ,
                        drop = FALSE]
        )
      
      ans =
        ans +
        block_mean %*% t(block_mean)
    }
    
    ans =
      ans * BB / (n - BB + 1)
    
    return(ans)
  }
  
  
  BB_values =
    min(2, floor(n^(1/5))):
    ceiling(2 * n^(1/5) * log(n))
  
  ans_values =
    lapply(
      BB_values,
      calculate_ans
    )
  
  best_BB_index =
    which.max(
      sapply(
        ans_values,
        function(x) mean(diag(x))
      )
    )
  
  best_BB =
    BB_values[best_BB_index]
  
  ans =
    ans_values[[best_BB_index]]
  
  
  ############################################################
  # 5. Spline representation of Gamma_hat
  ############################################################
  
  b1 =
    t(B_N) %*%
    ans %*%
    B_N / N^2
  
  V_N =
    t(B_N) %*%
    B_N / N
  
  b2 =
    solve(V_N) %*%
    b1 %*%
    solve(V_N)
  
  V_T =
    t(B_T) %*%
    B_T / n
  
  Sigma =
    kronecker(
      solve(V_T),
      b2
    )
  
  
  ############################################################
  # 6. Centered temporal basis
  #
  # B_T^c(tau)
  # = B_T(tau) - bar B_T
  ############################################################
  
  B_T_bar =
    colMeans(B_T)
  
  B_T_center =
    sweep(
      B_T,
      MARGIN = 2,
      STATS = B_T_bar,
      FUN = "-"
    )
  
  B_center =
    kronecker(
      B_T_center,
      B_N
    )
  
  
  ############################################################
  # 7. Estimated variance of centered Gaussian process
  ############################################################
  
  B_center_Sigma =
    B_center %*% Sigma
  
  diag_var_center =
    rowSums(
      B_center_Sigma *
        B_center
    )
  
  # numerical protection
  diag_var_center =
    pmax(
      diag_var_center,
      1e-12
    )
  
  diag_sd_center =
    sqrt(diag_var_center)
  
  
  ############################################################
  # 8. Corrected observed test statistic
  #
  # sqrt(T/J_T) / sqrt(Omega_tilde)
  # simplifies numerically to
  #
  # sqrt(T) * |m_hat-m0| / diag_sd_center
  ############################################################
  
  difference_surface =
    m_hat -
    m0_surface
  
  difference_vec =
    c(t(difference_surface))
  
  test_stat =
    max(
      abs(difference_vec) *
        sqrt(n) /
        diag_sd_center
    )
  
  
  ############################################################
  # 9. Gaussian simulation for corrected critical value
  ############################################################
  
  Gaussian_coef =
    mvrnorm(
      B_boot,
      mu = rep(0, ncol(Sigma)),
      Sigma = Sigma
    )
  
  Z_center =
    B_center %*%
    t(Gaussian_coef)
  
  Z_center_std =
    sweep(
      Z_center,
      MARGIN = 1,
      STATS = diag_sd_center,
      FUN = "/"
    )
  
  bootstrap_stat =
    apply(
      abs(Z_center_std),
      2,
      max
    )
  
  
  ############################################################
  # 10. Critical values
  ############################################################
  
  Q_TI =
    as.numeric(
      quantile(
        bootstrap_stat,
        probs = 1 - alpha,
        names = FALSE
      )
    )
  
  
  ############################################################
  # 11. Original rejection decision
  ############################################################
  
  reject_TI_original =
    as.numeric(
      test_stat > Q_TI
    )
  
  
  ############################################################
  # 12. Small-sample adjustment
  #
  # Adjust the rejection boundary by
  #
  #       (1 + 20/T)
  #
  # i.e.
  #
  # Q_TI_adjusted = (1 + 20/T) * Q_TI
  ############################################################
  
  adjustment_factor =
    1 + 20 / n
  
  Q_TI_adjusted =
    adjustment_factor *
    Q_TI
  
  
  ############################################################
  # Adjusted rejection decision
  ############################################################
  
  reject_TI_adjusted =
    as.numeric(
      test_stat > Q_TI_adjusted
    )
  
  
  ############################################################
  # Optional Monte Carlo p-value
  ############################################################
  
  p_value_TI =
    (1 +
       sum(
         bootstrap_stat >= test_stat
       )) /
    (B_boot + 1)
  
  
  ############################################################
  # Return:
  #
  # alpha = c(0.01, 0.05)
  #
  # 1. Original rejection at alpha = 0.01
  # 2. Original rejection at alpha = 0.05
  # 3. Adjusted rejection at alpha = 0.01
  # 4. Adjusted rejection at alpha = 0.05
  ############################################################
  
  result =
    c(
      reject_TI_original,
      reject_TI_adjusted
    )
  
  names(result) =
    c(
      paste0(
        "TI_original_alpha_",
        alpha
      ),
      paste0(
        "TI_adjusted_alpha_",
        alpha
      )
    )
  
  return(result)
}

