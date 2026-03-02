library(MASS)

reconstruct_matrix <- function(input_matrix) {
  input_matrix = (input_matrix+t(input_matrix))/2
  eigen_decomp <- eigen(input_matrix)
  eigenvalues <- Re(eigen_decomp$values) 
  eigenvectors <- eigen_decomp$vectors
  positive_eigenvalues <- eigenvalues[eigenvalues >= 0]
  positive_eigenvectors <- eigenvectors[, eigenvalues >= 0]
  reconstructed_matrix <- positive_eigenvectors %*% diag(positive_eigenvalues) %*% t(positive_eigenvectors)
  return(reconstructed_matrix)
}

StaionarySCB = function(Y,n,N,alpha,J_T,J_N,p_T,p_N,N0 = 50){
  lattice = seq(1/N,1,1/N)
  knotsN = seq(0,1,length=J_N+2)  
  B_N = bs(lattice,knots=knotsN[-c(1,(J_N+2))],degree=p_N-1,intercept=TRUE,Boundary.knots=knotsN[c(1,(J_N+2))])
  lattice_time = seq(1/n,1,1/n)
  knotsT = seq(0,1,length=J_T+2)  
  B_T = bs(lattice_time,knots=knotsT[-c(1,(J_T+2))],degree=p_T-1,intercept=TRUE,Boundary.knots=knotsT[c(1,(J_T+2))])
  B = kronecker(B_T,B_N)
  V = solve(t(B)%*%B)
  a1 = t(B)%*% c(t(Y))
  a2 = V %*% a1
  a3 = B %*% a2
  
  
  lattice.eval = seq(1/N0,1,1/N0)
  B_N.eval = bs(lattice.eval,knots=knotsN[-c(1,(J_N+2))],degree=p_N-1,intercept=TRUE,Boundary.knots=knotsN[c(1,(J_N+2))])
  B_T.eval = bs(lattice.eval,knots=knotsT[-c(1,(J_T+2))],degree=p_T-1,intercept=TRUE,Boundary.knots=knotsT[c(1,(J_T+2))])
  B.eval = kronecker(B_T.eval,B_N.eval)
  m_hat.eval =  matrix(B.eval%*%a2,N0,N0,T)
  m_hat =  matrix(a3,n,N,T) 
  
  ####################################################
  trajectory_coef = solve(t(B_N)%*%B_N)%*%t(B_N)%*%t(Y)
  pre_smoothing = t(sapply(1:n, function(t,B_T,B_N,trajectory_coef,a2){
    B_t =  matrix(B_T[t,],1)
    return(B_N%*%trajectory_coef[,t]-kronecker(B_t,B_N) %*% a2)
  },B_T,B_N,trajectory_coef,a2))
  calculate_ans <- function(BB) {
    L = floor(n/BB)
    ans = matrix(0, N, N)
    for (l in 1:L) {
      start_index = 1 + (l-1) * BB
      end_index = l * BB
      ans = ans + colMeans(pre_smoothing[start_index:end_index,]) %*%
        t(colMeans(pre_smoothing[start_index:end_index,]))
    }
    ans = ans * BB/L #- BB * colMeans(pre_smoothing[1:(L*BB),]) %*% t(colMeans(pre_smoothing[1:(L*BB),]))
    return(ans)
  }
  BB_values = min(2,floor(n^(1/5))): ceiling(n^(1/5)*log(n))
  ans_values = lapply(BB_values, calculate_ans)
  best_BB_index = which.max(sapply(ans_values, function(x) mean(diag(x))))
  best_BB = BB_values[best_BB_index]
  ans = ans_values[[best_BB_index]]
  ##############################################################

  b1 = t(B_N)%*% ans %*%B_N/N^2
  V_N = t(B_N)%*%B_N/N
  b2 = solve(V_N)%*% b1 %*%  solve(V_N)
  V_T = t(B_T)%*%B_T/n
  Sigma = kronecker(solve(V_T),b2)
  
  diag_sd =  sqrt(diag(B.eval%*%Sigma%*%t(B.eval))) 
  Our_bootstrap = mvrnorm(1000,rep(0,ncol(Sigma)),Sigma)
  Our_result = B.eval%*%t(Our_bootstrap)/diag_sd
  Q = quantile(apply(abs(Our_result),2,max),1-alpha);Q
  
  lattice.eval = seq(1/N0,1,1/N0)
  time_varying = 3*lattice.eval*exp(-lattice.eval)
  time_varying = matrix(rep(time_varying,each = N0),N0,N0,T)
  m_real = 5-4*(lattice.eval-0.5)^2
  m.eval <- t(sapply(1:N0, function(j) m_real + time_varying[j,]))
  
  flag = rep(0,length(alpha))
  for(i in 1:length(alpha)){
    flag[i] = as.numeric(all(m_hat.eval+Q[i]*matrix(diag_sd,N0,N0,T)/sqrt(n)>=m.eval)&&
                        all(m_hat.eval-Q[i]*matrix(diag_sd,N0,N0,T)/sqrt(n)<=m.eval)) 
  }
  
  return(flag)
}

