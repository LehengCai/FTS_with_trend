# This is the code of FPCA analysis.
library(MASS)
library(rgl)
library(zoo)

fixed_x0_SCB = function(Y,n,N,alpha,J_T,J_N,p_T,p_N,x0){ 
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
  m_hat =  matrix(a3,n,N,T) 
  m_hat_x0 = m_hat[,round(x0*N)]

  ####################################################
  trajectory_coef = solve(t(B_N)%*%B_N)%*%t(B_N)%*%t(Y)###uncentered smoothed 
  pre_smoothing = t(sapply(1:n, function(t,B_T,B_N,trajectory_coef,a2){
    B_t =  matrix(B_T[t,],1)
    return(B_N%*%trajectory_coef[,t]-kronecker(B_t,B_N) %*% a2)
  },B_T,B_N,trajectory_coef,a2))###demeaned smoothed 
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

  BB_values = min(2,floor(n^(1/5))): ceiling(2*n^(1/5))
  ans_values = lapply(BB_values, calculate_ans)
  
  # best_BB_index <- which.min( rollapply(sapply(ans_values, function(x) diag(x)[round(N*x0)] ),
  #                 width = 3, FUN = var, fill = NA, align = "center") )
  best_BB_index = which.max(sapply(ans_values, function(x) mean(diag(x)) ))
  best_BB = BB_values[best_BB_index]
  ans = ans_values[[best_BB_index]]

  ##############################################################
  b1 = t(B_N)%*% ans %*%B_N/N^2
  V_N = t(B_N)%*%B_N/N
  b2 = solve(V_N)%*% b1 %*%  solve(V_N)
  V_T = t(B_T)%*%B_T/n
  Gamma_hat = B_N %*% b2 %*% t(B_N)
  Sigma =  solve(V_T)*Gamma_hat[round(x0*N),round(x0*N)]
  diag_sd =  sqrt(diag(B_T%*%Sigma%*%t(B_T))) 
  Our_bootstrap = mvrnorm(1000,rep(0,ncol(Sigma)),Sigma)
  Our_result = B_T%*%t(Our_bootstrap)/diag_sd
  Q = quantile(apply(abs(Our_result),2,max),1-alpha);Q
  
  time_varying = lattice_time^2 
  mean_function = function(v){
    return(1.5* sin(3*pi*(v+ 0.5))+2*v^3)
  }
  m_real = mean_function(x0)+time_varying
  
  flag = rep(0,length(alpha))
  width = rep(0,length(alpha))
  for(i in 1:length(alpha)){
    flag[i] = as.numeric(all(m_hat_x0 + Q[i]* diag_sd /sqrt(n)>=m_real)&&
                        all(m_hat_x0 - Q[i]* diag_sd /sqrt(n)<=m_real)) 
    width[i] = mean(Q[i]* diag_sd /sqrt(n))*2
  }
  
  return(c(flag,width))
  
  plot(lattice_time,m_hat_x0,'l',ylim=c(-1.5,3))
  lines(lattice_time,lattice_time^2+mean_function(x0),col='red')
  lines(lattice_time,  m_hat_x0+Q[2]* diag_sd /sqrt(n),col='blue')
  lines(lattice_time,  m_hat_x0-Q[2]* diag_sd /sqrt(n),col='blue')
 
}


 