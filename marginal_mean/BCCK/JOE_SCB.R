library(splines)
library(MASS)

fixed_x0_SCB = function(Y,n,N,J_T,p_T, alpha, x0){

  lattice_time = seq(1/n,1,1/n)
  
  knotsT = seq(0,1,length=J_T+2)  
  B_T = bs(lattice_time,knots=knotsT[-c(1,(J_T+2))],degree=p_T-1,intercept=TRUE,Boundary.knots=knotsT[c(1,(J_T+2))])
  V_inv = solve(t(B_T)%*%B_T/n)
  
  x=Y[,round(N*x0)]
  m_hat_x0 = as.vector( B_T %*% solve(t(B_T)%*%B_T)%*%t(B_T)%*% x  )
  
  e_hat = x - m_hat_x0
  Sigma = mean(e_hat^2)*V_inv
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