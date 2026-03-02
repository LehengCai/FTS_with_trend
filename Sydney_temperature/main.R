setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls( ))

library(scatterplot3d)
library(plotly)
library(splines)
library(rgl)
library(MASS)

df = read.csv('temperature_data.csv') 
Y = as.matrix(df[,2:366])
n <- nrow(Y)      
N <- ncol(Y)      

 

par(lwd = 0.1)
idx <- seq(1960, 2008, by = 10)
k   <- length(idx)
cols <- colorRampPalette(c("blue4", "deepskyblue"))(k)


s3d <- scatterplot3d(
  x = rep(idx[1], N),
  y = 1:N,
  z = Y[idx[1]-1959+1,],
  type = "l",
  color = cols[1],
  box = FALSE,
  grid = FALSE,
  xlim = c(1959,2008),
  ylim = c(1,365),
  zlim = c(10,40),
  ylab = 'Day',xlab = 'Year',zlab = 'Temperature' 
)

for(j in 2:k){
  i <- idx[j]
  s3d$points3d(
    x = rep(i, N),
    y = 1:N,
    z = Y[i-1959+1, ],
    type = "l",
    col = cols[j]
  )
}



set.seed(123)
N = ncol(Y)
n = nrow(Y);n
J_N = 0
p_N = 4
J_T = 0
p_T = 4
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

N0 = 50
lattice.eval = seq(1/N0,1,1/N0)
B_N.eval = bs(lattice.eval,knots=knotsN[-c(1,(J_N+2))],degree=p_N-1,intercept=TRUE,Boundary.knots=knotsN[c(1,(J_N+2))])
B_T.eval = bs(lattice.eval,knots=knotsT[-c(1,(J_T+2))],degree=p_T-1,intercept=TRUE,Boundary.knots=knotsT[c(1,(J_T+2))])
B_N.eval = B_N
B_T.eval = B_T

B.eval = kronecker(B_T.eval,B_N.eval)
m_hat.eval =  matrix(B.eval%*%a2,n,N,T)
m_hat =  matrix(a3,n,N,T) 
persp3d(m_hat,col='blue')

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
  ans = ans * BB/L 
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
alpha = 0.05
Q = quantile(apply(abs(Our_result),2,max),1-alpha);Q

persp3d(m_hat.eval,zlim=c(10,40), col = "red", shade = 0.5)
persp3d(m_hat.eval+Q*matrix(diag_sd,n,N,T)/sqrt(n ) ,zlim=c(10,40),shade = 0.5)
persp3d(m_hat.eval-Q*matrix(diag_sd,n,N,T)/sqrt(n ) ,zlim=c(10,40), shade = 0.5)

write.csv(m_hat.eval,'z1.csv')
write.csv(m_hat.eval+Q*matrix(diag_sd,n,N,T)/sqrt(n ),'z2.csv')
write.csv(m_hat.eval-Q*matrix(diag_sd,n,N,T)/sqrt(n ),'z3.csv')

alpha = 0.01
Q = quantile(apply(abs(Our_result),2,max),1-alpha);Q

persp3d(m_hat.eval,zlim=c(10,40), col = "red", shade = 0.5)
persp3d(m_hat.eval+Q*matrix(diag_sd,n,N,T)/sqrt(n ) ,zlim=c(10,40),shade = 0.5)
persp3d(m_hat.eval-Q*matrix(diag_sd,n,N,T)/sqrt(n ) ,zlim=c(10,40), shade = 0.5)

write.csv(m_hat.eval,'z1.csv')
write.csv(m_hat.eval+Q*matrix(diag_sd,n,N,T)/sqrt(n ),'z22.csv')
write.csv(m_hat.eval-Q*matrix(diag_sd,n,N,T)/sqrt(n ),'z33.csv')



################# 
m_tilde = B_N%*%solve(t(B_N)%*%B_N)%*%t(B_N)%*%colMeans(Y)
m0 = matrix(rep(m_tilde,50),nrow=50,byrow = T)
write.csv(m0,'z4.csv')
persp3d(m0 ,zlim=c(10,40), col = 'red')

source("Wu.R")
result_kernel = Wu(Y,n,N,alpha=0.05)
persp3d(result_kernel$m_hat,zlim=c(10,40), col = "red", shade = 0.5)
persp3d(result_kernel$upper ,zlim=c(10,40),shade = 0.5)
persp3d(result_kernel$lower ,zlim=c(10,40), shade = 0.5)

write.csv(result_kernel$m_hat,'wu1.csv')
write.csv(result_kernel$upper,'wu2.csv')
write.csv(result_kernel$lower,'wu3.csv')
 
