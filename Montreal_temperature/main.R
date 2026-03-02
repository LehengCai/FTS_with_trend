library(fda)
data(MontrealTemp)

library(scatterplot3d)
par(lwd = 0.2)
y = MontrealTemp
s3d <- scatterplot3d(rep(1961,365), 1:365, y[1,], highlight.3d=TRUE, type="l", main=NULL, box=F, xlim=c(1961,1994),zlim=c(-32,32),ylim=c(1,365),
                     ylab = 'Day',xlab = 'Year',zlab = 'Temperature',
                     grid=F)
for(i in 1962:1994){
  s3d$points3d(rep(i,365), 1:365, y[i-1960,], highlight.3d=TRUE, type="l", main=NULL, box=F)
}


library(MASS)
library(splines)
library(rgl)
library(stringr)
Y = MontrealTemp
set.seed(123)
N = ncol(Y);N
n = nrow(Y);n
J_N = 4
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


persp3d(m_hat.eval,zlim=c(-30,30), col = "red", shade = 0.5)
persp3d(m_hat.eval+Q*matrix(diag_sd,n,N,T)/sqrt(n ) ,zlim=c(-30,30),shade = 0.5)
persp3d(m_hat.eval-Q*matrix(diag_sd,n,N,T)/sqrt(n ) ,zlim=c(-30,30), shade = 0.5)

# write.csv(m_hat.eval,'z1.csv')
# write.csv(m_hat.eval+Q*matrix(diag_sd,n,N,T)/sqrt(n ),'z2.csv')
# write.csv(m_hat.eval-Q*matrix(diag_sd,n,N,T)/sqrt(n ),'z3.csv')

#####################
m_tilde = B_N%*%solve(t(B_N)%*%B_N)%*%t(B_N)%*%colMeans(Y)
m0 = matrix(rep(m_tilde,34),nrow=34,byrow = T)
#write.csv(m0,'z4.csv')


set.seed(123)
Our_bootstrap = mvrnorm(10000,rep(0,ncol(Sigma)),Sigma)
Our_result = B.eval%*%t(Our_bootstrap)/diag_sd
alpha = 0.166 #0.05
Q = quantile(apply(abs(Our_result),2,max),1-alpha);Q
min(m_hat.eval+Q*matrix(diag_sd,n,N,T)/sqrt(n )-m0)
min(m0-(m_hat.eval-Q*matrix(diag_sd,n,N,T)/sqrt(n )))
