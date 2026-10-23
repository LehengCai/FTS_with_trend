library(scatterplot3d)
library(rainbow)
df <- Australiafertility
y <- df$y
library(plotly)
data <- data.frame(Date=rep(seq(as.Date("1921-01-01"),by="year",length.out=95),each=35),
                   Value=as.numeric(y),
                   age=rep(as.numeric(row.names(y)),95))
library(scatterplot3d)
par(lwd = 0.2)
s3d <- scatterplot3d(rep(1921,35), as.numeric(row.names(y)), y[,1], highlight.3d=TRUE, type="l", main=NULL, box=F, xlim=c(1921,2015),
                     ylab = 'Age',xlab = 'Year',zlab = 'Fertility Rates',
                     grid=F)
for(i in 1922:2015){
  s3d$points3d(rep(i,35), as.numeric(row.names(y)), y[,i-1920], highlight.3d=TRUE, type="l", main=NULL, box=F)
}


# df <- Australiafertility
# y <- df$y
# library(plotly)
# data <- data.frame(Date=rep(seq(as.Date("1921-01-01"),by="year",length.out=95),each=35),
#                    Value=as.numeric(y),
#                    age=rep(as.numeric(row.names(y)),95))
# library(scatterplot3d)
# par(lwd = 0.2)
# s3d <- scatterplot3d(rep(1921,35), as.numeric(row.names(y)), y[,1], highlight.3d=F, pch=20 , main=NULL, box=F, xlim=c(1921,2015),
#                      ylab = 'Age',xlab = 'Year',zlab = 'Fertility Rates',
#                      grid=F )
# for(i in 1922:2015){
#   s3d$points3d(rep(i,35), as.numeric(row.names(y)), y[,i-1920], highlight.3d=TRUE, pch=20, main=NULL, box=F)
# }
###################################################################################
###SCRs
 
library(stringr)
library(MASS)
library(splines)
library(rgl)
Y = df$y

set.seed(123)
Y = t(Y)
N = ncol(Y)
n = nrow(Y);n
J_N = 4
p_N = 4
J_T = 4
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


persp3d(m_hat.eval,zlim=c(-10,300), col = "red", shade = 0.5)
persp3d(m_hat.eval+Q*matrix(diag_sd,n,N,T)/sqrt(n ) ,zlim=c(-10,300),shade = 0.5)
persp3d(m_hat.eval-Q*matrix(diag_sd,n,N,T)/sqrt(n ) ,zlim=c(-10,300), shade = 0.5)
# write.csv(m_hat.eval,'z1.csv')
# write.csv(m_hat.eval+Q*matrix(diag_sd,n,N,T)/sqrt(n ),'z2.csv')
# write.csv(m_hat.eval-Q*matrix(diag_sd,n,N,T)/sqrt(n ),'z3.csv')
###################################################################################
#SCB for some fixed x=25,30,35
Y = df$y
set.seed(123)
Y = t(Y)
N = ncol(Y)
n = nrow(Y);n
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

##########fixed x=25,30,35
b1 = t(B_N)%*% ans %*%B_N/N^2
V_N = t(B_N)%*%B_N/N
b2 = solve(V_N)%*% b1 %*%  solve(V_N)
set.seed(123)
x0 = 35 ###25 30
x0 = (x0-15)/35
B_x0 = bs(x0,knots=knotsN[-c(1,(J_N+2))],degree=p_N-1,intercept=TRUE,Boundary.knots=knotsN[c(1,(J_N+2))])
V_T = t(B_T)%*%B_T/n
Sigma = solve(V_T)*as.numeric(B_x0%*%b2%*%t(B_x0))

diag_sd =  sqrt(diag(B_T.eval%*%Sigma%*%t(B_T.eval))) 
Our_bootstrap = mvrnorm(1000,rep(0,ncol(Sigma)),Sigma)
Our_result = B_T.eval%*%t(Our_bootstrap)/diag_sd
alpha = 0.05
Q = quantile(apply(abs(Our_result),2,max),1-alpha);Q
m_hat_x0 = as.vector(kronecker(B_T.eval,B_x0)%*%a2)

ggplot(NULL)+
  geom_line(aes(1921:2015,m_hat_x0),color='blue',size=0.8)+
  geom_line(aes(1921:2015,m_hat_x0+Q*diag_sd /sqrt(n)),color='red',linetype=2,size=0.8)+
  geom_line(aes(1921:2015,m_hat_x0-Q*diag_sd /sqrt(n)),color='red',linetype=2,size=0.8)+
  ylim(c(30,270))+
  xlab('Year')+ylab('Fertility Rates')+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+ theme(legend.position = 'none')

##########fixed tau=1930,1960,1990
b1 = t(B_N)%*% ans %*%B_N/N^2
V_N = t(B_N)%*%B_N/N
b2 = solve(V_N)%*% b1 %*%  solve(V_N)
set.seed(123)
tau0 = 1960 ###1930 1990
tau0 = (tau0-1921)/95
B_tau0 = bs(tau0,knots=knotsT[-c(1,(J_T+2))],degree=p_T-1,intercept=TRUE,Boundary.knots=knotsT[c(1,(J_T+2))])
V_T = t(B_T)%*%B_T/n
Sigma = as.numeric(B_tau0%*%solve(V_T)%*%t(B_tau0)) * b2 

diag_sd =  sqrt(diag(B_N.eval%*%Sigma%*%t(B_N.eval))) 
Our_bootstrap = mvrnorm(1000,rep(0,ncol(Sigma)),Sigma)
Our_result = B_N.eval%*%t(Our_bootstrap)/diag_sd
alpha = 0.05
Q = quantile(apply(abs(Our_result),2,max),1-alpha);Q
m_hat_tau0 = as.vector(kronecker(B_tau0,B_N.eval)%*%a2)

ggplot(NULL)+
  geom_line(aes(15:49,m_hat_tau0),color='blue',size=0.8)+
  geom_line(aes(15:49,m_hat_tau0+Q*diag_sd /sqrt(n)),color='red',linetype=2,size=0.8)+
  geom_line(aes(15:49,m_hat_tau0-Q*diag_sd /sqrt(n)),color='red',linetype=2,size=0.8)+
  ylim(c(-10,270))+
  xlab('Age')+ylab('Fertility Rates')+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+ theme(legend.position = 'none')

####################################################################################
#additivity testing
Y = df$y
set.seed(123)
Y = t(Y)
N = ncol(Y)
n = nrow(Y);n
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
B.eval = kronecker(B_T.eval,B_N.eval)
m_hat.eval =  matrix(B.eval%*%a2,N0,N0,T)
m_hat_dt.eval = matrix(rep(colMeans(m_hat.eval),N0),N0,N0,T)
m_hat_dx.eval = matrix(rep(rowMeans(m_hat.eval),N0),N0,N0,F)
m_hat_dtdx.eval = matrix(mean(m_hat.eval),N0,N0)
U = m_hat.eval - m_hat_dt.eval - m_hat_dx.eval + m_hat_dtdx.eval


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

b1 = t(B_N)%*% ans %*%B_N/N^2
V_N = t(B_N)%*%B_N/N
b2 = solve(V_N)%*% b1 %*%  solve(V_N)
V_T = t(B_T)%*%B_T/n
Sigma = kronecker(solve(V_T),b2)

#####################testing
B_T.int <- sapply(1:(p_T+J_T), function(i) {
  integrate(function(x){
    return( bs(x,knots=knotsT[-c(1,(J_T+2))],degree=p_T-1,intercept=TRUE,Boundary.knots=knotsT[c(1,(J_T+2))])[,i]
    )
  }, lower = 0, upper = 1)$value
})
B_T.int = matrix(rep(B_T.int,N0),nrow = N0,byrow = T)

B_N.int <- sapply(1:(p_N+J_N), function(i) {
  integrate(function(x){
    return( bs(x,knots=knotsN[-c(1,(J_N+2))],degree=p_N-1,intercept=TRUE,Boundary.knots=knotsN[c(1,(J_N+2))])[,i]
    )
  }, lower = 0, upper = 1)$value
})
B_N.int = matrix(rep(B_N.int,N0),nrow = N0,byrow = T)
c1 = kronecker(B_T.eval,B_N.eval)
c2 = kronecker(B_T.eval,B_N.int)
c3 = kronecker(B_T.int,B_N.eval)
c4 = kronecker(B_T.int,B_N.int)

diag_sd =  sqrt(diag(
  (c1-c2-c3+c4)%*%Sigma%*%t((c1-c2-c3+c4))
))
Our_bootstrap = mvrnorm(1000,rep(0,ncol(Sigma)),Sigma)
Our_result = (c1-c2-c3+c4)%*%t(Our_bootstrap)/diag_sd
Q = quantile(apply(abs(Our_result),2,max),1-alpha);Q
flag = rep(0,length(alpha))
for(i in 1:length(alpha)){
  flag[i] = as.numeric(any(U+Q[i]*matrix(diag_sd,N0,N0,T)/sqrt(n) < 0 ) ||
                         any(U-Q[i]*matrix(diag_sd,N0,N0,T)/sqrt(n)> 0))
}

flag
max(abs(U/ matrix(diag_sd,N0,N0,T) *sqrt(n)))
