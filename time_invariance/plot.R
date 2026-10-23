setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls( ))
library(Rcpp)
library(RcppEigen)
library(stringr)
source("Gen_determinData.R" )
source("Gen_randomData.R" )

N = 400;n = 50
Data = Gen_determinData(N)
m_real = Data[[1]]
phi = Data[[2]]
Sigma1 = 1
Sigma2 = 1.2*matrix(rep(Data[[3]],n),nrow = n,byrow = T)
k = length(phi)
Randomdata = Gen_randomdata(n,N,k) 
epsilon = Randomdata[[1]]  
zeta = Randomdata[[2]] 
rm(Randomdata)
xi = lapply(zeta, function(zeta){return( 0.8*zeta[1:n+1,]+0.6*zeta[1:n,])})
lattice_time = seq(1/n,1,1/n)
time_varying = lattice_time^2 * 0
time_varying = matrix(rep(time_varying,each = N),n,N,T)
R = lapply(xi, function(i){ t(sapply(1:n, function(j) m_real + time_varying[j,])) + 
    t(sapply(1:n, function(j) Reduce(`+`, map2(i[j,],phi, function(x,y) x*y))))})
Y_1 = sapply(R, function(j){return(lapply(epsilon, function(i,j) i + j, j = j))})
Y_2 = sapply(R, function(j){return(
  lapply(epsilon, function(i,j) Sigma2*i + j, j = j)
)})
Y = c(Y_1,Y_2) 
y = Y[[1]]

library(scatterplot3d)
par(lwd = 0.2)
s3d <- scatterplot3d(rep(1,400), 1:400, y[1,], highlight.3d=TRUE, type="l", main=NULL, box=F, xlim=c(1,50),zlim=c(0,7),ylim=c(1,400),
                     grid=F,xlab='',ylab='',zlab='')
for(i in 1:50){
  s3d$points3d(rep(i,400), 1:400, y[i,], highlight.3d=TRUE, type="l", main=NULL, box=F)
}

