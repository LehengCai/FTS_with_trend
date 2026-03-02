rm(list = ls( ))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(Rcpp)
library(RcppEigen)
library(stringr)
source("Gen_determinData.R")
source("test.R") 

 
set.seed(123)
n = c(400,800,1200)  
MCMC = 500 
alpha = c(0.01,0.05,0.1,0.2)  
 
f = function(n){
  N = 100
  type = "MA" 
  
  Data = Gen_determinData(N)
  result = matrix(rowMeans(sapply(1:MCMC,test,N = N,n = n, 
                                  alpha = alpha,Data = Data,
                                  type = type),na.rm = T),ncol=length(alpha),byrow =T) 
  cat('n:',n,' N:',N,'---',result)
  write.csv(result,file=str_c('result/Wu_T=',n,'N=',N,'.csv'))
}

A = sapply(n,f)