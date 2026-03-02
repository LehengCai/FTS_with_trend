
rm(list = ls( ))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(Rcpp)
library(RcppEigen)
library(stringr)
source("Gen_determinData.R",encoding = "utf-8")
source("test.R",encoding = "utf-8")
source("estimate.R")


x0 = 1/2
 
set.seed(123)
n = c(400,800,1200)
MCMC = 500 
alpha = c(0.01,0.05,0.1,0.2)  
 
f = function(n){
  N = 100
  type = "MA"  
  Data = Gen_determinData(N)
  result = matrix(rowMeans(sapply(1:MCMC,test,N = N,n = n, 
                          alpha = alpha,Data = Data, type = type,  x0=x0),na.rm = T),
                  ncol=2*length(alpha),byrow =T) 
  cat('n:',n,'N:',N,'---',result)
  write.csv(result,file=str_c('result/T=',n,'N=',N,'x0=',x0,'.csv'))
}
A = sapply(n,f)

