setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls( ))
library(Rcpp)
library(RcppEigen)
library(stringr)
source("Gen_determinData.R",encoding = "utf-8")
source("test.R",encoding = "utf-8") 
source("estimate.R")
 

set.seed(123)
p_N = 4;p_T = 3 
N = 100 
MCMC = 500 
alpha = c(0.01,0.05,0.1,0.2) 
J_N =  floor(1.5*N^(1/4)*log(log(N))) #1.2-1.75
f = function(n){
  J_T =  max(0,floor(n^(1/6))-p_T) 
  type = "MA"
  Data = Gen_determinData(N)
  result = matrix(rowMeans(sapply(1:MCMC,test,N = N,n = n, 
                          alpha = alpha,Data = Data, p_T = p_T, p_N = p_N, 
                          type = type,  J_T=J_T,J_N=J_N),na.rm = T),ncol=length(alpha),byrow =T) 
  cat('J_T:',J_T,'J_N:',J_N,'---',result)
  write.csv(result,file=str_c('result/T=',n,'N=',N,'JT=',J_T,'JN=',J_N,'pT=',p_T,'pN=',p_N,'.csv'))
}
n = c(400,800,1200)
A = sapply(n,f)

