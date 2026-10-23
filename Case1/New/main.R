setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls( ))
library(Rcpp)
library(RcppEigen)
library(stringr)
source("Gen_determinData.R",encoding = "utf-8")
source("test.R",encoding = "utf-8") 
source("estimate.R")


p_N = 4 
p_T = 3 
set.seed(123)
n = c(400,800,1200)
MCMC = 500
alpha = c(0.01,0.05,0.1,0.2) 
N = 100
J_N =  floor(1.5*N^(1/4)*log(log(N)))
f = function(n){
  J_T =  max(0,floor(n^(1/6))-p_T) 
  type = "MA"
  Data = Gen_determinData(N)
  # result = matrix(rowMeans(sapply(1:MCMC,test,N = N,n = n, 
  #                         alpha = alpha,Data = Data, p_T = p_T, p_N = p_N, 
  #                         type = type,  J_T=J_T,J_N=J_N),na.rm = T),
  #                 ncol=length(alpha),byrow =T)
  
  raw_result = sapply(1:MCMC,test,N = N,n = n,
    alpha = alpha,Data = Data,p_T = p_T,p_N = p_N,
    type = type,J_T = J_T,J_N = J_N
  )
  avg_result = rowMeans(raw_result, na.rm = TRUE)
  result = matrix(
    avg_result,
    nrow = 6,
    ncol = 4 * length(alpha),
    byrow = TRUE
  )
  
  rownames(result) = c(
    "Homo_Normal",
    "Homo_Uniform",
    "Homo_Laplace",
    "Hetero_Normal",
    "Hetero_Uniform",
    "Hetero_Laplace"
  )
  colnames(result) = c(
    
    paste0("Coverage_full_alpha_", alpha),
    
    paste0("Volume_full_alpha_", alpha),
    
    paste0("Coverage_interior_alpha_", alpha),
    
    paste0("Volume_interior_alpha_", alpha)
    
  )
  
  write.csv(result,file=str_c('result/T=',n,'N=',N,'.csv'))
}
A = sapply(n,f)

