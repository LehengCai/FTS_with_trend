
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
  raw_result = sapply(
    1:MCMC,
    test,
    N = N,
    n = n,
    alpha = alpha,
    Data = Data,
    type = type,
    x0 = x0
  )
  
  avg_result = rowMeans(
    raw_result,
    na.rm = TRUE
  )
  
  result = matrix(
    avg_result,
    nrow = 6,
    ncol = 3 * length(alpha),
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
    paste0("Coverage_interior_alpha_", alpha),
    paste0("Average_width_alpha_", alpha)
  )
  cat('n:',n,'N:',N,'---',result)
  write.csv(result,file=str_c('result/T=',n,'N=',N,'x0=',x0,'.csv'))
}
A = sapply(n,f)

