

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls( ))
library(Rcpp)
library(RcppEigen)
library(stringr)
source("Gen_determinData.R",encoding = "utf-8")
source("test.R",encoding = "utf-8")
source("estimate.R")


set.seed(123)
p_N = 3;p_T = 3
N = 400;n = 50
MCMC = 500
alpha = c(0.01,0.05,0.1,0.2)
J_N =  floor(1*N^(1/5))
f = function(n){
  J_T = max(0, floor(n^(1/6)) - p_T)
  type = "MA"

  Data = Gen_determinData(N)

  raw_result = sapply(
    1:MCMC,
    test,
    N = N,
    n = n,
    alpha = alpha,
    Data = Data,
    p_T = p_T,
    p_N = p_N,
    type = type,
    J_T = J_T,
    J_N = J_N
  )

  avg_result = rowMeans(
    raw_result,
    na.rm = TRUE
  )

  ##########################################################
  # 6 scenarios
  # each scenario:
  # 4 original coverage + 4 adjusted coverage
  ##########################################################

  result = matrix(
    avg_result,
    nrow = 6,
    ncol = 2 * length(alpha),
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
    paste0("Coverage_original_alpha_", alpha),
    paste0("Coverage_adjusted_alpha_", alpha)
  )

  print(result)

  write.csv(
    result,
    file = str_c(
      "result/T=", n,
      "N=", N,
      "JN=",J_N,
      ".csv"
    )
  )

  return(result)
}
A = sapply(n,f)