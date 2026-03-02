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
n = c(200,500,800,1200,1500,2000)
MCMC = 5
alpha = c(0.01,0.05,0.1,0.2) 
N = 100
J_N =  floor(1.5*N^(1/4)*log(log(N)))
f = function(n){
  J_T =  max(0,floor(n^(1/6))-p_T) 
  type = "MA"
  Data = Gen_determinData(N)
  result = rowMeans(sapply(1:MCMC,test,N = N,n = n, 
                          alpha = alpha,Data = Data, p_T = p_T, p_N = p_N, 
                          type = type,  J_T=J_T,J_N=J_N),na.rm = T) 
  return(result)
}
A = sapply(n,f);A


m = ncol(A)
library(ggplot2)
df <- data.frame(
  x = n,
  log_time = log(c(A[1, ], A[2, ])),
  Method = factor(rep(c("B-spline", "Kernel"), each = m))
)

p <- ggplot(df, aes(x = x, y = log_time, linetype = Method)) +
  geom_line(size = 0.8) +    
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(
    x = "Sample size (T)",
    y = "Running time (seconds, ln scale)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.line = element_line(size = 0.8),   
    axis.ticks = element_line(size = 0.8)    
  )+ylim(c(0,7))

print(p)
