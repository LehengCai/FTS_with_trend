#This is code of Generate data 
library(tidyverse)
library(GeneralizedHyperbolic)
Gen_randomdata = function(n , N ,k){ #样本个数，总观测格点个数,特征函数个数
  para1 = list(
    # list( mean = 0 , sd = 1), # 正态参数1
    # list( min = -1*sqrt(3), max = 1*sqrt(3)),#均匀分布参数1
    # list( mu = 0 ,alpha = 1/sqrt(2),beta = 1/sqrt(2)), #双指数分布参数1
    list( mean = 0 , sd = 0.5)# 正态参数2
    #list( min = -0.1*sqrt(3), max = 0.1*sqrt(3)),#均匀分布参数2
    #list( mu = 0 ,alpha = 0.1/sqrt(2),beta = 0.1/sqrt(2)) #双指数分布参数2
    
  )
  f=c( "rnorm" )#,"runif","rskewlap")#,"rnorm","runif","rskewlap")
  epsilon = invoke_map(f,para1,n*N)  #独立同分布的epsilon  
  g = function(X,n,N){
    return(matrix(X,n,N))
  }
  
  para2 = list(
    list(mean = 0, sd = 1),  #正态分布参数
    list(min = -sqrt(3), max = sqrt(3)), #均匀分布参数
    list(mu = 0 , alpha = sqrt(2)/2, beta = sqrt(2)/2) #双指数分布参数
    
  )
  #f2 = c("rnorm")
  f2 = c("rnorm","runif","rskewlap")
  
  M = 10
  xi = invoke_map(f2,para2,(n+M)*k) #产生独立同分布的xi (n*k个)##
  
  
  return(list(lapply(epsilon,g, n = n, N = N),
              lapply(xi,g, n = n+M ,N = k))) #t(N * n) t(k*n)矩阵 记得要转置 
  
}