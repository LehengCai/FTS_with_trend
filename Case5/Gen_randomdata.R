library(tidyverse)
library(GeneralizedHyperbolic)
Gen_randomdata = function(n , N ,k){ 
  para1 = list(
    list( mean = 0 , sd = 0.5)
  )
  f=c( "rnorm" )
  epsilon = invoke_map(f,para1,n*N) 
  g = function(X,n,N){
    return(matrix(X,n,N))
  }
  
  para2 = list(
    list(mean = 0, sd = 1), 
    list(min = -sqrt(3), max = sqrt(3)),
    list(mu = 0 , alpha = sqrt(2)/2, beta = sqrt(2)/2) 
    
  )
  f2 = c("rnorm","runif","rskewlap")
  
  M = 10
  xi = invoke_map(f2,para2,(n+M)*k)
  
  
  return(list(lapply(epsilon,g, n = n, N = N),
              lapply(xi,g, n = n+M ,N = k))) 
}