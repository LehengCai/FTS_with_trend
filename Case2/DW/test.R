source("Gen_randomdata.R")
source("Wu.R")

test = function(xxx, N ,n , alpha,Data,type){
  m_real = Data[[1]]
  phi = Data[[2]]
  
  Sigma1 = 1
  Sigma2 = 1.2*matrix(rep(Data[[3]],n),nrow = n,byrow = T)
  
  rm(Data) 
  
  k = length(phi)
  Randomdata = Gen_randomdata(n,N,k) 
  epsilon = Randomdata[[1]]  
  zeta = Randomdata[[2]] 
  rm(Randomdata)
  
  if(type =="independent"){
    xi = zeta
  }
  
  if(type == 'AR'){
    xi = lapply(zeta,function(zeta){
      X <- matrix(0,n,k)
      X[1,] <- rnorm(k)
      for (t in 2:n) {
        X[t,] <- -1/sqrt(2) * X[t - 1,] + zeta[t,]/sqrt(2)
      }
      return(X)
    })
  }
  
  lattice_time = seq(1/n,1,1/n)
  time_varying = 3*lattice_time*exp(-lattice_time) 
  time_varying = matrix(rep(time_varying,each = N),n,N,T)
  
  R = lapply(xi, function(i){ t(sapply(1:n, function(j) m_real + time_varying[j,])) + 
      t(sapply(1:n, function(j) Reduce(`+`, map2(i[j,],phi, function(x,y) x*y))))})
  
 
  
  Y_1 = sapply(R, function(j){return(lapply(epsilon, function(i,j) i + j, j = j))})
  Y_2 = sapply(R, function(j){return(
    lapply(epsilon, function(i,j) Sigma2*i + j, j = j)
  )})
  
  Y = c(Y_1,Y_2) 
  rm(epsilon)
  
  result = c(sapply(1:length(Y), function(i){Wu(Y[[i]],n,N,alpha)}))
  rm(Y)
  cat("N:",N,"n:",n,"times:",xxx," \r")
  return(result)
}
