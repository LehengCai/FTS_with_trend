source("Gen_randomdata.R")
source("Wu.R")

test = function(xxx, N ,n , alpha,Data, type){
  m_real = matrix(Data[[1]],n,N,byrow = T) 
  phi = Data[[2]]
 
  
  Sigma1 = 0.5
  Sigma2 = 0.5*matrix(Data[[3]],nrow = n,byrow = T)
  
  rm(Data)
  
  k = length(phi)
  Randomdata = Gen_randomdata(n,N,k) 
  epsilon = Randomdata[[1]]  
  zeta = Randomdata[[2]]
  rm(Randomdata)  
  
  if(type =="independent"){
    xi = zeta
  }
  
  if(type == "MA3"){
    xi = lapply(zeta, function(zeta){return( 0.8*zeta[1:n+3,]+0.4*zeta[1:n+2,]-
                                               0.4*zeta[1:n+1,]+0.2*zeta[1:n,]
    )})
  }
  
  R = lapply(xi, function(i){ m_real +
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