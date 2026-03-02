source("Gen_randomdata.R")
source("JOE_SCB.R")

test = function(xxx, N ,n , alpha,Data,p_T, type,J_T,x0){
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
  
  if(type == "MA"){
    xi = lapply(zeta, function(zeta){return( 0.8*zeta[1:n+1,]+0.6*zeta[1:n,]
  )})
  }

  lattice_time = seq(1/n,1,1/n)
  time_varying = lattice_time^2 
  time_varying = matrix(rep(time_varying,each = N),n,N,T)

  
  R = lapply(xi, function(i){ t(sapply(1:n, function(j) m_real + time_varying[j,])) +
      t(sapply(1:n, function(j) Reduce(`+`, map2(i[j,],phi, function(x,y) x*y))))})

  Y_1 = sapply(R, function(j){return(lapply(epsilon, function(i,j) i + j, j = j))})
  Y_2 = sapply(R, function(j){return(
    lapply(epsilon, function(i,j) Sigma2*i + j, j = j)
  )})
  
  Y = c(Y_1,Y_2) 
  rm(epsilon)
  # YY = Y
  # Y = YY[[1]]  
  
  result = c(sapply(1:length(Y), function(i){fixed_x0_SCB(Y[[i]],n,N,J_T,p_T,alpha,x0)}))
  rm(Y)
  cat("N:",N,"n:",n,"times:",xxx," \r")
  return(result)
}
