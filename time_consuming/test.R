source("Gen_randomdata.R")
source("SCR.R")
source("Wu.R")

test = function(xxx, N ,n , alpha,Data,p_T,p_N, type,J_T,J_N){
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
  
  if(type == "MA"){
    xi = lapply(zeta, function(zeta){return( 0.8*zeta[1:n+1,]+0.6*zeta[1:n,]
    )})
  }
  
  lattice_time = seq(1/n,1,1/n)
  time_varying = lattice_time^2 
  time_varying = matrix(rep(time_varying,each = N),n,N,T)
  #result <- t(sapply(1:n, function(j) m_real + time_varying[j,]))
  
  R = lapply(xi, function(i){ t(sapply(1:n, function(j) m_real + time_varying[j,])) +#rep(1,n)%*%t(m_real) + 
      t(sapply(1:n, function(j) Reduce(`+`, map2(i[j,],phi, function(x,y) x*y))))})

  Y_1 = sapply(R, function(j){return(lapply(epsilon, function(i,j) i + j, j = j))})
  Y_2 = sapply(R, function(j){return(
    lapply(epsilon, function(i,j) Sigma2*i + j, j = j)
  )})
  Y = c(Y_1,Y_2)
  rm(epsilon)

  result = c(sapply(1:length(Y_1), function(i){StaionarySCB(Y_1[[i]],n,N,alpha,J_T,J_N,p_T,p_N)}))
  result_Wu = c(sapply(1:length(Y_1), function(i){Wu(Y_1[[i]],n,N,alpha)}))
  rm(Y)
  cat("N:",N,"n:",n,"times:",xxx," \r")
  return(c(result,result_Wu))
}
