Gen_determinData = function(n, N){
  lattice = seq(1/N,1,1/N)
  mean_function = function(tau,x){
    return((1+tau)*(6*(x-1/2)^2+1))
  }
  
  lambda = 0.5^(c(0,1,2,3))
  phi_function_1 = function(s){
    return ((sqrt(lambda[1]))*sqrt(2)*cos(2*pi*s)) 
  }
  phi_function_2 = function(s){
    return ((sqrt(lambda[2]))*sqrt(2)*sin(2*pi*s))
  }
  phi_function_3 = function(s){
    return ((sqrt(lambda[3]))*sqrt(2)*cos(4*pi*s))  
  }
  phi_function_4 = function(s){
    return ((sqrt(lambda[4]))*sqrt(2)*sin(4*pi*s))  
  }
  
  sigma_function = function(tau,x){
    return(1.2*(0.2*tau+0.8)*(5-exp(-(x)))/(5+exp(-(x))))
  }
  phi_1 = phi_function_1(lattice)
  phi_2 = phi_function_2(lattice)
  phi_3 = phi_function_3(lattice)
  phi_4 = phi_function_4(lattice)
  
  
  
  m = mean_function(rep(seq(1/n,1,1/n),each = N), rep(lattice,n))
  Sigma = sigma_function(rep(seq(1/n,1,1/n),each = N), rep(lattice,n))
  
  phi = list(phi_1,phi_2,phi_3,phi_4) 
  return(list(m,phi,Sigma))
  
}