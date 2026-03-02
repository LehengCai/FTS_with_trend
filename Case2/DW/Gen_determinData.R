Gen_determinData = function(N){
  lattice = seq(1/N,1,1/N)
  mean_function = function(v){
     return(5-4*(v-0.5)^2)
  }
  lambda = c(1,2,3,4)^(-2)
 
  phi_function_1 = function(s){
    return ((sqrt(lambda[1]))*1) 
  }
  phi_function_2 = function(s){
    return ((sqrt(lambda[2]))*(2*s-1)*sqrt(3)) 
  }
  phi_function_3 = function(s){
    return ((sqrt(lambda[3]))*(6*s^2-6*s+1)*sqrt(5)) 
  }
  phi_function_4 = function(s){
    return ((sqrt(lambda[4]))*(20*s^3-30*s^2+12*s-1)*sqrt(7))
  }
  
  sigma_function = function(s){
    return((5-exp(-(s)))/(5+exp(-(s))))
  }
  phi_1 = phi_function_1(lattice)
  phi_2 = phi_function_2(lattice)
  phi_3 = phi_function_3(lattice)
  phi_4 = phi_function_4(lattice)
  
  m = mean_function(lattice) 
  Sigma = sigma_function(lattice)
  
  phi = list(phi_1,phi_2,phi_3,phi_4)
  return(list(m,phi,Sigma))
  
}