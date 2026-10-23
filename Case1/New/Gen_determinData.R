Gen_determinData = function(N){
  lattice = seq(1/N,1,1/N)
  mean_function = function(v){
    return(1.5* sin(3*pi*(v+ 0.5))+2*v^3)
  }
  lambda = c(2,1/2)
  phi_function_1 = function(v){
    return ((sqrt(2)*sqrt(lambda[1]))*-cos(pi*v-pi/2))  
  }
  phi_function_2 = function(v){
    return ((sqrt(2)*sqrt(lambda[2]))*sin(pi*v-pi/2)) 
  }

  
  sigma_function = function(s){
    return((5-exp(-(s)))/(5+exp(-(s))))
  }
  phi_1 = phi_function_1(lattice)
  phi_2 = phi_function_2(lattice)


  m = mean_function(lattice) 
  Sigma = sigma_function(lattice)
  
  phi = list(phi_1,phi_2 )
  return(list(m,phi,Sigma))
  
}
