Gen_determinData = function(N){
  lattice = seq(1/N,1,1/N)
  mean_function = function(v){
    return(1.5* sin(3*pi*(v+ 0.5))+2*v^3)
     
  }
  
  lambda <- 2 * (1:100)^(-2)
  
  make_phi <- function(k) {
    freq <- ceiling(k / 2) * 2 * pi
    trig <- if (k %% 2 == 1) sin else cos
    
    function(s) {
      sqrt(2 * lambda[k]) * trig(freq * s)
    }
  }

  phi_functions <- lapply(1:100, make_phi)
  phi <- lapply(phi_functions, function(f) f(lattice))
  

  sigma_function = function(s){
    return((5-exp(-(s)))/(5+exp(-(s))))
  }


  m = mean_function(lattice)  
  Sigma = sigma_function(lattice)
  
  return(list(m,phi,Sigma))
  
}
