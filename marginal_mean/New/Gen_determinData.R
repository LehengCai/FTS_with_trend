Gen_determinData = function(N){
  lattice = seq(1/N,1,1/N)
  mean_function = function(v){
    return(1.5* sin(3*pi*(v+ 0.5))+2*v^3)
    #return (sin(2*pi*(v-0.5)))   # 均值函数设定备选
  }
  lambda = c(2,1/2)
  phi_function_1 = function(v){
    return ((sqrt(2)*sqrt(lambda[1]))*-cos(pi*v-pi/2))  #lambda_1 =2
  }
  phi_function_2 = function(v){
    return ((sqrt(2)*sqrt(lambda[2]))*sin(pi*v-pi/2)) #lambda_2 = 1/2
  }
  # phi_function_3 = function(s){
  #   return ((sqrt(2)*sqrt(lambda[3]))*sin(4*pi*s))  #lambda_3 = 0.5555
  # }
  # phi_function_4 = function(s){
  #   return ((sqrt(2)*sqrt(lambda[4]))*cos(4*pi*s)) #lambda_4 = 0.3125
  # }
  
  sigma_function = function(s){
    return((5-exp(-(s)))/(5+exp(-(s))))
  }
  phi_1 = phi_function_1(lattice)
  phi_2 = phi_function_2(lattice)
  # phi_3 = phi_function_3(lattice)
  # phi_4 = phi_function_4(lattice)

  m = mean_function(lattice) # 截距(方便展示)#产生均值
  Sigma = sigma_function(lattice)
  
  phi = list(phi_1,phi_2 )#,phi_3,phi_4)#,phi_5,phi_6,phi_7,phi_8,phi_9,phi_10,
             #phi_11,phi_12,phi_13,phi_14,phi_15,phi_16,phi_17,phi_18,phi_19,phi_20)
  return(list(m,phi,Sigma))
  
}
