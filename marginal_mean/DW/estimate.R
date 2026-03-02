#This is the code of estimation
library(splines)
#library(pracma) # 算张量积更快
estimate = function (Ns,p,N){ 
  lattice = seq(1/N,1,1/N)
  knots2 = seq(0,1,length=Ns+2)  #节点位置
  B_p = bs(lattice,knots=knots2[-c(1,(Ns+2))],degree=p-1,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns+2))])
  # N*(N_s+p) 系数矩阵
  
  V = solve(t(B_p)%*%B_p) #经验内积
  m_pre  = MatMult(V,B_p) ##V%*%t(B_p)
  return(list(m_pre,B_p))
}