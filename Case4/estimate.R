library(splines)

estimate = function (Ns,p,N){ 
  lattice = seq(1/N,1,1/N)
  knots2 = seq(0,1,length=Ns+2) 
  B_p = bs(lattice,knots=knots2[-c(1,(Ns+2))],degree=p-1,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns+2))])
  V = solve(t(B_p)%*%B_p) 
  m_pre  = MatMult(V,B_p)
  return(list(m_pre,B_p))
}