# rm(list=ls(all=TRUE))
library(matrixStats)
library(foreach)
library(doParallel)
# library(rgl)

wapply1 <- function(x, width, by = 1, FUN = mean, ...)
{
  FUN <- match.fun(FUN)
  if (is.null(by)) by <- width
  
  lenX <- length(x)
  SEQ1 <- seq(1, lenX - width + 1, by = by)
  SEQ2 <- lapply(SEQ1, function(x) x:(x + width - 1))
  
  OUT <- lapply(SEQ2, function(a) FUN(x[a], ...))
  OUT <- base:::simplify2array(OUT, higher = TRUE)
  return(OUT)
}


######This version: Single Kernel for boundary, NW estimator for interior, and new MGCV
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
EKernel=function(t)
{
  return(0.75*(1-t^2)*(abs(t)<=1))
}

llinear=function(y,t,b)
{ n=length(y)
zz=((1:n)/n-t)/b
M1=sum(EKernel(zz)*zz)/n/b
M0=sum(EKernel(zz))/n/b
M2=sum(EKernel(zz)*zz^2)/n/b
return(sum((M2*EKernel(zz)-M1*EKernel(zz)*zz)/(M0*M2-M1^2)*y)/b/n)
}

HDtest=function(y_matrix,m=10,lit=2000,rep_time)##This version saves operation time but consumes space
{
  pp=dim(y_matrix)[1]
  TT=dim(y_matrix)[2]
  width=2*floor(m/2)
  SEQ1 <- seq(1, TT - width + 1, by = 1)
  SEQ2 <- lapply(SEQ1, function(x) x:(x + width - 1))
  SjmMatrix=matrix(0,nrow=pp,ncol=TT-width+1)
  for (uu in 1:(TT-width+1))
  { Temp_matrix=y_matrix[,SEQ2[[uu]]]
  SjmMatrix[,uu]=rowSums( Temp_matrix[,1:floor(m/2),drop=FALSE])-rowSums( Temp_matrix[,(floor(m/2)+1):width,drop=FALSE])
  }
  Critical=0
  for(gg in 1:lit)
  {  Noise=rnorm(TT+pp)
  rr=rep(0,ceiling(pp/rep_time))
  for( ss in 1: ceiling(pp/rep_time))
  { 
    temp_max=SjmMatrix[((((ss-1)*rep_time)+1):min(ss*rep_time,pp)),,drop=FALSE]%*%Noise[ss:(TT-width+ss)]##Noise[ss:(TT-m+ss)]
    rr[ss]=max(abs(temp_max))
  }
  Critical[gg]= max(rr)
  }
  norm_Const=sqrt(width*(TT-width+1))
  return(Critical/norm_Const)
}




fastlocal_constant=function(yy,bb,T=1/length(yy),endj=nn){
  nn=length(yy)
  
  a1=rep(0,nn)
  a2=rep(0,nn)
  a3=rep(0,nn)
  a4=rep(0,nn)
  a5=rep(0,nn)
  a6=rep(0,nn)
  
  b1=rep(0,nn)
  b2=rep(0,nn)
  b3=rep(0,nn)
  b4=rep(0,nn)
  b5=rep(0,nn)
  b6=rep(0,nn)
  
  TT=rep(0,nn)
  Tstart=nn*T
  
  TT[1]=T
  
  a1[1]=sum(yy*((1:nn)>=1)*((1:nn)<=floor(Tstart+nn*bb)) )
  a2[1]=sum((1:nn)*yy*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  a3[1]=sum((1:nn)^2*yy*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  a4[1]=sum((1:nn)^3*yy*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  
  
  b1[1]=sum(((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  b2[1]=sum((1:nn)*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  b3[1]=sum((1:nn)^2*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  b4[1]=sum((1:nn)^3*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  b5[1]=sum((1:nn)^4*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  Lend=-1*floor((nn*bb)-Tstart)
  Rend=floor(Tstart+(nn*bb))+1
  Endj=endj
  for(i in 2: Endj)
  {lend1=(Lend>=1)
  rend1= (Rend<=Endj)
  lend=Lend*lend1+(Lend<1)
  rend=Rend*rend1+(Rend>Endj)
  a1[i]=a1[i-1]-yy[lend]*(1)*lend1+yy[rend]*(1)*rend1
  a2[i]=a2[i-1]-yy[lend]*(lend)*lend1+yy[rend]*(rend)*rend1
  a3[i]=a3[i-1]-yy[lend]*(lend^2)*lend1+yy[rend]*(rend^2)*rend1
  a4[i]=a4[i-1]-yy[lend]*(lend^3)*lend1+yy[rend]*(rend^3)*rend1
  
  b1[i]=b1[i-1]-(1)*lend1+(1)*rend1
  b2[i]=b2[i-1]-(lend)*lend1+(rend)*rend1
  b3[i]=b3[i-1]-(lend^2)*lend1+(rend^2)*rend1
  b4[i]=b4[i-1]-(lend^3)*lend1+(rend^3)*rend1
  b5[i]=b5[i-1]-(lend^4)*lend1+(rend^4)*rend1
  
  T=T+1/nn
  Rend=Rend+1
  Lend=Lend+1
  TT[i]=T
  
  }
  
  Y_1=a1
  Y_2=bb^(-1)*(a2/nn-TT*a1)
  Y_3=bb^(-2)*(a3/nn^2-2*a2*TT/nn+TT^2*a1)
  Y_4=bb^(-3)*(a4/nn^3-3*a3*TT/nn^2+3*a2*TT^2/nn-TT^3*a1)
  
  
  
  
  R_1=b1
  R_2=bb^(-1)*(b2/nn-TT*b1)
  R_3=bb^(-2)*(b3/nn^2-2*b2*TT/nn+TT^2*b1)
  R_4=bb^(-3)*(b4/nn^3-3*b3*TT/nn^2+3*b2*TT^2/nn-TT^3*b1)
  R_5=bb^(-4)*(b5/nn^4-4*b4*TT/nn^3+6*b3*TT^2/nn^2-4*b2/nn*TT^3+TT^4*b1)
  
  
  K1=3/4*(R_1-R_3)/nn/bb
  
  K2=3/4*(R_2-R_4)/nn/bb
  
  K3=3/4*(R_3-R_5)/nn/bb
  
  
  Y1=3/4*(Y_1-Y_3)/nn/bb
  
  Y2=3/4*(Y_2-Y_4)/nn/bb
  
  coeff1=Y1/(K1)
  return(list(mean=coeff1[1:endj],time=TT[1:endj]))
  
}

Longrundim1_lconst=function(xx,m,tau,nn2=length(xx))
{nn=length(xx)
rolmean=wapply1(x=xx,width=m,by=1)
Diff=rolmean[1:(nn-2*m)]-rolmean[(m+1):(nn-m)]
tt=(1:nn2)/nn2
Q=rep(0,nn)

Q[(m+1):(nn-m)]=Diff^2*m/2
Longvar=fastlocal_constant(yy=Q,bb=tau,T=1/nn,endj=nn)$mean

ss=m
Longvar[1:ss]=Longvar[(ss+1)]
Longvar[(nn-ss+1):nn]=Longvar[(nn-ss)]
return(Longvar^0.5)
}
# ###New function######
# M1=function(u,t)
# {return((u+2*t)^2/2)
#   #return((1+u)*(6*(t-0.5)^2+1))
# }


KKernel=function(x)
{
  ss=0.75*(1-x^2)*(abs(x)<=1)
  return(ss)
}
HKernel=function(x)
{ 
  return((45/32-x^2*150/32+105*x^4/32)*(abs(x)<=1))
}

M_select=function(y_matrix, Mseq,rep_time)###Mseq is the set of m's
{
  pp=dim(y_matrix)[1]
  TT=dim(y_matrix)[2]
  MM=length(Mseq)
  SjmMatrix_tuning=matrix(0,nrow=pp,ncol=MM)
  for (um in 1:(MM))
  {
    width=2*floor(Mseq[um]/2)
    SEQ1 <- seq(1, TT - width + 1, by = 1)
    SEQ2 <- lapply(SEQ1, function(x) x:(x + width - 1))
    
    for (uu in 1:(TT-width+1))
    { Temp_matrix=y_matrix[,SEQ2[[uu]]]
    SjmMatrix_temp=rowSums( Temp_matrix[,1:floor(width/2),drop=FALSE])-rowSums( Temp_matrix[,floor(width/2):width,drop=FALSE])
    SjmMatrix_tuning[,um]= SjmMatrix_tuning[,um]+(SjmMatrix_temp)^2
    }
    SjmMatrix_tuning[,um]=  SjmMatrix_tuning[,um]/width/(TT-width+1)
  }
  zzz=0
  for (uu in 1:(MM-4))
  {
    zzz[uu]=mean(rowSds(SjmMatrix_tuning[,uu:(uu+4)]))
  }
  return(Mseq[(order(zzz)+2)[1]])
}


# 
# 
# LSAR_fts=function(uu,tp)
# { y=matrix(0,nrow=uu,ncol=length(tp))
# NN=100###Burnin
# Innov=rnorm(uu+NN)
# for (ss in 1:length(tp))
# {
#   for (i in 1:uu)
#   {ARcoeff=0.5*cos(pi/3*i/uu)-0.1*tp[ss]
#   innov=Innov[(i+NN):1] 
#   y[i,ss]=sum(innov*ARcoeff^(0:(i+NN-1)))
#   }
# }
# return(y)
# }
# 
# 
# LSAR1=function(nn)
# { y=rep(0,nn)
# NN=100###Burnin
# Innov=rnorm(nn+NN)
# for (i in 1:nn)
# {ARcoeff=0.5*cos(pi/3*i/nn)
# innov=Innov[(i+NN):1] 
# y[i]=sum(innov*ARcoeff^(0:(i+NN-1)))
# }
# return(y)
# }
# 
# LSARMA=function(nn)
# { y=rep(0,nn)
# NN=100###Burnin
# Innov1=rt(nn+NN+1,df=8)
# for (i in 1:nn)
# {ARcoeff=0.4*i/nn
# MAcoeff=0.3*(i/nn)^2
# innov=Innov1[(i+NN+1):2]-MAcoeff*Innov1[(i+NN):1] 
# y[i]=sum(innov*ARcoeff^(0:(i+NN-1)))
# }
# return(y)
# }

# M2=function(u,t)
# {
#   return((1+u^2)*(6*(t-0.5)^2*(1+1*(t>0.3))+1))
# }
# FData_simu_12=function(uu,nn,vv=1)
# { Time_Series=LSAR_fts(uu=uu,tp=(1:nn)/nn)/vv##now it is a matrix
# Signal_matrix=matrix(0,nrow=uu,ncol=nn)
# Data_matrix=matrix(0,nrow=uu,ncol=nn)
# for (ss in 1:uu)
# { temp=M2(u=ss/uu,t=1:nn/nn)
# Signal_matrix[ss,]=temp
# Data_matrix[ss,]=temp+Time_Series[ss,]*(1+0.5*sin(pi*(1:nn)/nn))##the second part is the basis which could be changed
# }
# return(list(Data=Data_matrix, Signal=Signal_matrix, Time_Series=Time_Series ))
# }
# 
# FData_simu_22=function(uu,nn,vv=1)
# { Time_Series1=LSAR1(uu)/vv
# Time_Series2=LSARMA(uu)/vv
# Signal_matrix=matrix(0,nrow=uu,ncol=nn)
# Data_matrix=matrix(0,nrow=uu,ncol=nn)
# zz=(1:nn)/nn
# for (ss in 1:uu)
# { temp=M2(u=ss/uu,t=1:nn/nn)
# Signal_matrix[ss,]=temp
# Data_matrix[ss,]=temp+Time_Series1[ss]*(2*zz-1)+Time_Series2[ss]*(6*zz^2-6*zz+1)##the second part is the basis which could be changed
# }
# return(list(Data=Data_matrix, Signal=Signal_matrix, Time_Series=cbind(Time_Series1,Time_Series2)  ))
# }
# 
# 
# 
# FData_simu=function(uu,nn,vv=1)
# { Time_Series=LSAR_fts(uu=uu,tp=(1:nn)/nn)/vv###This is a matrix
# Signal_matrix=matrix(0,nrow=uu,ncol=nn)
# Data_matrix=matrix(0,nrow=uu,ncol=nn)
# for (ss in 1:uu)
# { temp=M1(u=ss/uu,t=1:nn/nn)
# Signal_matrix[ss,]=temp
# Data_matrix[ss,]=temp+Time_Series[ss,]*(1+0.5*sin(pi*(1:nn)/nn))##the second part is the basis which could be changed
# }
# return(list(Data=Data_matrix, Signal=Signal_matrix, Time_Series=Time_Series ))
# }
# 
# FData_simu_2=function(uu,nn,vv=1)
# { Time_Series1=LSAR1(uu)/vv
# Time_Series2=LSARMA(uu)/vv
# Signal_matrix=matrix(0,nrow=uu,ncol=nn)
# Data_matrix=matrix(0,nrow=uu,ncol=nn)
# zz=(1:nn)/nn
# for (ss in 1:uu)
# { temp=M1(u=ss/uu,t=1:nn/nn)
# Signal_matrix[ss,]=temp
# Data_matrix[ss,]=temp+Time_Series1[ss]*(2*zz-1)+Time_Series2[ss]*(6*zz^2-6*zz+1)##the second part is the basis which could be changed
# }
# return(list(Data=Data_matrix, Signal=Signal_matrix, Time_Series=cbind(Time_Series1,Time_Series2)  ))
# }


Fix_u_SCB_matrix_matrix=function(Residual_Matrix,n,b)
{ startv=ceiling(n*u-n*b); endv=floor(n*u+n*b)
Kvector=HKernel((startv:endv/n-u)/b); Test_matrix_u=Residual_Matrix[startv:endv,]*Kvector
return(t(Test_matrix_u))
}



SCB_matrix=function(Residual_Matrix,b)
{ pp=dim(Residual_Matrix)[1]
n=pp
TT=dim(Residual_Matrix)[2]
upbound=2*ceiling(n*b)-1
Matrix_SCB=matrix(0,nrow=TT*(n-2*ceiling(n*b)+1),ncol=upbound)
residual_matrix=Residual_Matrix
for(ss in 1 : upbound)
{
  for(kk in 0: (n-2*ceiling(n*b)) ) 
  {
    Matrix_SCB[(kk*TT+1):(kk*TT+TT),ss]=residual_matrix[(ss+kk),]*HKernel(((ss)/n-ceiling(n*b)/n)/b)
  }
}
return(Matrix_SCB)
}


SCB_matrix_boundary=function(Residual_Matrix,b)
{ pp=dim(Residual_Matrix)[1]
n=pp
TT=dim(Residual_Matrix)[2]
upbound=ceiling(n*b)

Matrix_SCB_left=matrix(0,nrow=TT*(ceiling(n*b)-1),ncol=upbound)
residual_matrix=Residual_Matrix
for(ss in 1 : (upbound))
{ 
  for(kk in 0: (ceiling(n*b)-2) ) 
  { ggtemp=(1:n/n-(kk+1)/n)/b
  NWconstant=sum(SideKl(ggtemp))
  Matrix_SCB_left[(kk*TT+1):(kk*TT+TT),ss]=residual_matrix[(ss+kk),]*SideKl(((ss-1)/n)/b)/NWconstant
  }
}
Matrix_SCB_right=matrix(0,nrow=TT*(ceiling(n*b)),ncol=upbound)
for(ss in (n-2*upbound+2) : (n-upbound+1))
{
  for(kk in 0: (ceiling(n*b)-1)) 
  {ggtemp=(1:n/n-(n-upbound+kk+1)/n)/b
  NWconstant=sum(SideKr(ggtemp))
  Matrix_SCB_right[(kk*TT+1):(kk*TT+TT),ss-(n-2*(upbound)+2)+1]=residual_matrix[(ss+kk),]*SideKr(((ss)/n-(n-ceiling(n*b)+1)/n)/b)/NWconstant
  }
}
return(list(Matrix_SCB_left=Matrix_SCB_left,Matrix_SCB_right=Matrix_SCB_right))
}





#################################New Single Kernel################
# if(1<0)
# {
#   
#   SideKl=function(x)
#   {
#     return((-140*x^3+240*x^2-120*x+16)*(x>=0)*(x<=1))
#   }
#   SideKr=function(x){return(SideKl(-1*x))}
#   
#   
#   
#   SideKl=function(x)
#   {
#     return((-2520*x^3+4200*x^2-2100*x+300)*x*(1-x)*(x>=0)*(x<=1))
#   }
#   SideKr=function(x){return(SideKl(-x))}
#   
# }


SideKl=function(x)
{
  return((420*x^2-480*x+120)*x*(1-x)*(x>=0)*(x<=1))
}
SideKr=function(x){return(SideKl(-x))}
##################################################################
# Meanestimate=function(y,t,n,b)################NewMeanestimate#################
# {
#   zz=(1:n/n-t)/b
#   if(t<b){
#     yy=sum(SideKl(zz)*y)/sum(SideKl(zz))}
#   else if(t>1-b){
#     yy=sum(SideKr(zz)*y)/sum(SideKr(zz))}
#   else{yy=sum(HKernel(zz)*y)/sum(HKernel(zz))}
#   return(yy)
# }

Meanestimate=function(y,t,n,b)################NewMeanestimate#################
{
  zz=(1:n/n-t)/b
  yy=sum(HKernel(zz)*y)/sum(HKernel(zz)) 
  return(yy)
}

Curveestimate=function(y,b,tt=1:length(y)/length(y))############New Curveestimate#############
{ 
  n_y=length(y)
  esti=rep(0,length(tt))
  for (i in 1:length(tt))
  {
    esti[i]=Meanestimate(y=y,t=tt[i],n=n_y,b=b)
  }
  return(esti)
}




Gcvtrace=function(b,n)
{ DA=rep(0,n)
for(tt in 1:n)
{
  t=tt/n
  zz=(1:n/n-t)/b
  if(t<b){
    qq=(SideKl(0))/sum(SideKl(zz))}
  else if(t>1-b){
    qq=(SideKr(0))/sum(SideKr(zz))}
  else{qq=(HKernel(0))/sum(HKernel(zz))}
  DA[tt]=qq
}
return(DA)
}
###############################New MGCV#########
bandselect=function(yy,b,n)
{ res=rep(0,n)
for (tt in 1:n)
{t=tt/n
y_fit=Meanestimate(y=yy,t=t,b=b,n=n)
res[tt]=(yy-y_fit)[tt]
}
Est_trace=Gcvtrace(b,n)
Upp=sum(res^2)
Low=(1-sum(Est_trace)/n)^2
return(Upp/Low)
}

Meanbandselect=function(y,bb)
{ NN=dim(y)[2]
UU=dim(y)[1]
selectValue=0
for(u in 1:NN)
{
  selectValue[u]= bandselect(yy=y[,u],b=bb,n=UU)
}
return(max(selectValue))
}


##############################
##NN: the number of grids in function (could choose UU^(1/2); UU: number of functional time series (or Time of functional time series))
##New version: delete aa
GCVnorm=function(b,n)
{
  DAnorm=rep(0,n)
  for(tt in 1:n)
  {
    t=tt/n
    zz=(1:n/n-t)/b
    if(t<b){
      qq=sum(SideKl(zz))}
    else if(t>1-b){
      qq=sum(SideKr(zz))}
    else{qq=sum(HKernel(zz))}
    DAnorm[tt]=qq
  }
  return(DAnorm)
}
#DAnorm=GCVnorm(b=0.1,n=500)
Gcvtrace_usenorm=function(b=0.1,n=500,DAnorm)
{ DA=rep(0,n)
for(tt in 1:n)
{
  t=tt/n
  zz=(1:n/n-t)/b
  if(t<b){
    qq=(SideKl(0))}
  else if(t>1-b){
    qq=(SideKr(0))}else{qq=(HKernel(0))}
  DA[tt]=qq
}
return(DA/DAnorm)
}

bandselect_usenorm=function(yy,b,n,DAnorm)
{ res=rep(0,n)
for (tt in 1:n)
{t=tt/n
y_fit=Meanestimate(y=yy,t=t,b=b,n=n)
res[tt]=(yy-y_fit)[tt]
}
Est_trace=Gcvtrace_usenorm(b,n,DAnorm)
Upp=sum(res^2)
Low=(1-sum(Est_trace)/n)^2
return(Upp/Low)
}

Meanbandselect_usenorm=function(y,bb,DAnorm)
{ NN=dim(y)[2]
UU=dim(y)[1]
selectValue=0
for(u in 1:NN)
{
  selectValue[u]= bandselect_usenorm(yy=y[,u],b=bb,n=UU,DAnorm)
}
return(max(selectValue))
}


Meanbandselect_usenorm=function(y,bb,DAnorm)
{ NN=dim(y)[2]
UU=dim(y)[1]
selectValue=0
for(u in 1:NN)
{
  selectValue[u]= bandselect_usenorm(yy=y[,u],b=bb,n=UU,DAnorm)
}
return(max(selectValue))
}

round_up_2 <- function(x) {
  ceiling(x * 100) / 100
}

Wu = function(Y,n,N,alpha=c(0.01,0.05,0.1,0.2),
              B_error=0.15,boots=500,
              Mseq=seq(15,30,by=1), band.choose='F',bseq=seq(0.1,0.3,by=0.02)){
  
  DAnorm=matrix(0,nrow=length(bseq), ncol=n)
  
  for(bs in 1:length(bseq)){
    DAnorm[bs,]=GCVnorm(b=bseq[bs],n=n)
  }
  
  GCV=0
  if(band.choose=='T'){
    for(bs in 1:length(bseq)){
      GCV[bs]=Meanbandselect_usenorm(y=Y,bb=bseq[bs],DAnorm=DAnorm[bs,])
    }
    berror=bseq[order(GCV)[1]]
    bb=berror
  }else{
    # berror=round_up_2(0.6*n^(-1/5)) 0.75*n^(-1/8)
    # bb=round_up_2(0.6*n^(-1/5)) n^(-1/6) 0.85*n^(-1/7)
    berror=round_up_2(n^(-1/6)) 
    bb=round_up_2(n^(-1/6)) 
  }
  
  Estimat_FDA=matrix(0,nrow=n,ncol=N)
  Estimat_Lrv=matrix(0,nrow=n,ncol=N)
  for(u in 1:N){
    Estimat_FDA[,u]= Curveestimate(Y[,u],b=bb)
    Estimat_Lrv[,u]= Longrundim1_lconst(xx=Y[,u],m=round(n^(2/7)),tau=n^(-1/7),nn2=n)
  }
  
  lattice = seq(1/N,1,1/N)
  lattice_time = seq(1/n,1,1/n)
  mean_function = function(tau,x){
    return((1+tau)*(6*(x-1/2)^2+1))
  }
  m.eval = matrix(mean_function(rep(lattice_time,each = N), rep(lattice,n)),
                  n,N,byrow=T)
  
  
  # max_dev=max(abs((Estimat_FDA - m.eval )/
  #                   Estimat_Lrv ))*sqrt(bb*n)/sqrt(2)
  # 
  # Residual_Matrix=Y-Estimat_FDA
  # Matrix_SCB=SCB_matrix(Residual_Matrix/Estimat_Lrv,b=bb)
  # mm = 30
  # # mm=M_select(y_matrix= Matrix_SCB, Mseq=Mseq,rep_time=N)
  # Critical=HDtest((Matrix_SCB),m=mm,lit=boots,rep_time =N)
  # return(c(max_dev<quantile(Critical,1-alpha)))
  
  ############################################################
  # Critical values
  ############################################################
  Residual_Matrix = Y - Estimat_FDA
  Matrix_SCB = SCB_matrix(Residual_Matrix / Estimat_Lrv, b = bb)
  
  mm = 30
  Critical = HDtest(Matrix_SCB, m = mm, lit = boots, rep_time = N)
  
  Q = as.numeric(
    quantile(Critical, probs = 1 - alpha, names = FALSE)
  )
  
  ############################################################
  # Full domain [0,1]^2
  ############################################################
  
  # standardized maximal deviation over [0,1]^2
  max_dev_full =
    max(abs((Estimat_FDA - m.eval) / Estimat_Lrv)) *
    sqrt(bb * n) / sqrt(2)
  
  coverage_full = as.numeric(max_dev_full < Q)
  
  # Half-width of the SCR at every grid point:
  # Q * Estimat_Lrv * sqrt(2 / (n * h)), where h = bb
  #
  # For each alpha, the volume is
  # integral_[0,1]^2 (Upper - Lower) d tau d x.
  volume_full = sapply(Q, function(q) {
    width_matrix =
      2 * q * Estimat_Lrv * sqrt(2 / (bb * n))
    
    # area([0,1]^2) = 1
    mean(width_matrix)
  })
  
  ############################################################
  # Interior domain [h,1-h] x [0,1]
  ############################################################
  
  # This is the same interior grid used by SCB_matrix()
  left_index  = ceiling(bb * n)
  right_index = n - ceiling(bb * n)
  
  interior_index = left_index:right_index
  
  max_dev_interior =
    max(
      abs(
        (Estimat_FDA[interior_index, , drop = FALSE] -
           m.eval[interior_index, , drop = FALSE]) /
          Estimat_Lrv[interior_index, , drop = FALSE]
      )
    ) * sqrt(bb * n) / sqrt(2)
  
  coverage_interior = as.numeric(max_dev_interior < Q)
  
  # Volume over [h,1-h] x [0,1].
  # Its area is 1 - 2h = 1 - 2*bb.
  volume_interior = sapply(Q, function(q) {
    width_matrix =
      2 * q *
      Estimat_Lrv[interior_index, , drop = FALSE] *
      sqrt(2 / (bb * n))
    
    mean(width_matrix) * (1 - 2 * bb)
  })
  
  return(c(
    coverage_full,
    volume_full,
    coverage_interior,
    volume_interior
  ))
}
#Wu(Y,n,N)

