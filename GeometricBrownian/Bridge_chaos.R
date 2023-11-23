###Exact OU bridge from a to b
#dXt=-thetaOU_Xtdt+sigmaOU_dW
alpha=0.2
sigma=0.3
delta=1/100
n=1000
nb=12005
a=0.2
b=0.3
# number of bridges
M=100

Xms=t(as.matrix(read.csv(file = 'C:/Users/Usuario1/Desktop/Propagador/Propagator_order12_1000MB.csv'))[,-1])
dim(Xms)
Xms_aux=Xms[1:1001,]
BGBM=BridgesGBM(Xms,a,b,alpha,sigma,n,delta,nb,M)
M_MM=BGBM$MB
BChaos_II=BGBM$CB
plot(density(M_Y[,500]),col="red")
plot(density(BChaos_II[,500]),col="red")
lines(density(M_MM[,500]))
dev.new()
plot(TiempoC,M_MM[1,],type="l")
lines(TiempoC,BChaos_II[1,],col="red")
lines(TiempoC,BChaos_II[3,],col="blue")
lines(TiempoC,BChaos_II[4,],col="green")
lines(TiempoC,BChaos_II[5,],col="orange")
lines(TiempoC,M_MM_II[5,],col="orange")

obs=500
qqplot(BChaos_II[,obs],M_MM[,obs],ylab="Exact bridge",xlab="Wiener chaos approximation bridge", cex.lab=0.5,main="Quantile-Quantile comparison", cex.main=0.5)
abline(0,1)

plot(density(M_MM[,obs]),col="red")
lines(density(BChaos[,obs]))
qqplot(BChaos_1_1[,obs],BChaos_1_1000[,obs],ylab="Exact bridge",xlab="Wiener chaos approximation bridge", cex.lab=0.5,main="Quantile-Quantile comparison", cex.main=0.5)
abline(0,1)

ks.test(BChaos[,obs],M_MM[,obs])


##########
Bri_Chaos_GBM=function(Xms,a,b,alpha,sigma,n,delta,nb,M)
{
  
  
  TF=delta*n
  TiempoC<-seq(0, TF, length.out = n+1)
  TL<-length(TiempoC)

  
  M_Y=mat.or.vec(M,TL)
  
  
  Yms=mat.or.vec(nb+1,TL)

  #### propagator for |m|=0
 
 
  int1=c((TiempoC[TL]-TiempoC[1:(TL-1)])*Integrate_Xms(1/(TF-TiempoC[-TL]),Xms[1,]),0)
  
  for(i in 1:TL)
  {
    Yms[1,i]= a+(b-a)*TiempoC[i]/TF+int1[i]
  }
  
  
  
  
  #### propagator for |m|>0
  
  for(i in 1:nb){
    
      Yms[i+1,]=c((TiempoC[TL]-TiempoC[1:(TL-1)])*Integrate_Xms(1/(TF-TiempoC[-TL]),Xms[1+i,]),0)
  }
  
  
  for (k in 1:M){
    print(k)
    BM=My_SBM(1000,TiempoC)
    Xis=rbind(fXis_1(1000,TiempoC,BM),fXis_2_j(c(1,2),TiempoC,BM[1:2,]),fXis_2_j(c(1,3),TiempoC,BM[1:2,]),fXis_2_j(c(2,3),TiempoC,BM[1:2,]),fXis_2(1000,TiempoC,BM),fXis_3_12(c(1,2),TiempoC,BM[1:2,]),fXis_3_21(c(1,2),TiempoC,BM[1:2,]),fXis_3(1000,TiempoC,BM),fXis_4(1000,TiempoC,BM),fXis_5(1000,TiempoC,BM),fXis_6(1000,TiempoC,BM),fXis_7(1000,TiempoC,BM),fXis_8(1000,TiempoC,BM),fXis_9(1000,TiempoC,BM),fXis_10(1000,TiempoC,BM),fXis_11(1000,TiempoC,BM),fXis_12(1000,TiempoC,BM))
    M_Y[k,]=apply(Yms[2:(nb+1),]*Xis,2,sum)+Yms[1,]
  }
  return(M_Y)
}

#generando Brownianos estandar
My_SBM=function(n,TiempoC)
{
  MBM<-mat.or.vec(n,length(TiempoC)) 
  l=length(TiempoC)
  delta=TiempoC[2]-TiempoC[1]
  for(i in 1:n){ 
    for (j in 2:l)
    {
      MBM[i,j]<- MBM[i,j-1]+rnorm(1,0,sqrt(delta))
    }
  }
  return(MBM)
  
}
#integral con respecto a Xms
Integrate_Xms=function(path,Xms)
{
  npoints=length(path)
  integral=numeric(npoints)
  for( i in 2:npoints)
  {
    integral[i]=(Xms[i]-Xms[i-1])*(path[i-1]+path[i])/2
  }
  integral=cumsum(integral)
  return(integral)
  
}

### Calcula las Xis
fXis_1=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=Ints
  
  
  return(Xis)
}

fXis_2=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)

  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(2)*H2(Ints)
  
  
  return(Xis)
}

fXis_3=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
 
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(3))*H3(Ints)
  
  
  return(Xis)
}

fXis_4=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(4))*H4(Ints)
  
  
  return(Xis)
}

fXis_5=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)

  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(5))*H5(Ints)
  
  
  return(Xis)
}

fXis_6=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
 
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(6))*H6(Ints)
  
  
  return(Xis)
}



fXis_7=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
 
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(7))*H7(Ints)
  
  
  return(Xis)
}


fXis_8=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)

  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(8))*H8(Ints)
  
  
  return(Xis)
}


fXis_9=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
 
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(9))*H9(Ints)
  
  
  return(Xis)
}


fXis_10=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
 
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(10))*H10(Ints)
  
  
  return(Xis)
}


fXis_11=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)

  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(11))*H11(Ints)
  
  
  return(Xis)
}


fXis_12=function(nb,TiempoC,BM)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(12))*H12(Ints)
  
  
  return(Xis)
}
#integral de ito por trapecio acumulado
Int_Trap_Acum_2=function(senos,cosenos,j,t,W)
{
  npoints=length(senos)
  integral=numeric(npoints)
  
  delta=t[2]-t[1]
  TF=tail(t,1)
  for( i in 2:npoints)
  {
    integral[i]=(senos[i]*W[i]-senos[i-1]*W[i-1])-(delta*j*pi)*(W[i-1]*cosenos[i-1]+W[i-1]*cosenos[i-1])/(2*TF)
  }
  
  integral[]=sum(integral)
  return(integral)
  
}

#calcula senos
fsenos=function(n,TiempoC)
{
  TL<-length(TiempoC)
  TF=tail(TiempoC,1)
  senos<-mat.or.vec(n,TL) 
  for(i in 1:n){   
    senos[i,]=sqrt(2/TF)*sin(i*pi*TiempoC/TF)
    
  }
  return(senos)
  
}

#calcula cosenos
fcosenos=function(n,TiempoC)
{
  TL<-length(TiempoC)
  TF=tail(TiempoC,1)
  cosenos<-mat.or.vec(n,TL) 
  for(i in 1:n){   
    cosenos[i,]=sqrt(2/TF)*cos(i*pi*TiempoC/TF)
    
  }
  return(cosenos)
  
}


#calcula senos
fsenos_j=function(js,TiempoC)
{
  TL<-length(TiempoC)
  TF=tail(TiempoC,1)
  n=length(js)
  senos<-mat.or.vec(n,TL) 
  for(i in 1:n){   
    senos[i,]=sqrt(2/TF)*sin(js[i]*pi*TiempoC/TF)
    
  }
  return(senos)
  
}

#calcula cosenos
fcosenos_j=function(js,TiempoC)
{
  TL<-length(TiempoC)
  TF=tail(TiempoC,1)
  n=length(js)
  cosenos<-mat.or.vec(n,TL) 
  for(i in 1:n){   
    cosenos[i,]=sqrt(2/TF)*cos(js[i]*pi*TiempoC/TF)
    
  }
  return(cosenos)
  
}
####hermite polynomial

H2=function(x)
{
  h=((x^2)-1)/sqrt(2)
  return(h)
}

H3=function(x)
{
  h=(x*((x^2)-3))/sqrt(factorial(3))
  return(h)
}

H4=function(x)
{
  h=(x^4-6*(x^2)+3)/sqrt(factorial(4))
  return(h)
}

H5=function(x)
{
  h=(x*((x^4)-10*(x^2)+15))/sqrt(factorial(5))
  return(h)
}

H6=function(x)
{
  h=((x^6)-15*(x^4)+45*(x^2)-15)/sqrt(factorial(6))
  return(h)
}


H7=function(x)
{
  h=(x*((x^6)-21*(x^4)+105*(x^2)-105))/sqrt(factorial(7))
  return(h)
}

H8=function(x)
{
  h=((x^8)-28*(x^6)+210*(x^4)-420*(x^2)+105)/sqrt(factorial(8))
  return(h)
}

H9=function(x)
{
  h=(x*((x^8)-36*(x^6)+378*(x^4)-1260*(x^2)+945))/sqrt(factorial(9))
  return(h)
}
H10=function(x)
{
  h=((x^10)-45*(x^8)+630*(x^6)-3150*(x^4)+4725*(x^2)-945)/sqrt(factorial(10))
return(h)
}
H11=function(x)
{
  h=(x*(((x^10)-55*(x^8)+990*(x^6)-6930*(x^4)+1735*(x^2)-10395)))/sqrt(factorial(11))
return(h)
}

H12=function(x)
{
  h=((x^12)-66*(x^10)+1485*(x^8)-13860*(x^6)+51975*(x^4)-62370*(x^2)+10395)/sqrt(factorial(12))
return(h)
}

fXis_2_j=function(js,TiempoC,BM)
{
  nb=length(js)
  senos=fsenos_j(js,TiempoC)
  cosenos=fcosenos_j(js,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  Xis=1
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],js[j],TiempoC,BM[j,])
    Xis=Xis*Ints[j,]
  }
  return(Xis)
}


fXis_3_12=function(js,TiempoC,BM)
{
  nb=length(js)
  senos=fsenos_j(js,TiempoC)
  cosenos=fcosenos_j(js,TiempoC)
  
  TF=tail(TiempoC,1)
  BM=My_SBM(nb,TiempoC)
  
  Ints=BM
  
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],js[j],TiempoC,BM[j,])
    
  }
  Xis=Ints[1,]*H2(Ints[2,])
  return(Xis)
}


fXis_3_21=function(js,TiempoC,BM)
{
  nb=length(js)
  senos=fsenos_j(js,TiempoC)
  cosenos=fcosenos_j(js,TiempoC)
  
  TF=tail(TiempoC,1)
  Ints=BM
  
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],js[j],TiempoC,BM[j,])
    
  }
  Xis=sqrt(2)*H2(Ints[1,])*Ints[2,]
  return(Xis)
}

####Radon-Nykodim

RNBGM=function(TiempoC,b,sigma,BM,Ys)
{
  n=length(TiempoC)
  path1=(b-Ys)/((tail(TiempoC,1)-TiempoC)*sigma*Ys)
   
  I1=Integrate_Xms(path1,BM)
  I2=Integrate_Xms(path1^2,TiempoC)
  I1[n]=I2[n]=0
  RN=exp(I1-I2/2)
  
}
##### generate a path of GBM with milstein scheme
# a= initial value
# delta= discretization
# n=number of points
# alpha= drift parameter
# sigma=diffusion parameter
path_GBM=function(a,delta,n,alpha,sigma)
{
  X=numeric(n+1)
  X[1]=a
  inc_brow=rnorm(n,0,sqrt(delta))
  for(i in 1:n)
  {
    X[i+1]=X[i]+alpha*X[i]*delta+sigma*X[i]*inc_brow[i]+(1/2)*sigma*X[i]*sigma*(inc_brow[i]^2-delta)
  }
  return(X)
  
}

##### puentes M&M

Bridge_MM_GBM=function(a,b,delta,n,alpha,sigma)
{
  X=path_GBM(a,delta,n,alpha,sigma)
  bridge=numeric(n+1)
  ban=0
  while(ban==0){
    Y=rev(path_GBM(b,delta,n,alpha,sigma))
    if(X[1]<=Y[1])
    {
      
      for(i in 2:(n+1))
      {
        
        if(X[i]>Y[i])
        {
          
          bridge[1:(i-1)]=X[1:(i-1)]
          bridge[i:(n+1)]=Y[i:(n+1)]
          ban=1
          break
        }
      }
    }
    else{
      
      for(i in 2:(n+1))
      {
        if(X[i]<Y[i])
        {
          
          bridge[1:(i-1)]=X[1:(i-1)]
          bridge[i:(n+1)]=Y[i:(n+1)]
          ban=1
          break
        }
      }
      
    }
  }
  
  return(bridge)
  
}


M_bridges_MM_GBM=function(M,a,b,delta,n,alpha,sigma)
{
  MOU=matrix(0,nrow=M,ncol=(n+1))
  for(i in 1:M){
    print(i)
    MOU[i,]=Bridge_MM_GBM(a,b,delta,n,alpha,sigma)
  }
  return(MOU)
}

#generate M paths 
Mpath_GBM=function(M,a,delta,n,alpha,sigma)
{
  MX=matrix(0,nrow = M,ncol=(n+1))
  
  for(i in 1:M)
  {
    MX[i,]=path_GBM(a,delta,n,alpha,sigma)
  }
  return(MX)
  
}

BridgesGBM=function(Xms,a,b,alpha,sigma,n,delta,nb,M){
  

Chaos_bridges=Bri_Chaos_GBM(Xms,a,b,alpha,sigma,n,delta,nb,M)

  MM_bridges=M_bridges_MM_GBM(M,a,b,delta,n,alpha,sigma)
  
 
  return(list(CB=Chaos_bridges,MB=MM_bridges)) 
}
