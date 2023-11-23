

###Exact OU bridge from a to b
#dXt=-thetaOU_Xtdt+sigmaOU_dW
alpha=0.2
sigma=0.3
delta=1/1000
n=1000
nb=16
a=0.2
b=0.3
# number of bridges
M=1000

Xms=t(as.matrix(read.csv(file = 'propagator.csv'))[,-1])
BChaos=Bri_Chaos_GBM(Xms,a,b,alpha,sigma,n,delta,nb,M)
plot( TiempoC,M_MM[1,],type="l")
lines(TiempoC,BChaos[2,],col="red")
lines(TiempoC,BChaos[3,],col="blue")
lines(TiempoC,BChaos[4,],col="green")
lines(TiempoC,BChaos[5,],col="orange")
lines(TiempoC,M_MM[5,],col="orange")

obs=800
qqplot(BChaos[,obs],M_MM[,obs],ylab="Exact bridge",xlab="Wiener chaos approximation bridge", cex.lab=0.5,main="Quantile-Quantile comparison", cex.main=0.5)
abline(0,1)


##########
Bri_Chaos_GBM=function(Xms,a,b,alpha,sigma,n,delta,nb,M)
{
  
  
  TF=delta*n
  TiempoC<-seq(0, TF, length.out = n+1)
  TL<-length(TiempoC)
  
  OUchaos=NULL
  
  M_Y=mat.or.vec(M,TL)
  
  
  Yms=mat.or.vec(117,TL)

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
    
  
    Xis=rbind(fXis_1(100,TiempoC),fXis_2_j(c(1,2),TiempoC),fXis_2_j(c(1,3),TiempoC),fXis_2_j(c(2,2),TiempoC),fXis_2(3,TiempoC),fXis_3_12(c(1,2),TiempoC),fXis_3_21(c(1,2),TiempoC),fXis_3(2,TiempoC),fXis_4(2,TiempoC),fXis_5(2,TiempoC),fXis_6(2,TiempoC))
   
 
    M_Y[k,]=apply(Yms[2:117,]*Xis,2,sum)+Yms[1,]
    
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
    integral[i]=(Xms[i]-Xms[i-1])*(path[i-1]+path[i])/2.0
  }
  integral=cumsum(integral)
  
  
  return(integral)
  
}










### Calcula las Xis

fXis_1=function(nb,TiempoC)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  BM=My_SBM(nb,TiempoC)
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(2/TF)*Ints
  
  
  return(Xis)
}

fXis_2=function(nb,TiempoC)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  BM=My_SBM(nb,TiempoC)
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(2)*H2(sqrt(2/TF)*Ints)
  
  
  return(Xis)
}

fXis_3=function(nb,TiempoC)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  BM=My_SBM(nb,TiempoC)
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(3))*H3(sqrt(2/TF)*Ints)
  
  
  return(Xis)
}

fXis_4=function(nb,TiempoC)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  BM=My_SBM(nb,TiempoC)
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(4))*H4(sqrt(2/TF)*Ints)
  
  
  return(Xis)
}

fXis_5=function(nb,TiempoC)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  BM=My_SBM(nb,TiempoC)
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(5))*H5(sqrt(2/TF)*Ints)
  
  
  return(Xis)
}

fXis_6=function(nb,TiempoC)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  BM=My_SBM(nb,TiempoC)
  
  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }
  
  Xis=sqrt(factorial(6))*H6(sqrt(2/TF)*Ints)
  
  
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
    integral[i]=(senos[i]*W[i]-senos[i-1]*W[i-1])-(delta*j*pi/TF)*(W[i]*cosenos[i]+W[i-1]*cosenos[i-1])/2
  }
  
  integral=cumsum(integral)
  return(integral)
  
}

#calcula senos
fsenos=function(n,TiempoC)
{
  TL<-length(TiempoC)
  TF=tail(TiempoC,1)
  senos<-mat.or.vec(n,TL) 
  for(i in 1:n){   
    senos[i,]=sin(i*pi*TiempoC/TF)
    
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
    cosenos[i,]=cos(i*pi*TiempoC/TF)
    
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
    senos[i,]=sin(js[i]*pi*TiempoC/TF)
    
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
    cosenos[i,]=cos(js[i]*pi*TiempoC/TF)
    
  }
  return(cosenos)
  
}
####hermite polynomial

H2=function(x)
{
  h=(x^2-1)/sqrt(2)
  return(h)
}

H3=function(x)
{
  h=(x*(x^2-3))/sqrt(factorial(3))
  return(h)
}

H4=function(x)
{
  h=(x^4-6*x^2+3)/sqrt(factorial(4))
  return(h)
}

H5=function(x)
{
  h=(x*(x^4-10*x^2+15))/sqrt(factorial(5))
  return(h)
}

H6=function(x)
{
  h=(x^6-15*x^4+45*x^2-15)/sqrt(factorial(6))
  return(h)
}


fXis_2_j=function(js,TiempoC)
{
  nb=length(js)
  senos=fsenos_j(js,TiempoC)
  cosenos=fcosenos_j(js,TiempoC)
  
  TF=tail(TiempoC,1)
  BM=My_SBM(nb,TiempoC)
  
  Ints=BM
  Xis=1
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],js[j],TiempoC,BM[j,])
    Xis=Xis*sqrt(2/TF)*Ints[j,]
  }
  
 
  
  
  return(Xis)
}


fXis_3_12=function(js,TiempoC)
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
  
  Xis=sqrt(2)*sqrt(2/TF)*Ints[1,]*H2(sqrt(2/TF)*Ints[2,])
  
  
  return(Xis)
}


fXis_3_21=function(js,TiempoC)
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
  
  Xis=sqrt(2)*H2(sqrt(2/TF)*Ints[1,])*sqrt(2/TF)*Ints[2,]
  
  
  return(Xis)
}





