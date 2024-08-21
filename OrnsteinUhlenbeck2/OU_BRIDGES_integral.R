
####OU_bridge####
## Parametros del OU
thetaOU<-0.5
sigmaOU<-1.0


delta=1/1000
### Parametros del puente :
a<-0

b<-0

##numero de brownianos
nb=1000
#numero de puentes
M=1000
#numero de puntos por puente
n=1/delta
#ejemplo


exa=M_bridges_Exact_OU(M,a,b,delta,n,thetaOU,sigmaOU)
chaos=BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,nb,M)

obs=n/2

qqplot(exa[,obs],chaos[,obs])
abline(0,1)





BChaos1=BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,1,M)
BChaos10=BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,10,M)
BChaos100=BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,100,M)
BChaos500=BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,500,M)
BChaos1000=BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,1000,M)
BChaos1000_10=BridgeChaos(a,b,thetaOU,sigmaOU,n,1/10,nb,M)
BChaos1000_100=BridgeChaos(a,b,thetaOU,sigmaOU,n,1/100,nb,M)
BChaos1000_1000=BridgeChaos(a,b,thetaOU,sigmaOU,n,1/1000,nb,M)
dev.new(4,4)
plot( TiempoC,BChaos1[1,],type="l",main="Ornstein Uhlenbeck bridges",xlab="Time",cex.lab=0.5,ylab="Diffusion bridges",cex.main=0.5,ylim=c(-3,3))
lines(TiempoC,BChaos10[1,],col="red")
lines(TiempoC,BChaos100[1,],col="blue")
lines(TiempoC,BChaos500[1,],col="green")
lines(TiempoC,BChaos1000[1,],col="orange")
legend(x = "topright",  
       cex = 0.5,# Position
       legend = c("1", "10","100","500","1000"),  # Legend texts
       lty = c(1, 1,1,1,1),           # Line types
       col = c("black","red","blue","green","orange"),           # Line colors
       lwd = 1) 
lines(M_OUE1[3,],col="blue")
lines(M_OUE1[4,],col="green")

EB_OU=M_bridges_MM_OU(M,a,b,delta,n,thetaOU,sigmaOU)


obs=500

qqplot(BChaos1000[,obs],EB_OU,ylab="Exact bridge",xlab="Wiener chaos approximation bridge", cex.lab=0.5,main="Quantile-Quantile comparison", cex.main=0.5)
abline(0,1)
qqplot(BChaos1000_100[,obs],M_OUE100[,obs],ylab="Exact bridge",xlab="Wiener chaos approximation bridge", cex.lab=0.5,main="Quantile-Quantile comparison", cex.main=0.5)
abline(0,1)
qqplot(BChaos1000_10[,obs],M_OUE10[,obs],ylab="Exact bridge",xlab="Wiener chaos approximation bridge", cex.lab=0.5,main="Quantile-Quantile comparison", cex.main=0.5)
abline(0,1)


plot(density(BChaos[,obs]))
lines(density(M_OUE1[,obs]),col="red")


##########
BridgeChaos=function(a,b,thetaOU,sigmaOU,n,delta,nb,M)
{
  
  
  TF=delta*n
  TiempoC<-seq(0, TF, length.out = n+1)
  TL<-length(TiempoC)
  
  OUchaos=NULL

M_Y=mat.or.vec(M,TL)
#X_iniOU<-rnorm(1,0,sqrt(sigmaOU^2/(2*thetaOU)))
#### propagator for |m|=0
X_0 <- a*exp(-thetaOU*TiempoC)
Y_0=X_0
int1=c((TiempoC[TL]-TiempoC[1:(TL-1)])*Integrate_Xms(1/(TF-TiempoC[-TL]),X_0),0)

for(i in 1:TL)
{
  Y_0[i]= a+(b-a)*TiempoC[i]/TF+int1[i]
}




#### propagator for |m|=1
cosenos=fcosenos(nb,TiempoC)
senos=fsenos(nb,TiempoC)
Y_1=NULL
for(i in 1:nb){
  
  X_1=sqrt(2/TF)*sigmaOU*(i*pi*exp(-thetaOU*TiempoC)-i*pi*cosenos[i,]+thetaOU*TF*senos[i,])/((i*pi)^2+(thetaOU*TF)^2)
  Y_1=rbind(Y_1,c((TiempoC[TL]-TiempoC[1:(TL-1)])*Integrate_Xms(1/(TF-TiempoC[-TL]),X_1),0))
  
}


for (k in 1:M){
  print(k)
  

Xis=fXis(nb,TiempoC)


M_Y[k,]=apply(Y_1*Xis[,TL],2,sum)+Y_0

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
      MBM[i,j]<- MBM[i,j-1]+rnorm(1,0,1)
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

fXis=function(nb,TiempoC)
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



#### Exact bridge
Exact_OU=function(a,b,delta,n,thetaOU,sigmaOU)
{
  Z=numeric(n+1)
  
  X=path_X(a,delta,n,thetaOU,sigmaOU)
  ts=seq(0,delta*n,by=delta)
  for(i in 1:(n+1))
  {
    Z[i]=X[i]+(b-X[n+1])*(exp(thetaOU*ts[i])-exp(-thetaOU*ts[i]))/(exp(thetaOU*ts[n+1])-exp(-thetaOU*ts[n+1]))
  }
  return(Z)
}



M_bridges_Exact_OU=function(M,a,b,delta,n,thetaOU,sigmaOU)
{
  MOU=matrix(0,nrow=M,ncol=(n+1))
  for(i in 1:M){
    MOU[i,]=Exact_OU(a,b,delta,n,thetaOU,sigmaOU)
  }
  return(MOU)
}

### path OU
path_X=function(a,delta,n,thetaOU,sigmaOU)
{
  X=numeric(n+1)
  X[1]=a
  W=rnorm(n,0,sqrt((sigmaOU^2)*(1-exp(-2*thetaOU*delta))/(2*thetaOU)))
  for(i in 1:n)
  {
    X[i+1]=exp(-thetaOU*delta)*X[i]+W[i] 
  }
  return(X)
}

###AQUI

##### puentes M&M

Bridge_MM=function(a,b,delta,n,theta,sigma)
{
  X=path_OU(a,delta,n,theta,sigma)
  bridge=numeric(n+1)
  ban=0
  while(ban==0){
    Y=rev(path_OU(b,delta,n,theta,sigma))
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


M_bridges_MM_OU=function(M,a,b,delta,n,theta,sigma)
{
  MOU=matrix(0,nrow=M,ncol=(n+1))
  for(i in 1:M){
    MOU[i,]=Bridge_MM(a,b,delta,n,theta,sigma)
  }
  return(MOU)
}



path_OU=function(a,delta_bri,nbri,theta,sigma)
{
  X=numeric(nbri+1)
  X[1]=a
  for(i in 1:nbri)
  {
    X[i+1]=X[i]-theta*X[i]*delta_bri+sigma*rnorm(1,0,sqrt(delta_bri))
  }
  return(X)
  
}

