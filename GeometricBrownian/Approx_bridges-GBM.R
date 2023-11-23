
###Exact OU bridge from a to b
#dXt=-thetaOU_Xtdt+sigmaOU_dW
alpha=0.2
sigma=0.3
delta=1/1000
n=1000
nbri=1000
a=0.2
b=0.3
# number of bridges
M=1000

GBM=path_GBM(M,delta,n,alpha,sigma)

plot(MGBM[1,],type="l")

MGBM=Mpath_GBM(M,a,delta,n,alpha,sigma)


M_MM1=M_bridges_MM_GBM(M,a,b,delta,n,alpha,sigma)

datos_Chaos_GBM=as.matrix(read.csv(file = '/Users/fernandobaltazar-larios/Dropbox/Working_papers/Pancho/wiener_chaos/GBM/Datos.csv'))[,-1]
data=as.matrix(read.csv(file = '/Users/fernandobaltazarlarios/Dropbox/Working_papers/Pancho/wiener_chaos/GBM/Datos2.csv'))[,-1]

plot(MGBM[3,],type="l")
lines([2,])
lines(M_MM[1,])
lines(M_MM[4,])
lines(M_MM[5,])
plot(datos_Chaos_GBM[,1],type="l",col="purple")
lines(data[,2],col="red")
lines(data[,3],col="blue")
lines(data[,4],col="green")
lines(data[,5],col="yellow")


qqplot(MGBM[,900],data[900,])
abline(0,1)




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
