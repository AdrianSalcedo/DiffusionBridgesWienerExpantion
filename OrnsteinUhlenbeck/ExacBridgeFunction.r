### Exact bridge
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

