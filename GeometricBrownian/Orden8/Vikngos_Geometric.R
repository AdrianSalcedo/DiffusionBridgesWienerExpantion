path_GBM=function(a,delta,n,alpha,sigma)
{
  X=numeric(n+1)
  X[1]=a
  inc_brow=sqrt(delta)*rnorm(n,0,1)
  for(i in 1:n)
  {
    X[i+1]=X[i]+alpha*X[i]*delta+sigma*X[i]*inc_brow[i]+(1/2)*sigma*X[i]*sigma*(inc_brow[i]^2-delta)
  }
  return(X)
}

M_bridges_MM_GBM=function(M,a,b,delta,n,alpha,sigma)
{
  MOU=matrix(0,nrow=M,ncol=(n+1))
  for(i in 1:M){
    #print(i)
    MOU[i,]=Bridge_MM_GBM(a,b,delta,n,alpha,sigma)
  }
  return(MOU)
}

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





