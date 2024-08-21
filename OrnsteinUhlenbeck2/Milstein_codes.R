Milstein_GBM <- function(x0,xf,alpha,sigma,n,delta,TiempoC)
{
  Sol = matrix(0,1,n+1)
  Wt = matrix(0,1,n+1)
  for (i in 2:(n+1)){
    Wt[i] = Wt[i-1] + sqrt(delta)*rnorm(1,0,1)
  }
  Sol[1] = x0
  for(k in 2:(n+1))
  {
    DeltaW = Wt[k]-Wt[k-1]
    Sol[k] = Sol[k-1] +(alpha*Sol[k-1]+Sol[k-1]*(log(xf)-log(Sol[k-1])-(alpha-(sigma^2)/2)*(TF-TiempoC[k-1]))/(TF-TiempoC[k-1]))*delta +sigma*Sol[k-1]*DeltaW +(1/2)*(sigma^2)*Sol[k-1]*(DeltaW^2-delta)
  }
  return(Sol)
}

Milstein_OU <- function(M,x0,xf,alpha,sigma,n,delta,TiempoC)
{
  Chau_OU <- matrix(0,nrow=M,n+1)
  for(j in 1:M){
    Sol = matrix(0,1,n+1)
    Wt = matrix(0,1,n+1)
    for (i in 2:(n+1)){
      Wt[i] = Wt[i-1] + sqrt(delta)*rnorm(1,0,1)
    }
    Sol[1] = x0
    for(k in 2:(n+1))
    {
      DeltaW = Wt[k]-Wt[k-1]
      Sol[k] = Sol[k-1] +(-alpha*Sol[k-1]+((2*alpha*exp(-alpha*(TF-TiempoC[k-1])))/(1-exp(-2*alpha*(TF-TiempoC[k-1]))))*(xf-Sol[k-1]*exp(-alpha*(TF-TiempoC[k-1]))))*delta +sigma*DeltaW
    }
    Chau_OU[j,] <- Sol
  }
  return(Chau_OU)
}



