Bri_Chaos_GBM=function(Xms,a,b,alpha,sigma,n,delta,nb,M)
{
  #### propagator for |m|>0
  for(i in 1:(Lp-1)){
    Yms[i+1,]=c((TiempoC[TL]-TiempoC[1:(TL-1)])*Integrate_Xms(1/(TF-TiempoC[-TL]),Xms[1+i,]),0)
  }
  #### propagator for |m|=0
  int1=c((TiempoC[TL]-TiempoC[1:(TL-1)])*Integrate_Xms(1/(TF-TiempoC[-TL]),Xms[1,]),0)
  for(i in 1:TL)
  {
    Yms[1,i]= a+(b-a)*TiempoC[i]/TF+int1[i]
  }
  return(Yms)
}
Gen_Bridges=function(Xms,a,b,alpha,sigma,n,delta,nb,M,BM){
  Yms <- Bri_Chaos_GBM(Xms,a,b,alpha,sigma,n,delta,nb,M) ##Genera propagador del puente
  for (k in 1:M){
    BM<-rnorm(nb,0,1)
    print(k)
    Xis=c(fXis_1(nb,TiempoC,BM),fXis_2_j(c(1,2),TiempoC,BM),fXis_2_j(c(1,3),TiempoC,BM),fXis_2_j(c(2,3),TiempoC,BM),fXis_2(nb,TiempoC,BM),fXis_3_12(c(1,2),TiempoC,BM),fXis_3_21(c(1,2),TiempoC,BM),fXis_3(nb,TiempoC,BM),fXis_4(nb,TiempoC,BM),fXis_5(nb,TiempoC,BM),fXis_6(nb,TiempoC,BM),fXis_7(nb,TiempoC,BM),fXis_8(nb,TiempoC,BM))
    M_Y[k,]=apply(Yms[2:Lp,]*Xis,2,sum)+Yms[1,]
  }
  Chaos_bridges = M_Y #borrar solo de pruebas
  return(M_Y)
}
#integral con respecto a Xms
Integrate_Xms=function(path,Xms)
{
  npoints=length(path)
  integral=numeric(npoints)
  for( i in 2:npoints)
  {
    integral[i]=path[i-1]*(Xms[i]-Xms[i-1])#*(path[i-1]+path[i])/2
  }
  integral=cumsum(integral)
  
  
  return(integral)
  
}
