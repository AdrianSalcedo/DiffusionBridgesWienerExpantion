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