library(sde)
library(plotly)
library(ggplot2)
library(MASS)
library(deSolve)
library(stats4)
library(fda)
library(splines)
library(fds)
library(Matrix)
library(rainbow)
library(pcaPP)
library(RCurl)
library(foreach)

simu1<-as.matrix(read.csv("Propagator_1000MB.csv",header = T))
#simuSol<-read.csv("Solution_plot.csv",header = T)
#Tiempo<-read.csv("Time_Brownian_motion2023.csv",header = T)
TiempoC<-simu1[,1]
simu<-Propagador <- simu1[,2:1018]

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
  
 ### Parametros del puente :
theta<-0.3
eta<-0.2

T1<-1
t0<-0


M<-2
TL<-length(TiempoC)

M_Y<-M_X<-mat.or.vec(TL,M)
Vector_xi<-c()
#M_Y<-mat.or.vec(TL,M)
Y<-mat.or.vec(TL,1)
#W <- rep(0, 100)
vec_i = c()
for (i in 1:1000) {
  vec_i[i] = i
  
}

start_time <- Sys.time()
for(k in 1:M)
{
  
  BM<-My_SBM(1000,TiempoC) 
  #W<-BM-BM
  #M_X[1,k]<-0.1 #### initial condition
  W <- rep(0, 1000)
    
  for(t in 2:TL)
  {
      z <-BM[1:1000,t]
      zr <-BM[1:1000,t-1]

      W <- W + sqrt(2)*sin(vec_i*pi*TiempoC[t-1])*(z-zr)
      Vector_xi[1]<-1
      Vector_xi[2:1001] <- W
      Vector_xi[1002]<- (W[1]) * (W[2])
      Vector_xi[1003]<- (W[1]) * (W[3])
      Vector_xi[1004]<- (W[2]) * (W[3])
      Vector_xi[1005]<- sqrt(2)*((W[1]) ^ 2 - 1) / sqrt(2)
      Vector_xi[1006]<- sqrt(2)*((W[2]) ^ 2 - 1) / sqrt(2)
      Vector_xi[1007]<- sqrt(2)*((W[3]) ^ 2 - 1) / sqrt(2)
      Vector_xi[1008]<- sqrt(2)*(W[1]) * ((W[2]) ^ 2 - 1) / sqrt(2)
      Vector_xi[1009]<- sqrt(2)*(((W[1]) ^ 2 - 1)  * (W[2]) / sqrt(2))
      Vector_xi[1010]<- sqrt(6)*(W[1]) * (((W[1]) ^ 2 - 3) / sqrt(6))
      Vector_xi[1011]<- sqrt(6)*(W[2]) * (((W[2]) ^ 2 - 3) / sqrt(6))
      Vector_xi[1012]<-sqrt(24)* (W[1] ^ 4 - 6 *W[1] ^ 2 + 3) / (2 * sqrt(6))
      Vector_xi[1013]<-sqrt(24)* (W[2] ^ 4 - 6 *W[2] ^ 2 + 3) / (2 * sqrt(6))
      Vector_xi[1014]<- sqrt(120)*(W[1]*(15 - 10*W[1]^2 + W[1]^4))/(2*sqrt(30))
      Vector_xi[1015]<- sqrt(120)*(W[2]*(15 - 10*W[2]^2 + W[2]^4))/(2*sqrt(30))
      Vector_xi[1016]<- sqrt(720)*(-15 + 45*W[1]^2 - 15*W[1]^4 + W[1]^6)/(12*sqrt(5))
      Vector_xi[1017]<- sqrt(720)*(-15 + 45*W[2]^2 - 15*W[2]^4 + W[2]^6)/(12*sqrt(5))
#############  HERMITE polynomials
   # H1 = x
   # H2 = (x² - 1) / sqrt(2)
   # H3 = x (x8² - 3) / sqrt(6)
   # H4 (x⁴ - 6 x² + 3) / (2 sqrt(6))
   # H5 = x (x⁴ - 10 x² + 15) / 2 sqrt(30))
   # H6 = (-15 + 45 x^2 - 15 x^4 + x^6)/(12 sqrt(5))
   ##################################
 suma<-Ym<-simu-simu

  for(h in 2:1017) #### for |m|>0
  {
    suma[1,h]<-Ym[1,h]<-simu[1,h]

    for(i in 2:TL)
    {
       suma[i,h]<- suma[i-1,h] + (simu[i,h]-simu[i-1,h])/(T1-TiempoC[i-1])
      Ym[i,h]<- (T1-TiempoC[i])*suma[i,h]
    }
  }

   #### for |m|=0
   suma[1,1]<-simu[1,1]- simu[1,1]

   Ym[1,1]<-eta

    for(i in 2:TL)
    {
     suma[i,1]<- suma[i-1,1] + (simu[i,1]-simu[i-1,1])/(T1-TiempoC[i-1])
     Ym[i,1]<- eta + (theta-eta ) *(TiempoC[i]/T1) + (T1-TiempoC[i])*suma[i-1,1]

    }
   M_Y[t,k] <- t(Ym[t,])%*%Vector_xi
  } ## end for t
  M_Y[1,k] <- t(Ym[1,])%*%Vector_xi
  print(k)
} #### end for k
end_time <- Sys.time()
end_time - start_time

Yz<-rowMeans(M_Y)
plot( M_X, type = "l")#, lty = 1:5)
plot(TiempoC, Yz, type = "l", lty = 1:5)
matplot(TiempoC,M_Y, type = "l", lwd=1)
df <- data.frame(M_Y)
write.csv(df, "Datos2_theta_0.3_eta_0.2_2023.csv")
