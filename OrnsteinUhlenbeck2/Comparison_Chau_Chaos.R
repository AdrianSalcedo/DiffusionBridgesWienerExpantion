require(stats)
#source('ChaosBridgeFunction.R')
source('Mi_main.R')
source('ExacBridgeFunction.r')
source('Milstein_codes.r')
####OU_bridge####
## Parametros del OU
thetaOU<-0.5
alpha <- thetaOU
sigmaOU<-1.0
sigma <- sigmaOU
delta=1/1000
### Parametros del puente:
a<- 0
b<- 0 # test aprobado 0.4
##numero de brownianos
nb=1000
#numero de puentes
M=1000
#numero de puntos por puente
n=1/delta
T0 = 0
TF=delta*n
TiempoC<-seq(T0, TF, length.out = n+1)
TL<-length(TiempoC)
#########
obs=501
p1 =0
p2=0
while (p1<0.05){
chaos<- BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,nb,M)
###################
Chau_OU <- Milstein_OU(M,a,b,alpha,sigma,n,delta,TiempoC)
p1=ks.test(Chau_OU[,obs],chaos[,obs])$p.value
print(p1)
par(mfrow=c(1,1))
qqplot(Chau_OU[,obs],chaos[,obs],ylab="Chau Method",xlab="Wiener Chaos approximation bridge")
abline(0,1)
}
 # pathw1 = 'D:\\DiffusionBridgesWienerExpantion\\DiffusionBridgesWienerExpantion\\OrnsteinUhlenbeck\\DATA_QQPlot\\Chau_Chaos'
 # Data<- rbind(TiempoC,chaos)
 # write.csv(Data,file.path(pathw1,"Chaos_OU_a0b2alpha0.5sigma1.0_BM100.csv"))
 # Data2<- rbind(TiempoC,Chau_OU)
 # write.csv(Data2,file.path(pathw1,"Chau_OU_a0b2alpha0.5sigma1.0_BM100.csv"))
chaos<- BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,nb,M)
while (p2<0.05){
  chaos<- BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,nb,M)
  Vikingos_OU = M_bridges_MM_OU(M,a,b,delta,n,thetaOU,sigmaOU) 
  p2=ks.test(Vikingos_OU[,obs],chaos[,obs])$p.value
  print(p2)
  par(mfrow=c(1,2))
  qqplot(Vikingos_OU[,obs],chaos[,obs],ylab="Vikingos Method",xlab="Wiener Chaos approximation bridge")
  abline(0,1)
  plot(TiempoC,Vikingos_OU[M/2,],type = 'l',col="red")
  lines(TiempoC,chaos[M/2,],type = 'l',col="blue")
}
pathw2 = 'D:\\DiffusionBridgesWienerExpantion\\DiffusionBridgesWienerExpantion\\OrnsteinUhlenbeck\\DATA_QQPlot\\Viki_Chaos'
Data<- rbind(TiempoC,chaos)
write.csv(Data,file.path(pathw2,"Chaos_OU_a0b2alpha0.5sigma1.0_BM1000.csv"))
Data3<- rbind(TiempoC,Vikingos_OU)
write.csv(Data3,file.path(pathw2,"Vikingos_OU_a0b2alpha0.5sigma1.0_BM1000.csv"))
# 
