####### My main
source('HermiteFunction.R')
source('BridgeChaosFunction.R')
source('IntegralItoTrapecio.R')
source('Vikngos_Geometric.R')
library(ggplot2)
#####
sigma=1
alpha=(sigma^2)/2
delta=1/1000
n=1000
nb=1000
a=exp(1)
b=exp(2)
# number of bridges
M=1000
#Xms=t(as.matrix(read.csv(file = 'C:/Users/DELL/Desktop/Propagador/Propagator_order12_100MB.csv'))[,-1])
Xms=t(as.matrix(read.csv(file = 'D:/DiffusionBridgesWienerExpantion/DiffusionBridgesWienerExpantion/GeometricBrownian/Propagator_order15_exp1_alphasigmasigma_sigma1_MB1000.csv'))[,-1])
Lp = dim(Xms)[[1]]

TF=delta*n
TiempoC<-seq(0, TF, length.out = n+1)
TL<-length(TiempoC)
M_Y=mat.or.vec(M,TL)
Yms=mat.or.vec(Lp,TL)
Wt <- matrix(0,nrow=M,nb+1)
Chaos_bridges <- Gen_Bridges(Xms,a,b,alpha,sigma,n,delta,nb,M,Wt)
Data<- rbind(TiempoC,Chaos_bridges)
write.csv(Data,"GBM_exp1_BM1000.csv")
MM_bridges<- M_bridges_MM_GBM(M,a,b,delta,n,alpha,sigma)
#########
obs=500
qqplot(Chaos_bridges[,obs],GBM,ylab="Exact bridge",xlab="Wiener chaos approximation bridge", cex.lab=0.5,main="Quantile-Quantile comparison", cex.main=0.5)
abline(0,1)
ks.test(Chaos_bridges[,obs],GBM)
###########

matplot(TiempoC,t(Chaos_bridges),type = "l")
matplot(TiempoC,t(MM_bridges),type = "l")
#BGBM = Bri_Chaos_GBM(Xms,a,b,alpha,sigma,n,delta,nb,M)
#BGBM=BridgesGBM(Xms,a,b,alpha,sigma,n,delta,nb,M)
#BGBM = Bri_Chaos_GBM(Xms,a,b,alpha,sigma,n,delta,nb,M)

