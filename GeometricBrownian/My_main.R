####### My main
source('HermiteFunction.R')
source('BridgeChaosFunction.R')
source('IntegralItoTrapecio.R')
source('Vikngos_Geometric.R')
library(ggplot2)
#####
alpha=0.2
sigma=0.3
delta=1/1000
n=1000
nb=100
a=0.2
b=0.3
# number of bridges
M=10
#Xms=t(as.matrix(read.csv(file = 'C:/Users/DELL/Desktop/Propagador/Propagator_order12_1000MB.csv'))[,-1])
Xms=t(as.matrix(read.csv(file = 'C:/Users/Usuario1/Desktop/Propagador/Propagator_order12_dt0.001MB100.csv'))[,-1])
Lp = dim(Xms)[[1]]

TF=delta*n
TiempoC<-seq(0, TF, length.out = n+1)
TL<-length(TiempoC)
M_Y=mat.or.vec(M,TL)
Yms=mat.or.vec(Lp,TL)
Wt <- matrix(rnorm(nb,0,1),nrow=M)
Chaos_bridges <- Gen_Bridges(Xms,a,b,alpha,sigma,n,delta,nb,M,Wt)
Data<- rbind(TiempoC,Chaos_bridges)
write.csv(Data,"GBM_step0.001_BM100.csv")
MM_bridges<- M_bridges_MM_GBM(M,a,b,delta,n,alpha,sigma)
#########
obs=500
qqplot(Chaos_bridges[,obs],MM_bridges[,obs],ylab="Exact bridge",xlab="Wiener chaos approximation bridge", cex.lab=0.5,main="Quantile-Quantile comparison", cex.main=0.5)
abline(0,1)
ks.test(Chaos_bridges[,obs],MM_bridges[,obs])
###########

matplot(TiempoC,t(Chaos_bridges),type = "l")
matplot(TiempoC,t(MM_bridges),type = "l")
#BGBM = Bri_Chaos_GBM(Xms,a,b,alpha,sigma,n,delta,nb,M)
#BGBM=BridgesGBM(Xms,a,b,alpha,sigma,n,delta,nb,M)
#BGBM = Bri_Chaos_GBM(Xms,a,b,alpha,sigma,n,delta,nb,M)

