#p_test <- 0
#while(p_test<0.05){
require(stats)
#source('ChaosBridgeFunction.R')
source('Mi_main.R')
source('ExacBridgeFunction.r')
source('Milstein_codes.r')
####OU_bridge####
## Parametros del OU
thetaOU<-0.5
sigmaOU<-1.0
delta=1/1000
### Parametros del puente:
a<- 0
b<- 0 # test aprobado 0.4
##numero de brownianos
nb=100
#numero de puentes
M=1000
#numero de puntos por puente
n=1/delta
T0 = 0
TF=delta*n
TiempoC<-seq(T0, TF, length.out = n+1)
TL<-length(TiempoC)

#ejemplo
exa<- M_bridges_Exact_OU(M,a,b,delta,n,thetaOU,sigmaOU)
#start_time_v <- Sys.time()
#Vikingos<-M_bridges_MM_OU(M,a,b,delta,n,thetaOU,sigmaOU)
#end_time_v <- Sys.time()
#total_time_v <- end_time_v - start_time_v
#total_time_v
#system.time({M_bridges_MM_OU(M,a,b,delta,n,thetaOU,sigmaOU)})
#########
start_time_c <- Sys.time()
chaos<- BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,nb,M)
end_time_c <- Sys.time()
total_time_c <- end_time_c - start_time_c
total_time_c
#system.time({BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,nb,M)})
#Data<- rbind(TiempoC,chaos)
#write.csv(Data,"UBB_B0.80.5_mb0.csv")
obs=n/2
#par(mfrow=c(1,2))
#########Vikingos comparation ####################
#qqplot(exa[,obs],Vikingos[,obs])
#abline(0,1)
qqplot(exa[,obs],chaos[,obs],ylab="Exact bridge",xlab="Wiener chaos approximation bridge")
abline(0,1)
##################################################
#test_vikingos <-ks.test(Vikingos[,obs],exa[,obs])
#test_vikingos
####### acepp mayor a 0.05 ##########
test_chaos <-ks.test(exa[,obs],chaos[,obs])
test_chaos 
print(test_chaos["p.value"])
print(total_time_c)
p_test<-test_chaos[["p.value"]]
#}
###
par(mfrow=c(1,1))
plot(TiempoC,chaos[1,],type="l",main="Ornstein Uhlenbeck bridges",
      xlab="Time",cex.lab=0.5,ylab="Diffusion bridges",cex.main=0.5,
      ylim=c(-4,4))
lines(TiempoC,exa[1,],col="red")
write.csv(chaos,"ChaosWiner10000.csv")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


BChaos1=BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,1,M)
BChaos10=BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,10,M)
BChaos100=BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,100,M)
BChaos500=BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,500,M)
BChaos1000=BridgeChaos(a,b,thetaOU,sigmaOU,n,delta,1000,M)
BChaos1000_10=BridgeChaos(a,b,thetaOU,sigmaOU,n,1/10,nb,M)
BChaos1000_100=BridgeChaos(a,b,thetaOU,sigmaOU,n,1/100,nb,M)
BChaos1000_1000=BridgeChaos(a,b,thetaOU,sigmaOU,n,1/1000,nb,M)
dev.new(4,4)
plot( TiempoC,BChaos1[1,],type="l",main="Ornstein Uhlenbeck bridges",xlab="Time",cex.lab=0.5,ylab="Diffusion bridges",cex.main=0.5,ylim=c(-3,3))
lines(TiempoC,BChaos10[1,],col="red")
lines(TiempoC,BChaos100[1,],col="blue")
lines(TiempoC,BChaos500[1,],col="green")
lines(TiempoC,BChaos1000[1,],col="orange")
legend(x = "topright",  
       cex = 0.5,# Position
       legend = c("1", "10","100","500","1000"),  # Legend texts
       lty = c(1, 1,1,1,1),           # Line types
       col = c("black","red","blue","green","orange"),           # Line colors
       lwd = 1) 
lines(M_OUE1[3,],col="blue")
lines(M_OUE1[4,],col="green")

EB_OU=M_bridges_MM_OU(M,a,b,delta,n,thetaOU,sigmaOU)


obs=500

qqplot(BChaos1000[,obs],EB_OU,ylab="Exact bridge",xlab="Wiener chaos approximation bridge", cex.lab=0.5,main="Quantile-Quantile comparison", cex.main=0.5)
abline(0,1)
qqplot(BChaos1000_100[,obs],M_OUE100[,obs],ylab="Exact bridge",xlab="Wiener chaos approximation bridge", cex.lab=0.5,main="Quantile-Quantile comparison", cex.main=0.5)
abline(0,1)
qqplot(BChaos1000_10[,obs],M_OUE10[,obs],ylab="Exact bridge",xlab="Wiener chaos approximation bridge", cex.lab=0.5,main="Quantile-Quantile comparison", cex.main=0.5)
abline(0,1)


plot(density(BChaos[,obs]))
lines(density(M_OUE1[,obs]),col="red")