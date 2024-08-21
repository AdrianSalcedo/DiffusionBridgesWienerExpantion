Chaos_OU_00 = read.csv(file ='DATA_QQPlotChaos_OU_a0b0alpha0.5sigma1.0_BM1000.csv')
Chau_OU_00 = read.csv(file ='DATA_QQPlotVikingos_OU_a0b0alpha0.5sigma1.0.csv')
Chaos_OU_01 = read.csv(file ='D:/DiffusionBridgesWienerExpantion/DiffusionBridgesWienerExpantion/OrnsteinUhlenbeck/Chaos_OU_a0b1.0alpha0.5sigma1.0_BM1000.csv')
Chau_OU_01 = read.csv(file ='D:/DiffusionBridgesWienerExpantion/DiffusionBridgesWienerExpantion/OrnsteinUhlenbeck/Chau_OU_a0b1.0alpha0.5sigma1.0.csv')
Chaos_OU_02 = read.csv(file ='D:/DiffusionBridgesWienerExpantion/DiffusionBridgesWienerExpantion/OrnsteinUhlenbeck/Chaos_OU_a0b2.0alpha0.5sigma1.0_BM1000.csv')
Chau_OU_02 = read.csv(file ='D:/DiffusionBridgesWienerExpantion/DiffusionBridgesWienerExpantion/OrnsteinUhlenbeck/Chau_OU_a0b2.0alpha0.5sigma1.0.csv')
Chaos_OU_0.80.5 = read.csv(file ='D:/DiffusionBridgesWienerExpantion/DiffusionBridgesWienerExpantion/OrnsteinUhlenbeck/Chaos_OU_a0.8b0.5alpha0.5sigma1.0_BM1000.csv')
Chau_OU_0.80.5 = read.csv(file ='D:/DiffusionBridgesWienerExpantion/DiffusionBridgesWienerExpantion/OrnsteinUhlenbeck/Chau_OU_a0.8b0.5alpha0.5sigma1.0.csv')

Tiempo = Chaos_OU_00[1,]
Chaos_OU_00 = Chaos_OU_00[2:1001,]
Chau_OU_00 = Chau_OU_00[2:1001,]
Chaos_OU_01 = Chaos_OU_01[2:1001,]
Chau_OU_01 = Chau_OU_01[2:1001,]
Chaos_OU_02 = Chaos_OU_02[2:1001,]
Chau_OU_02 = Chau_OU_02[2:1001,]
Chaos_OU_0.80.5 = Chaos_OU_0.80.5[2:1001,]
Chau_OU_0.80.5 = Chau_OU_0.80.5[2:1001,]
obs=501
p1=ks.test(Chaos_OU_00[,obs],Chau_OU_00[,obs])$p.value
p2=ks.test(Chaos_OU_01[,obs],Chau_OU_01[,obs])$p.value
p3=ks.test(Chaos_OU_02[,obs],Chau_OU_02[,obs])$p.value
p4=ks.test(Chaos_OU_0.80.5[,obs],Chau_OU_0.80.5[,obs])$p.value
print(c(p1,p2,p3,p4))
par(mfrow=c(1,1))
qqplot(Chau_OU_00[,obs],Chaos_OU_00[,obs],ylab="Sorensen and Bladt approach",xlab="Wiener Chaos approach")
abline(0,1)
qqplot(Chau_OU_01[,obs],Chaos_OU_01[,obs],ylab="Chau method",xlab="Wiener Chaos approximation method")
abline(0,1)
qqplot(Chau_OU_02[,obs],Chaos_OU_02[,obs],ylab="Chau method",xlab="Wiener Chaos approximation method")
abline(0,1)
qqplot(Chau_OU_0.80.5[,obs],Chaos_OU_0.80.5[,obs],ylab="Chau method",xlab="Wiener Chaos approximation method")
abline(0,1)
