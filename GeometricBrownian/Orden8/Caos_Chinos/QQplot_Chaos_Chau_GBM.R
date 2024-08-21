Chaos_OU_00 = read.csv(file ='Chaos_GBM_prob8_a0.2b0.3_BM1000.csv')
Chau_OU_00 = read.csv(file ='ChauM_GBM_a0.2b0.3.csv')
Chaos_OU_01 = read.csv(file ='Chaos_GBM_prob8_a1b1_alpha0.2_sigma0.3_BM1000.csv')
Chau_OU_01 = read.csv(file ='ChauM_GBM_prob8_a1b1_alpha0.2_sigma0.3_BM1000.csv')
Chaos_OU_02 = read.csv(file ='Chaos_GBM_prob8_a1b0.1_alpha0.2_sigma0.3_BM1000.csv')
Chau_OU_02 = read.csv(file ='ChauM_GBM_prob8_a1b0.1_alpha0.2_sigma0.3_BM1000.csv')
Chaos_OU_0.80.5 = read.csv(file ='Chaos_GBM_prob8_a3b5_alpha0.2_sigma0.3_BM1000_1.csv')
Chau_OU_0.80.5 = read.csv(file ='ChauM_GBM_prob8_a3b5_alpha0.2_sigma0.3_BM1000_1.csv')

Tiempo = Chaos_OU_00[1,]
Chaos_OU_00 = Chaos_OU_00[2:1001,]
Chau_OU_00 = Chau_OU_00[2:1001,]
Chaos_OU_01 = Chaos_OU_01[2:501,]
Chau_OU_01 = Chau_OU_01[2:501,]
Chaos_OU_02 = Chaos_OU_02[2:1001,]
Chau_OU_02 = Chau_OU_02[2:1001,]
Chaos_OU_0.80.5 = Chaos_OU_0.80.5[2:501,]
Chau_OU_0.80.5 = Chau_OU_0.80.5[2:501,]
obs=501
p1=ks.test(Chaos_OU_00[,obs],Chau_OU_00[,obs])$p.value
p2=ks.test(Chaos_OU_01[,obs],Chau_OU_01[,obs])$p.value
p3=ks.test(Chaos_OU_02[,obs],Chau_OU_02[,obs])$p.value
p4=ks.test(Chaos_OU_0.80.5[,obs],Chau_OU_0.80.5[,obs])$p.value
print(c(p1,p2,p3,p4))
par(mfrow=c(1,1))
qqplot(Chau_OU_00[,obs],Chaos_OU_00[,obs],ylab="Lyons and Zheng's approach",xlab="Wiener Chaos approach")
abline(0,1)
qqplot(Chau_OU_01[,obs],Chaos_OU_01[,obs],ylab="Chau method",xlab="Wiener Chaos approximation method")
abline(0,1)
qqplot(Chau_OU_02[,obs],Chaos_OU_02[,obs],ylab="Chau method",xlab="Wiener Chaos approximation method")
abline(0,1)
qqplot(Chau_OU_0.80.5[,obs],Chaos_OU_0.80.5[,obs],ylab="Chau method",xlab="Wiener Chaos approximation method")
abline(0,1)
