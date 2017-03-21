#!/usr/bin/Rscript 

rm(list=ls())
ls()

library(abc)
library('weights')
require(laeken)
library('grid')
require(hexbin)

argv<-commandArgs(TRUE)

model<-argv[1]

#iso
nlinesFul=as.numeric(strsplit(system("wc -l si.simul.ABC.txt", intern=T), " ")[[1]][1])
#load simul
M_SI=matrix(scan(model), byrow=T, nrow=nlinesFul) 

means <- colMeans(M_SI, na.rm=TRUE)
for (j in 1:ncol(M_SI)){
     M_SI[is.na(M_SI[, j]), j] <- means[j]
 }

#exclude parameters (contained in the first columns) of the simulation table to keep only summary statistics
param.si=c(2:5)
target.stat=c(12:17,24:45)

obs=read.table("target_sumstats.txt", header=T)
#drop out Ho, NA,(NA=Number of allele not missing data!!):
obs2=obs[,-c(1:6,13:18)]

abc.posterior <- abc(target  = obs2, param = M_SI[,param.si], sumstat = M_SI[,target.stat], tol = 0.001, transf=c("logit", "logit", "logit", "logit"), 
	logit.bounds = rbind(range(M_SI[, 2]), range(M_SI[, 3]), range(M_SI[, 4]), range(M_SI[, 5])),  hcorr=T, method  = "neuralnet", numnet=20, sizenet=5)


write.table(abc.posterior$weights,"weight.si",quote=F,row.names=F,col.names=F)
write.table(abc.posterior$adj.values,"adjval.si",quote=F,row.names=F)

z<-summary(abc.posterior)

write.table(z,"abc.si",col.names=T,row.names=F,quote=F)
write.table(z[4,],"abc.mean.si",col.names=F,row.names=F,quote=F)

#plot the graphic

pdf(file="pior.posterior.si.pdf",width=12,height=8)
par(mfrow=c(2,3))

hist(M_SI[,2],
     breaks=seq(0,max(M_SI[,2]+0.1),0.5),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="theta1/thetaRef",
     xlab=expression((theta[1])),
     ylab="probability density")
hist((abc.posterior$unadj.values[,1] ), breaks=seq(0, max(M_SI[,2]+0.1),0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.posterior$adj.values[,1]), breaks=seq(0,max(M_SI[,2]+0.1),0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights)
box()

legend(
 legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
pch = c(15,15,15),
x = "topleft",
cex = 1,
bty ="n")
box()

hist(M_SI[,3],
     breaks=seq(0,max(M_SI[,3]+0.1),0.5),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="theta2/thetaRef",
     xlab=expression((theta[2])),
     ylab="probability density")
hist((abc.posterior$unadj.values[,2] ), breaks=seq(0,max(M_SI[,3]+0.1),0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.posterior$adj.values[,2]), breaks=seq(0,max(M_SI[,3]+0.1),0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights)

box()

hist(M_SI[,4],
     breaks=seq(0,max(M_SI[,4]+0.1),0.05),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="thetaAncien/thetaRef",
     xlab=expression((theta[Ancien])),
     ylab="probability density")
hist(abc.posterior$unadj.values[,3],   breaks=seq(0,max(M_SI[,4]+0.1) ,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,3], breaks=seq(0,max(M_SI[,4]+0.1) ,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()

hist(M_SI[,5],
	breaks=seq(0,max(M_SI[,5]+0.1),0.5),
	col="lightgrey",
	freq=F,
	ylim=c(0,0.5),
	main="divergence time( 4N generations)",
	xlab= expression(tau),
	ylab="probab density")
hist(abc.posterior$unadj.values[,4],   breaks=seq(0,max(M_SI[,5]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,4], breaks=seq(0,max(M_SI[,5]+0.1) ,0.5), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()

dev.off()

########################################################################
#Sans R et GW:
obs2=obs2[,-c(13:24)]
target.stat=c(12:17,24:29,42:45)

abc.posterior <- abc(target  = obs2, param = M_SI[,param.si], sumstat = M_SI[,target.stat], tol = 0.001, transf=c("logit", "logit", "logit", "logit"), 
	logit.bounds = rbind(range(M_SI[, 2]), range(M_SI[, 3]), range(M_SI[, 4]), range(M_SI[, 5])),  hcorr=T, method  = "neuralnet", numnet=20, sizenet=5)


write.table(abc.posterior$weights,"weight.si_v2",quote=F,row.names=F,col.names=F)
write.table(abc.posterior$adj.values,"adjval.si_v2",quote=F,row.names=F)

z<-summary(abc.posterior)

write.table(z,"abc.si_v2",col.names=T,row.names=F,quote=F)
write.table(z[4,],"abc.mean.si_v2",col.names=F,row.names=F,quote=F)

pdf(file="pior.posterior.si_v2.pdf",width=12,height=8)
par(mfrow=c(2,3))

hist(M_SI[,2],
     breaks=seq(0,max(M_SI[,2]+0.1),0.5),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="theta1/thetaRef",
     xlab=expression((theta[1])),
     ylab="probability density")
hist((abc.posterior$unadj.values[,1] ), breaks=seq(0, max(M_SI[,2]+0.1),0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.posterior$adj.values[,1]), breaks=seq(0,max(M_SI[,2]+0.1),0.1), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights)
box()

legend(
 legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
pch = c(15,15,15),
x = "topleft",
cex = 1,
bty ="n")
box()

hist(M_SI[,3],
     breaks=seq(0,max(M_SI[,3]+0.1),0.5),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="theta2/thetaRef",
     xlab=expression((theta[2])),
     ylab="probability density")
hist((abc.posterior$unadj.values[,2] ), breaks=seq(0,max(M_SI[,3]+0.1),0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.posterior$adj.values[,2]), breaks=seq(0,max(M_SI[,3]+0.1),0.1), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights)

box()

hist(M_SI[,4],
     breaks=seq(0,max(M_SI[,4]+0.1),0.5),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="thetaAncien/thetaRef",
     xlab=expression((theta[Ancien])),
     ylab="probability density")
hist(abc.posterior$unadj.values[,3],   breaks=seq(0,max(M_SI[,4]+0.1) ,0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,3], breaks=seq(0,max(M_SI[,4]+0.1) ,0.1), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()

hist(M_SI[,5],
	breaks=seq(0,max(M_SI[,5]+0.1),1),
	col="lightgrey",
	freq=F,
	ylim=c(0,0.5),
	main="divergence time( 4N generations)",
	xlab= expression(tau),
	ylab="probab density")
hist(abc.posterior$unadj.values[,4],   breaks=seq(0,max(M_SI[,5]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,4], breaks=seq(0,max(M_SI[,5]+0.1) ,0.5), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()

dev.off()

