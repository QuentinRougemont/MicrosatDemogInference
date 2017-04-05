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
pop1 <-argv[2]
pop2 <-argv[3]
pop3 <-argv[4]

colon_count = paste(" awk -F ' ' '{print NF }' ", model , " |head -1 " )

#nlinesFul=as.numeric(strsplit(system("wc -l si.simul.ABC.txt", intern=T), " ")[[1]][1])
ncol = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])
#load simul
#m.si=matrix(scan(model), byrow=T, nrow=nlinesFul) 
m.si=matrix(scan(model), byrow=T, ncol=ncol) 

means <- colMeans(m.si, na.rm=TRUE)
for (j in 1:ncol(m.si)){
     m.si[is.na(m.si[, j]), j] <- means[j]
 }

#exclude parameters (contained in the first columns) of the simulation table to keep only summary statistics
param.si=c(2:8)
#target.stat=c(12:17,24:45)
target.stat=c(17:24,33:64,73:ncol(m.si))

obs=read.table("target_sumstats.txt", header=T)
#drop out Ho, NA,(NA=Number of allele not missing data!!):
obs2=obs[,-c(1:8,17:24,57:64)]

abc.posterior <- abc(target  = obs2, param = m.si[,param.si], sumstat = m.si[,target.stat], tol = 0.001, transf=c(rep("logit",length(param.si) )), 
	logit.bounds = rbind(range(m.si[, 2]), range(m.si[, 3]), range(m.si[, 4]), range(m.si[, 5]), range(m.si[,6]), range(m.si[, 7]), range(m.si[, 8]) ),  hcorr=T, method  = "neuralnet", numnet=20, sizenet=5)

model2 = strsplit(model,"ABC.txt")

write.table(abc.posterior$weights,
	paste("param.weight", model2, sep="."), quote=F,row.names=F,col.names=F)
write.table(abc.posterior$adj.values,
	paste("param.adj.value", model2, sep="."), quote=F,row.names=F)

write.table(summary(abc.posterior),
	paste("param.abc", model2,sep=".") ,
	col.names=c("Npop1","Npop2", "Npop3", "Npop2_3", "Nancestral", "Tsplit", "Tsplit2_3"),
	row.names=c("min","ic2.5","med","mean","mode","ic97.5", "max") ,quote=F)
#write.table(summary(abc.posteriori)[4,] 
#	,paste("param.abc.mean",model2, sep=".") ,col.names=F,row.names=F,quote=F)

#plot the graphic
pdf(file=paste("pior.posterior",model2, "pdf", sep="") ,width=12, height=8)
par(mfrow=c(2,4))
hist(m.si[,2],
     breaks=seq(0,max(m.si[,2]+0.1),0.5),
     col="gray94",
     freq=F,
     ylim=c(0,1),
     main=paste("theta", pop1, "/thetaRef", sep=" "),
     xlab=expression(theta[1]),
     ylab="probability density")
#hist((abc.posterior$unadj.values[,1] ), breaks=seq(0, max(m.si[,2]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.posterior$adj.values[,1]), breaks=seq(0,max(m.si[,2]+0.1),0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights)
box()
legend(
 legend = c("prior distribution","adjust posterior (neural net)" ), # "adjusted posterior(neural net)"),
col = c("gray94",rgb(1,0.3,0,0.6) ) , #rgb(1,0.3,0,0.6)),
pch = c(15,15), #15
x = "topright",
cex = 1,
bty ="n")
box()
hist(m.si[,3],
     breaks=seq(0,max(m.si[,3]+0.1),0.5),
     col="gray94",
     freq=F,
     ylim=c(0,1),
     main=paste("theta", pop2, "/thetaRef", sep=" "),
     xlab=expression(theta[2]),
     ylab="probability density")
#hist((abc.posterior$unadj.values[,2] ), breaks=seq(0,max(m.si[,3]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.posterior$adj.values[,2]), breaks=seq(0,max(m.si[,3]+0.1),0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights)
box()
hist(m.si[,4],
     breaks=seq(0,max(m.si[,4]+0.1),0.5),
     col="gray94",
     freq=F,
     ylim=c(0,0.2),
     main=paste("theta", pop3, "/thetaRef", sep=" "),
     xlab=expression(theta[3]),
     ylab="probability density")
#hist(abc.posterior$unadj.values[,3],   breaks=seq(0,max(m.si[,4]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,3], breaks=seq(0,max(m.si[,4]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
hist(m.si[,5],
	breaks=seq(0,max(m.si[,5]+0.1),0.5),
	col="gray94",
	freq=F,
	ylim=c(0,0.2),
	main=paste("theta", pop3, "+", pop2, "(4N gen)", sep=" "),
	xlab= expression(theta[4]),
	ylab="probab density")
#hist(abc.posterior$unadj.values[,4],   breaks=seq(0,max(m.si[,5]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,4], breaks=seq(0,max(m.si[,5]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()

hist(m.si[,6],
	breaks=seq(0,max(m.si[,6]+0.1),0.5),
	col="gray94",
	freq=F,
	ylim=c(0,0.2),
	main="theta Ancestral(4N gen)",
	xlab= expression(theta[ancestral]),
	ylab="probab density")
#hist(abc.posterior$unadj.values[,5],   breaks=seq(0,max(m.si[,5]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,5], breaks=seq(0,max(m.si[,5]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
hist(m.si[,7],
	breaks=seq(0,max(m.si[,7]+0.1),0.5),
	col="gray94",
	freq=F,
	ylim=c(0,0.5),
	main="divergence time (4N gen)",
	xlab= expression(tau),
	ylab="probab density")
#hist(abc.posterior$unadj.values[,6],   breaks=seq(0,max(m.si[,7]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,6], breaks=seq(0,max(m.si[,7]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
hist(m.si[,8],
	breaks=seq(0,max(m.si[,8]+0.1),0.5),
	col="gray94",
	freq=F,
	ylim=c(0,0.5),
	main=paste("divergence time",pop3,"-",pop2,"(4N gen)", sep=" "),
	xlab= expression(tau),
	ylab="probab density")
#hist(abc.posterior$unadj.values[,7],   breaks=seq(0,max(m.si[,8]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,7], breaks=seq(0,max(m.si[,8]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
dev.off()
