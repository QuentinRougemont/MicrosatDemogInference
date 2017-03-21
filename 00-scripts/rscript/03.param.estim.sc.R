#!/usr/bin/Rscript 

library(abc)
library('weights')
require(laeken)
library('grid')
require(hexbin)

argv<-commandArgs(TRUE)

model<-argv[1]

#iso
nlinesFul=as.numeric(strsplit(system("wc -l im.simul.ABC.txt", intern=T), " ")[[1]][1])
#load simul
model="sc.simul.ABC.txt"
M_SC=matrix(scan(model), byrow=T, nrow=nlinesFul) 

means <- colMeans(M_SC, na.rm=TRUE)
for (j in 1:ncol(M_SC)){
     M_SC[is.na(M_SC[, j]), j] <- means[j]
 }

#exclude parameters (contained in the first columns) of the simulation table to keep only summary statistics
#delete Ho, et NA (NA=Number of allele not missing data!!)
param.sc=c(2:8)
target.stat=c(14:19,26:47)

obs=read.table("target_sumstats.txt", header=T)
#drop out Ho, NA,(NA=Number of allele not missing data!!):
obs2=obs[,-c(1:6,13:18)]

abc.posterior <- abc(target  = obs2, param = M_SC[,param.sc], sumstat = M_SC[,target.stat], tol = 0.001, transf=c("logit", "logit", "logit", "logit","logit", "logit","logit" ), 
	logit.bounds = rbind(range(M_SC[, 2]), range(M_SC[, 3]), range(M_SC[, 4]), range(M_SC[, 5]),range(M_SC[, 6]),range(M_SC[, 7]), range(M_SC[,8])),  hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)

write.table(abc.posterior$weights,"weight.sc",quote=F,row.names=F,col.names=F)
write.table(abc.posterior$adj.values,"adjval.sc",quote=F,row.names=F)

z<-summary(abc.posterior)

write.table(z,"abc.sc",col.names=T,row.names=F,quote=F)
write.table(z[4,],"abc.mean.sc",col.names=F,row.names=F,quote=F)

#plot graphe
pdf(file="pior.posterior.sc.pdf",width=12,height=8)
par(mfrow=c(2,4))

hist(M_SC[,2],
     breaks=seq(0,max(M_SC[,2]+0.1),0.5),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="theta1/thetaRef",
     xlab=expression((theta[1])),
     ylab="probability density")
hist((abc.posterior$unadj.values[,1] ),  breaks=seq(0,max(M_SC[,2]+0.1),0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.posterior$adj.values[,1]), breaks=seq(0,max(M_SC[,2]+0.1),0.1), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights)
box()

legend(
 legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
pch = c(15,15,15),
x = "topleft",
cex = 1,
bty ="n")
box()

hist(M_SC[,3],
     breaks=seq(0,max(M_SC[,3]+0.1),0.5),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="theta2/thetaRef",
     xlab=expression((theta[2])),
     ylab="probability density")
hist((abc.posterior$unadj.values[,2] ),  breaks=seq(0,max(M_SC[,3]+0.1),0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.posterior$adj.values[,2]), breaks=seq(0,max(M_SC[,3]+0.1),0.1), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights)

box()

hist(M_SC[,4],
     breaks=seq(0,max(M_SC[,4]+0.1),0.05),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="thetaAncien/thetaRef",
     xlab=expression((theta[Ancien])),
     ylab="probability density")
hist(abc.posterior$unadj.values[,3],   breaks=seq(0,max(M_SC[,4]+0.1),0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,3], breaks=seq(0,max(M_SC[,4]+0.1),0.1), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()

hist(M_SC[,5],
	breaks=seq(0,max(M_SC[,5]+0.1),1),
	col="lightgrey",
	freq=F,
	ylim=c(0,2),
	main="Effective migration rate ",
	xlab= "M1<-2",
	ylab="probability density")
hist(abc.posterior$unadj.values[,4],   breaks=seq(0,max(M_SC[,5]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,4], breaks=seq(0,max(M_SC[,5]+0.1),0.5), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()
#
hist(M_SC[,6],
	breaks=seq(0,max(M_SC[,6]+0.1),1),
	col="lightgrey",
	freq=F,
	ylim=c(0,0.5),
	main="Effective migration rate ",
	xlab= "M2<-1",
	ylab="probab density")
hist(abc.posterior$unadj.values[,5],   breaks=seq(0,max(M_SC[,5]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,5], breaks=seq(0,max(M_SC[,5]+0.1),0.5), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()
                  
hist(M_SC[,7],
	breaks=seq(0,max(M_SC[,7]+0.1),1),
	col="lightgrey",
	freq=F,
	ylim=c(0,0.5),
	main="time of secondary contact( 4N generations)",
	xlab= expression(tau),
	ylab="probab density")
hist(abc.posterior$unadj.values[,6],   breaks=seq(0,max(M_SC[,7]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,6], breaks=seq(0,max(M_SC[,7]+0.1),0.5), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
hist(M_SC[,8],
        breaks=seq(0,max(M_SC[,7]+0.1),1),
        col="lightgrey",
        freq=F,
        ylim=c(0,0.5),
        main="divergence time( 4N generations)",
        xlab= expression(tau),
        ylab="probab density")
hist(abc.posterior$unadj.values[,6],   breaks=seq(0,max(M_SC[,7]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,6], breaks=seq(0,max(M_SC[,7]+0.1),0.5), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()

dev.off()

#sans GW et R
target.stat=c(14:19,26:31,44:47)

obs=read.table("target_sumstats.txt", header=T)
#drop out Ho, NA,(NA=Number of allele not missing data!!):
obs2=obs2[,-c(13:24)]


abc.posterior <- abc(target  = obs2, param = M_SC[,param.sc], sumstat = M_SC[,target.stat], tol = 0.001, transf=c("logit", "logit", "logit", "logit","logit", "logit" ,"logit"), 
	logit.bounds = rbind(range(M_SC[, 2]), range(M_SC[, 3]), range(M_SC[, 4]), range(M_SC[, 5]),range(M_SC[, 6]),range(M_SC[, 7]), range(M_SC[,8])),  hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)

write.table(abc.posterior$weights,"weight.sc_v2",quote=F,row.names=F,col.names=F)
write.table(abc.posterior$adj.values,"adjval.sc_v2",quote=F,row.names=F)

z<-summary(abc.posterior)

write.table(z,"abc.sc_v2",col.names=T,row.names=F,quote=F)
write.table(z[4,],"abc.mean.sc_v2",col.names=F,row.names=F,quote=F)

pdf(file="pior.posterior.sc_v2.pdf",width=12,height=8)
par(mfrow=c(2,3))

hist(M_SC[,2],
     breaks=seq(0,max(M_SC[,2]+0.1),0.1),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="theta1/thetaRef",
     xlab=expression((theta[1])),
     ylab="probability density")
hist((abc.posterior$unadj.values[,1] ),  breaks=seq(0,max(M_SC[,2]+0.1),0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.posterior$adj.values[,1]), breaks=seq(0,max(M_SC[,2]+0.1),0.1), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights)
box()

legend(
 legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
pch = c(15,15,15),
x = "topleft",
cex = 1,
bty ="n")
box()

hist(M_SC[,3],
     breaks=seq(0,max(M_SC[,3]+0.1),0.1),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="theta2/thetaRef",
     xlab=expression((theta[2])),
     ylab="probability density")
hist((abc.posterior$unadj.values[,2] ),  breaks=seq(0,max(M_SC[,3]+0.1),0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.posterior$adj.values[,2]), breaks=seq(0,max(M_SC[,3]+0.1),0.1), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights)

box()

hist(M_SC[,4],
     breaks=seq(0,max(M_SC[,4]+0.1),0.1),
     col="grey",
     freq=F,
     ylim=c(0,1),
     main="thetaAncien/thetaRef",
     xlab=expression((theta[Ancien])),
     ylab="probability density")
hist(abc.posterior$unadj.values[,3],   breaks=seq(0,max(M_SC[,4]+0.1),0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,3], breaks=seq(0,max(M_SC[,4]+0.1),0.1), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()

hist(M_SC[,5],
	breaks=seq(0,max(M_SC[,5]+0.1),1),
	col="lightgrey",
	freq=F,
	ylim=c(0,2),
	main="Effective migration rate ",
	xlab= "M1<-2",
	ylab="probability density")
hist(abc.posterior$unadj.values[,4],   breaks=seq(0,max(M_SC[,5]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,4], breaks=seq(0,max(M_SC[,5]+0.1),0.5), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()
#
hist(M_SC[,6],
	breaks=seq(0,max(M_SC[,6]+0.1),1),
	col="lightgrey",
	freq=F,
	ylim=c(0,0.5),
	main="Effective migration rate ",
	xlab= "M2<-1",
	ylab="probab density")
hist(abc.posterior$unadj.values[,5],   breaks=seq(0,max(M_SC[,5]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,5], breaks=seq(0,max(M_SC[,5]+0.1),0.5), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()

hist(M_SC[,7],
        breaks=seq(0,max(M_SC[,7]+0.1),1),
        col="lightgrey",
        freq=F,
        ylim=c(0,0.5),
        main="time of secondary contact( 4N generations)",
        xlab= expression(tau),
        ylab="probab density")
hist(abc.posterior$unadj.values[,6],   breaks=seq(0,max(M_SC[,7]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,6], breaks=seq(0,max(M_SC[,7]+0.1),0.5), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
hist(M_SC[,8],
        breaks=seq(0,max(M_SC[,7]+0.1),1),
        col="lightgrey",
        freq=F,
        ylim=c(0,0.5),
        main="divergence time( 4N generations)",
        xlab= expression(tau),
        ylab="probab density")
hist(abc.posterior$unadj.values[,6],   breaks=seq(0,max(M_SC[,7]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,6], breaks=seq(0,max(M_SC[,7]+0.1),0.5), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()

dev.off()
