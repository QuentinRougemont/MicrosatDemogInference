#!/usr/bin/Rscript 

rm(list=ls())
ls()

library(abc)
library('weights')
require(laeken)
library('grid')
require(hexbin)

argv<-commandArgs(TRUE)

model  <-argv[1]
header <-argv[2]
pop1   <-argv[3]#al_atl
pop2   <-argv[4]#af_atl
pop3   <-argv[5]#af_rhone
pop4   <-argv[6]#af_corse
colon_count = paste(" awk -F ' ' '{print NF }' ", model , " |head -1 " )

#nlinesFul=as.numeric(strsplit(system("wc -l si.simul.ABC.txt", intern=T), " ")[[1]][1])
ncol = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])
#load simul
#m.si=matrix(scan(model), byrow=T, nrow=nlinesFul) 
m.si=matrix(scan(model), byrow=T, ncol=ncol)
header=read.table(header)
means <- colMeans(m.si, na.rm=TRUE)
for (j in 1:ncol(m.si)){
         m.si[is.na(m.si[, j]), j] <- means[j]
 }

#exclude parameters (contained in the first columns) of the simulation table to keep only summary statistics
param.si=c(2:11)
#target.stat=c(12:17,24:45)
target.stat=c(22,24,26,28,30,42,44,46,48,50,62,64,66,68,70,72,74,76,78,80,92,94,96,98,100,102,104,106,108,110,112,114)

obs=read.table("target_sumstats.txt", header=T)
#drop out Ho, NA,(NA=Number of allele not missing data!!):
obs2=obs[,c(11,13,15,17,19,31,33,35,37,39,51,53,55,57,59,61,63,65,67,69,81,83,85,87,89,91,93,95,97,99,101,103)]

bound=NULL
log.bound=NULL

for(i in 2:length(param.si)+1)
{
            bound=(range(m.si[,i]))
        log.bound=rbind(log.bound,bound)
}

abc.posterior <- abc(target  = obs2, param = m.si[,param.si], sumstat = m.si[,target.stat], tol = 0.001, transf=c(rep("logit",length(param.si) )),
                             logit.bounds = log.bound,  hcorr=T, method  = "neuralnet", numnet=20, sizenet=5)

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
#       ,paste("param.abc.mean",model2, sep=".") ,col.names=F,row.names=F,quote=F)

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
     lab="probability density")
#hist(abc.posterior$unadj.values[,3],   breaks=seq(0,max(m.si[,4]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,3], breaks=seq(0,max(m.si[,4]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
hist(m.si[,5],
        breaks=seq(0,max(m.si[,5]+0.1),0.5),
        col="gray94",
        freq=F,
        ylim=c(0,0.2),
        main=paste("theta", pop4, "(4N gen)", sep=" "),
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
        main=paste("theta", pop4, "+" ,pop3,sep=" "),
        xlab= expression(theta[5]),
        ylab="probab density")
#hist(abc.posterior$unadj.values[,5],   breaks=seq(0,max(m.si[,5]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,5], breaks=seq(0,max(m.si[,5]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
hist(m.si[,7],
        breaks=seq(0,max(m.si[,7]+0.1),0.5),
        col="gray94",
        freq=F,
        ylim=c(0,0.5),
        main=paste("thetaAncestral1",sep=" "),
        xlab= expression(theta[Ancestral1]),
        ylab="probab density")
#hist(abc.posterior$unadj.values[,6],   breaks=seq(0,max(m.si[,7]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,6], breaks=seq(0,max(m.si[,7]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
hist(m.si[,8],
        breaks=seq(0,max(m.si[,8]+0.1),0.5),
        col="gray94",
        freq=F,
        ylim=c(0,0.5),
        main=paste("thetaAncestral"),,
        #main=paste("divergence time",pop3,"-",pop2,"(4N gen)", sep=" "),
        xlab= expression(theta[ancestral]),
        #xlab= expression(tau),
        ylab="probab density")
#hist(abc.posterior$unadj.values[,7],   breaks=seq(0,max(m.si[,8]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,7], breaks=seq(0,max(m.si[,8]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
hist(m.si[,9],
        breaks=seq(0,max(m.si[,9]+0.1),0.5),
        col="gray94",
        freq=F,
        ylim=c(0,0.5),
        main=paste("divergence time",pop4,"-",pop3,"(4N gen)", sep=" "),
        xlab= expression(tau),
        ylab="probab density")
#hist(abc.posterior$unadj.values[,7],   breaks=seq(0,max(m.si[,8]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,8], breaks=seq(0,max(m.si[,9]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
hist(m.si[,10],
        breaks=seq(0,max(m.si[,10]+0.1),0.5),
        col="gray94",
        freq=F,
        ylim=c(0,0.5),
        main=paste("divergence time",pop3,"-",pop2,"(4N gen)", sep=" "),
        xlab= expression(tau),
        ylab="probab density")
wtd.hist(abc.posterior$adj.values[,9], breaks=seq(0,max(m.si[,10]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
hist(m.si[,11],
        breaks=seq(0,max(m.si[,11]+0.1),0.5),
        col="gray94",
        freq=F,
        ylim=c(0,0.5),
        main=paste("divergence time",pop2,"-",pop1,"(4N gen)", sep=" "),
        xlab= expression(tau),
        ylab="probab density")
#hist(abc.posterior$unadj.values[,7],   breaks=seq(0,max(m.si[,8]+0.1) ,0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,10], breaks=seq(0,max(m.si[,11]+0.1) ,0.5), col=rgb(1,0.3,0,0.6), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
dev.off()
