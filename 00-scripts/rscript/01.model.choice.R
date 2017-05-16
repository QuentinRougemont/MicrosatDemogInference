#!/usr/bin/Rscript

#source('/home/ubuntu/abc/01.script/cv4abc.R')
#source('/home/quentin/script/rscript/cv4abc.R')
source('../00-scripts/rscript/cv4abc.R')
argv<-commandArgs(TRUE)

model.1 <-argv[1]
model.2 <-argv[2]
model.3 <-argv[3]
model.4 <-argv[4]

nlinesFul=as.numeric(strsplit(system("wc -l si.simul.ABC.txt", intern=T), " ")[[1]][1])

#load simulations
im=matrix(scan(model.1), byrow=T, nrow=nlinesFul)
si=matrix(scan(model.2), byrow=T, nrow=nlinesFul)
sc=matrix(scan(model.3), byrow=T, nrow=nlinesFul)
am=matrix(scan(model.4), byrow=T, nrow=nlinesFul)

means <- colMeans(si, na.rm=TRUE)
for (j in 1:ncol(si)){ si[is.na(si[, j]), j] <- means[j] }
means1 <- colMeans(im, na.rm=TRUE)
for (j in 1:ncol(im)){ im[is.na(im[, j]), j] <- means1[j] }
means1 <- colMeans(sc, na.rm=TRUE)
for (j in 1:ncol(sc)){ sc[is.na(sc[, j]), j] <- means1[j] }
means1 <- colMeans(am, na.rm=TRUE)
for (j in 1:ncol(am)){ am[is.na(am[, j]), j] <- means1[j] }

#exclude parameters (contained in the first columns) of the simulation table to keep only summary statistics
param.im=c(2:7)
param.si=c(2:5)
param.sc=c(2:8)
param.am=c(2:8)

si=si[,-param.si]
im=im[,-param.im]
am=am[,-param.am]
sc=sc[,-param.sc]
#
#Prepare full model
M.full=rbind(si, im ,am, sc )
#
target.stat=c(8:13,20:41)

target <- as.numeric(read.table("target_sumstats.txt", skip=1))
#drop out Ho, NA,(NA=Number of allele not missing data!!):

x  <- as.factor(rep(c(1:4),each=nlinesFul))
#all stats
obs <-  matrix( rep(target[-c(1:6,13:18)],10), byrow=T , nrow= 10)

#res=model_selection_abc_nnet(target=obs, x=x, sumstat=M.full[,target.stat], tol=1000/(nrow(M.full)), noweight=F, rejmethod=F, nb.nnet=50, size.nnet=15, output="obs")

res=model_selection_abc_nnet(target=obs, x=x, sumstat=M.full[,target.stat], tol=0.001, noweight=F, rejmethod=F, nb.nnet=50, size.nnet=15, output="obs_tol0.001")

###Différent choix de stat:
obs2=obs2[-c(13:24)]
target.stat=c(8:13,20:25,38:41)
res=model_selection_abc_nnet(target=obs, x=x, sumstat=M.full[,target.stat], tol=0.001, noweight=F, rejmethod=F, nb.nnet=50, size.nnet=15, output="obs.16stats.tol0.001")

res=model_selection_abc_nnet(target=obs, x=x, sumstat=M.full[,target.stat], tol=1000/(nrow(M.full)), noweight=F, rejmethod=F, nb.nnet=50, size.nnet=15, output="obs.16stats")

