#!/usr/bin/Rscript

#source('/home/ubuntu/abc/01.script/cv4abc.R')
#source('/home/quentin/script/rscript/cv4abc.R')
source('00-scripts/rscript/cv4abc.R')

argv<-commandArgs(TRUE)

lowbnd  <- as.numeric(argv[1])
numlines <- 200 #100
uperbnd <- lowbnd + numlines - 1
target  <- argv[2]

model.1 <-argv[3]
model.2 <-argv[4]
#model.1 <- "02.results/ABC.I.txt"
#model.2 <- "02.results/ABC.IM.txt"

nlinesFul=as.numeric(strsplit(system("wc -l 02-results/si.simul.ABC.txt", intern=T), " ")[[1]][1])

#load simulations
im=matrix(scan(model.1), byrow=T, nrow=nlinesFul) 
si=matrix(scan(model.2), byrow=T, nrow=nlinesFul) 

means <- colMeans(si, na.rm=TRUE)
for (j in 1:ncol(si)){ si[is.na(si[, j]), j] <- means[j] }
means1 <- colMeans(im, na.rm=TRUE)
for (j in 1:ncol(im)){ im[is.na(im[, j]), j] <- means1[j] }

#exclude parameters (contained in the first columns) of the simulation table to keep only summary statistics
param.im=c(2:7)
param.si=c(2:5)

#Prepare full model
M.full=rbind(si[-c(lowbnd:uperbnd),-param.si],im[-c(lowbnd:uperbnd),-param.im])
#target.stat=c(8:13,20:41)

nlineful2 <- nlinesFul-numlines

# Comparaison All
x <- as.factor(rep(c(1:2),each=nlineful2))
target.stat=c(8:13,20:25,38:41)

z <- M.full[,target.stat]

targetdata = eval(parse(text = target))
obs = targetdata[c(lowbnd:uperbnd),target.stat]

res=model_selection_abc_nnet(target=obs, x=as.vector(M.full[,1]), sumstat=z, tol=500/(nrow(M.full)), noweight=F, rejmethod=F, nb.nnet=50, size.nnet=15, output=paste("02-results/robustess_",target, "_",lowbnd,sep=""))
