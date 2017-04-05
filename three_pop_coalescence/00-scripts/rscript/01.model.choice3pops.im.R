#!/usr/bin/env Rscript

#source('/home/ubuntu/abc/01.script/cv4abc.R')
source('../00-scripts/rscript/cv4abc.R')
argv<-commandArgs(TRUE)

model_type <- argv[1]

model.1 <-argv[2]
model.2 <-argv[3]
#model.3 <-argv[4]
#model.4 <-argv[5]
#odel_type <- argv[4]

#nlinesFul=as.numeric(strsplit(system("wc -l si.ssi.2ul.ABC.txt", intern=T), " ")[[1]][1])
colon_count = paste(" awk -F ' ' '{print NF }' ", model.1 , " |head -1 " )
ncol = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])
  
#load simulations
si.1=matrix(scan(model.1), byrow=T, ncol=ncol)
colon_count = paste(" awk -F ' ' '{print NF }' ", model.2 , " |head -1 " )
ncol = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])
si.2=matrix(scan(model.2), byrow=T, ncol=ncol)
#si.3=matrix(scan(model.3), byrow=T, ncol=ncol)
#si.4=matrix(scan(model.4), byrow=T, ncol=ncol)

means <- colMeans(si.1, na.rm=TRUE)
for (j in 1:ncol(si.1)){ si.1[is.na(si.1[, j]), j] <- means[j] }
means1 <- colMeans(si.2, na.rm=TRUE)
for (j in 1:ncol(si.2)){ si.2[is.na(si.2[, j]), j] <- means1[j] }
#means1 <- colMeans(si.3, na.rm=TRUE)
#for (j in 1:ncol(si.3)){ si.3[is.na(si.3[, j]), j] <- means1[j] }
#means1 <- colMeans(si.4, na.rm=TRUE)
#for (j in 1:ncol(si.4)){ si.4[is.na(si.4[, j]), j] <- means1[j] }

#exclude parameters (contained in the first columns) of the simulation table to keep only summary statistics
param.si=c(2:8)
param.im=c(2:14)

#si.1=si.1[,-param.si]
#si.2=si.2[,-param.si]
#si.3=si.3[,-param.si]
if(model_type=="si")
	{
	si.1=si.1[,-param.si]
	si.2=si.2[,-param.si]
	si.3=si.3[,-param.si]
} else if (model_type=="im") {
	si.1=si.1[,-param.im]
	si.2=si.2[,-param.im]
	#si.3=si.3[,-param.im]
	#si.4=si.4[,-param.im]
} else if (model_type=="im_si") {
	si.1=si.1[,-param.si]
	si.2=si.2[,-param.im]
}

#Prepare full model
if(model_type=="si"){
	M.full = rbind(si.1, si.2, si.3 )
} else if (model_type=="im") {
	M.full = rbind(si.1, si.2 )#, si.3, si.4)
} else if (model_type=="im_si") {
	M.full = rbind(si.1, si.2) 
}

target.stat=c(10:17,26:57,66:ncol(M.full))
target <- as.numeric(read.table("target_sumstats.txt", skip=1))

if(model_type=="si"){
	x  <- as.factor(rep(c(1:3),each=nrow(si.1)))
	} else {
	x  <- as.factor(rep(c(1:2),each=nrow(si.1)))
	}

obs <-  matrix( rep(target[-c(1:8,17:24,57:64)],10), byrow=T , nrow= 10)

model1 <- strsplit(model.1, ".ABC.txt")
model2 <- strsplit(model.2, ".ABC.txt")

res=model_selection_abc_nnet(target=obs, x=x, sumstat=M.full[,target.stat], tol=0.001, noweight=F, rejmethod=F, nb.nnet=50, size.nnet=15, output=paste("modelchoice",model1, "vs",  model2, "tol0.001", sep=".") )
