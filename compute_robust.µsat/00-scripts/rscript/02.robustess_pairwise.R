#!/usr/bin/Rscript

source('00-scripts/rscript/cv4abc.R')
#arguments
argv<-commandArgs(TRUE)

lowbnd   <- as.numeric(argv[1])
numlines <- 100
uperbnd  <- lowbnd + numlines - 1
target   <- argv[2]
target2  <- argv[3]
model.1 <-argv[4]
model.2 <-argv[5]
model.3 <-argv[6]
model.4 <-argv[7]

#count simulations
nlinesFul <- as.numeric(strsplit(system("wc -l 02-results/am.simul.ABC.txt", intern=T), " ")[[1]][1])
colon_count = paste(" awk -F ' ' '{print NF }' ", "02-results/si.simul.ABC.txt" , " |head -1 " )
ncol_si = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])
colon_count = paste(" awk -F ' ' '{print NF }' ", "02-results/am.simul.ABC.txt" , " |head -1 " )
ncol_am = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])
colon_count = paste(" awk -F ' ' '{print NF }' ", "02-results/sc.simul.ABC.txt" , " |head -1 " )
ncol_sc = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])
colon_count = paste(" awk -F ' ' '{print NF }' ", "02-results/im.simul.ABC.txt" , " |head -1 " )
ncol_im = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])
#load simulations
im <- matrix(scan(model.1), byrow=T, ncol=ncol_im) 
si <- matrix(scan(model.2), byrow=T, ncol=ncol_si) 
sc <- matrix(scan(model.3), byrow=T, ncol=ncol_sc) 
am <- matrix(scan(model.4), byrow=T, ncol=ncol_am) 

#impute missing by means (very few missing so it doesn't matter much)
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

nlinesFul=1e6
si = si[1:nlinesFul,]
am = am[1:nlinesFul,]
sc = sc[1:nlinesFul,]
im = im[1:nlinesFul,]
si = si[-c(lowbnd:uperbnd),-param.si]
im = im[-c(lowbnd:uperbnd),-param.im]
am = am[-c(lowbnd:uperbnd),-param.am]
sc = sc[-c(lowbnd:uperbnd),-param.sc]

#Prepare full model
#M.full <-rbind(si, sc ) #ici utiliser evaleparse
targetdata0 = eval(parse(text = target))
targetdata1 = eval(parse(text = target2))

M.full <-rbind( targetdata0,
		targetdata1 )
#M.full=rbind(si, im ,am, sc )
#target.stat=c(8:13,20:41)

nlineful2 <- nlinesFul-numlines
# Comparaison All
x <- as.factor(rep(c(1:2),each=nlineful2))
target.stat=c(8:13,20:41)
z <- M.full[,target.stat]

targetdata = eval(parse(text = target))

obs = targetdata[c(lowbnd:uperbnd),target.stat]

res <- model_selection_abc_nnet(target=obs, x=x, sumstat=z, tol=0.00025, noweight=F, rejmethod=F, nb.nnet=50, size.nnet=15, output=paste("02-results/robustess.",target, "_vs_", target2,lowbnd,sep=""))
