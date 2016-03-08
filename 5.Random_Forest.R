#######################################################
#Performing RF analysis using RandomForestSRC packages
#######################################################
#script by Quentin Rougemont quentinrougemont@orange.fr


rm(list=ls())
library(parallel)
library(randomForestSRC)
x=detectCores()-1
options(rf.cores=x, mc.cores=x)
nlinesFul=as.numeric(strsplit(system("wc -l M.IM", intern=T), " ")[[1]][1])
M_I=matrix(scan("M.I"), byrow=T, nrow=nlinesFul) 
#AM
M_AM=matrix(scan("M.AM"), byrow=T, nrow=nlinesFul) 
#IM
M_IM=matrix(scan("M.IM"), byrow=T, nrow=nlinesFul) 
#SC
M_SC=matrix(scan("M.SC"), byrow=T, nrow=nlinesFul) 
#PAN
M_PAN=matrix(scan("M.PAN"), byrow=T, nrow=nlinesFul) 


#replace missing data by means

colSums(is.na(M_I))
means <- colMeans(M_I, na.rm=TRUE)
for (j in 1:ncol(M_I)){
     M_I[is.na(M_I[, j]), j] <- means[j]
 }
colSums(is.na(M_AM))
means2 <- colMeans(M_AM, na.rm=TRUE)
for (j in 1:ncol(M_AM)){
     M_AM[is.na(M_AM[, j]), j] <- means2[j]
 }
colSums(is.na(M_IM))
means1 <- colMeans(M_IM, na.rm=TRUE)
for (j in 1:ncol(M_IM)){
     M_IM[is.na(M_IM[, j]), j] <- means1[j]
 }
colSums(is.na(M_SC))
means3 <- colMeans(M_SC, na.rm=TRUE)
for (j in 1:ncol(M_SC)){
     M_SC[is.na(M_SC[, j]), j] <- means3[j]
 }
colSums(is.na(M_PAN)) 
means4 <- colMeans(M_PAN, na.rm=TRUE)
 for (j in 1:ncol(M_PAN)){
     M_PAN[is.na(M_PAN[, j]), j] <- means4[j]
 }

#Rename the first column according to its code (I=0, IM=1, AM=2, SC=3, PAN=4). 
M_I[which(M_I[,1]==2),1] <- 0
M_IM[which(M_IM[,1]==2),1] <- 1
M_SC[which(M_SC[,1]==2),1] <- 3
M_PAN[which(M_PAN[,1]==2),1] <-4

#exclude parameters (contained in the first columns) of the simulation table to keep only summary statistics
I.stat=M_I[,-c(2:5)]
AM.stat=M_AM[,-c(2:8)]
IM.stat=M_IM[,-c(2:7)]
SC.stat=M_SC[,-c(2:8)]
PAN.stat=M_PAN[,-c(2)]

#delete Ho, et NA (NA=Number of allele not missing data!!)
I.stat2= I.stat[,-c(2:7,14:19)]
IM.stat2=IM.stat[,-c(2:7,14:19)]
AM.stat2= AM.stat[,-c(2:7,14:19)]
SC.stat2= SC.stat[,-c(2:7,14:19)]
PAN.stat2=PAN.stat[,-c(2:7,14:19)]

#Header (not compulsory for these analysis) 
colnames(I.stat2)=c("model","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")
colnames(AM.stat2)=c("model","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")
colnames(IM.stat2)=c("model","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")
colnames(SC.stat2)=c("model","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")
colnames(PAN.stat2)=c("model","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")

#upload observed data
obs=read.table("target_sumstats.txt", header=T)
#drop out Ho, NA,(NA=Number of allele not missing data!!):
obs2=obs[,-c(1:6,13:18)]

#We will used a training set (5% of the data) 
M.full0=rbind(I.stat2[1:25000,2:29], AM.stat2[1:25000,2:29], IM.stat2[1:25000,2:29],SC.stat2[1:25000,2:29],PAN.stat2[1:25000,2:29]) #trained set 
model_name=as.factor(c(rep("I",5e4),rep("AM",5e4),rep("IM",5e4),rep("SC",5e4),rep("PAN",5e4)))
M.full1=data.frame(model_name,M.full0)
X=as.data.frame(M.full1)
attach(X)
rf.obj<-rfsrc(model_name ~., data=X, ntree = 1000,
			bootstrap = c("by.root"),importance = c("permute"), #this can be personalized, do not forget to try different methods (see manual of the RFsrc Package)
			na.action = c("na.omit"), nimpute = 1, forest=T)

sink("RF_res.txt")
print(rf.obj)
sink()

pdf(file="RF_RES.pdf", 12,12)
plot(rf.obj)
dev.off()

#Observed data
rf.pred.dat<-predict(object=rf.obj, newdata=obs2, #importance = c("permute"),
				na.action = c("na.omit"),
				outcome = c("train"), proximity = FALSE, membership = TRUE, statistics = FALSE)

sink("rf.prediction.data.txt")
print(rf.pred.dat$predicted)
print(rf.pred.dat$class)
sink()
