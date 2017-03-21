#####################
#ABC for each river #
#####################
#clear R session
rm(list=ls())
ls()

library(abc)

#read the data
nlinesFul=as.numeric(strsplit(system("wc -l M.IM", intern=T), " ")[[1]][1])
M_AM=matrix(scan("M.AM"), byrow=T, nrow=nlinesFul) 
M_PAN=matrix(scan("M.PAN"), byrow=T, nrow=nlinesFul) 
M_IM=matrix(scan("M.IM"), byrow=T, nrow=nlinesFul) 
M_SC=matrix(scan("M.SC"), byrow=T, nrow=nlinesFul) 
M_I=matrix(scan("M.I"), byrow=T, nrow=nlinesFul) 


#Replace missing data by means
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


#Prepare full model
M.full=rbind(I.stat2, IM.stat2, AM.stat2, SC.stat2, PAN.stat2)
colnames(M.full) <- c("model","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")
#summary(M.full)

obs=read.table("target_sumstats.txt", header=T) #upload table of observed data
#drop out Ho, NA,(NA=Number of allele not missing data!!):
obs2=obs[,-c(1:6,13:18)]
length(obs2)


##############################
####ABC model selection ######
##############################




#############################################################
#######
#######Performing Cross validation###########################
#######
#############################################################
###edit to perform cross-validation one simply has to use a set of simulated data for each model (eg: 1000 randomly chosen data from the model pan) and to run the model selection procedure again
#example for a pairwise comparison 
#we test how well the simulated data under PAN model is classified when compared to the AM model
M.full=rbind(PAN.stat2, AM.stat2)
targetPAN=PAN.stat2[ 1 ,c(2:29)]
M.full1 = M.full[-c( 1 ),]

mod_selection_result1 <- postpr(target=targetPAN, index = as.vector(M.full1[,1]), sumstat = M.full1[,c(2:29)],tol  = 500/2e6,  method  = "neuralnet", numnet=50, sizenet=15) 
targetPAN=PAN.stat2[ 2 ,c(2:29)]
mod_selection_result2 <- postpr(target=targetPAN2, index = as.vector(M.full1[,1]), sumstat = M.full1[,c(2:29)],tol  = 500/2e6,  method  = "neuralnet", numnet=50, sizenet=15) 
#etc....

sink("crossvalPANvsAM.txt")
summary(mod_selection_result1  )
summary(mod_selection_result2 )
.....
summary(mod_selection_result1000)
sink()

#You can then put this in a loop to perform 1000 times the analysis
#I personnaly prefer using a R or bash script that will automatically generate a R script to perform 1000 model_selection and is largely faster than a R loop.
# /!\ This must be repeated for each pairwise comparison PAN model vs AM model, AM model vs PAN model, IM vs PAN, PAN vs IM, SC vs IM, and so on....



###########################################################
#############ABC parameter estimation ####################"
library(abc)

mod_names=c("model","mean_H_pop1","var_H_pop1","mean_H_pop2","var_H_pop2","mean_H_total","var_H_total","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_A_pop1","var_A_pop1","mean_A_pop2","var_A_pop2","mean_A_total","var_A_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")
length(mod_names)
colnames(I.stat2)=c("model","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")

colnames(AM.stat2)=c("model","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")

colnames(IM.stat2)=c("model","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")

colnames(SC.stat2)=c("model","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")

colnames(PAN.stat2)=c("model","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")


abc.postSC <- abc(target  = obs2, param = M_SC[,c(2:8)], sumstat = SC.stat2[,c(2:17)], tol = 1000/5e6, transf=c("logit", "logit", "logit", "logit","logit", "logit", "logit" ), 
	logit.bounds = rbind(range(M_SC[, 2]), range(M_SC[, 3]), range(M_SC[, 4]), range(M_SC[, 5]),range(M_SC[, 6]),range(M_SC[, 7]), range(M_SC[, 8])),  hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)
	
abc.postIM <- abc(target  = obs2, param = M_IM[,c(2:7)], sumstat = IM.stat2[,c(2:17)], tol = 1000/5e6, transf=c("logit", "logit", "logit", "logit","logit", "logit" ), 
	logit.bounds = rbind(range(M_IM[, 2]), range(M_IM[, 3]), range(M_IM[, 4]), range(M_IM[, 5]),range(M_IM[, 6]),range(M_IM[, 7])),  hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)

#Compute mean, mode, 95%HPD
sink("summary_posteriorABC.SC5000.txt")
print(summary(abc.postSC))
sink()


sink("summary_posteriorABC.IM5000.txt")
print(summary(abc.postIM))
sink()


#Plot some of the parameter estimation #####################
####################
# Secondary Contact#
####################

pdf(file="SC_figureHist5000.pdf",width=20,height=16)
par(mfrow=c(3,3))

hist(M_SC[,2],
     breaks=seq(0,4,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="theta1/thetaRef",
     xlab=expression((theta[1])),
     ylab="probability density")
hist((abc.postSC$unadj.values[,1] ), breaks=seq(0,4,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.postSC$adj.values[,1]), breaks=seq(0,4,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postSC$weights)
box()
#bre=expression(italic("AA river"))
#text(8,6,bre,cex=1)
legend(
 legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
pch = c(15,15,15),
x = "topleft",
cex = 1,
bty ="n")
box()

hist(M_SC[,3],
     breaks=seq(0,4,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="theta2/thetaRef",
     xlab=expression((theta[2])),
     ylab="probability density")
hist((abc.postSC$unadj.values[,2] ), breaks=seq(0,4,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.postSC$adj.values[,2]), breaks=seq(0,4,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postSC$weights)
#bre=expression(italic("Bresle river"))
#text(8,6,bre,cex=1)
#legend(
# legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
box()

hist(M_SC[,4],
     breaks=seq(0,4,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="thetaAnc/thetaRef",
     xlab=expression((theta[Ancien])),
     ylab="probability density")
hist(abc.postSC$unadj.values[,3], breaks=seq(0,4,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.postSC$adj.values[,3], breaks=seq(0,4,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postSC$weights/sum(abc.postSC$weights))
#bre=expression(italic("Bresle river"))
#text(8,6,bre,cex=1)
#legend( legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
box()

hist(M_SC[,5],
	breaks=seq(0,20,0.05),
	col="lightgrey",
	freq=F,
	ylim=c(0,2),
	main="Effective migration rate (from Lp to Lf)",
	xlab= "M1<-2",
	ylab="probability density")
hist(abc.postSC$unadj.values[,4], breaks=seq(0,20,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.postSC$adj.values[,4], breaks=seq(0,20,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postSC$weights/sum(abc.postSC$weights))
#legend( legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
box()

hist(M_SC[,6],
	breaks=seq(0,20,0.05),
	col="lightgrey",
	freq=F,
	ylim=c(0,2),
	main="Effective migration rate (from Lf to Lp)",
	xlab= "M1<-2",
	ylab="probab density")
hist(abc.postSC$unadj.values[,5], breaks=seq(0,20,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.postSC$adj.values[,5], breaks=seq(0,20,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postSC$weights/sum(abc.postSC$weights))
#legend( legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
box()
                  
hist(M_SC[,7],
	breaks=seq(0,25,0.05),
	col="lightgrey",
	freq=F,
	ylim=c(0,2),
	main="divergence time( 4N generations)",
	xlab= expression(tau),
	ylab="probab density")
hist(abc.postSC$unadj.values[,6], breaks=seq(0,25,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.postSC$adj.values[,6], breaks=seq(0,25,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postSC$weights/sum(abc.postSC$weights))
#legend( legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
box()

hist(M_SC[,8],
	breaks=seq(0,25,0.05),
	col="lightgrey",
	freq=F,
	ylim=c(0,1.5),
	main="Time Secondary Contact( 4N generations)",
	xlab= expression(tau[SC]),
	ylab="probab density")
hist(abc.postSC$unadj.values[,7], breaks=seq(0,25,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.postSC$adj.values[,7], breaks=seq(0,25,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postSC$weights/sum(abc.postSC$weights))
#legend( legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
box()
dev.off()


###### Isolation with migration ########################################
pdf(file="IM_figureHist5000.pdf",width=16,height=12)
par(mfrow=c(3,2))
hist(M_IM[,2],
     breaks=seq(0,4,0.05),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="theta1/thetaRef",
     xlab=expression((theta[1])),
     ylab="probability density")
hist((abc.postIM$unadj.values[,1] ), breaks=seq(0,4,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.postIM$adj.values[,1]), breaks=seq(0,4,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postIM$weights)
box()
legend(
 legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
pch = c(15,15,15),
x = "topleft",
cex = 1,
bty ="n")
box()

hist(M_IM[,3],
     breaks=seq(0,4,0.05),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="theta2/thetaRef",
     xlab=expression((theta[2])),
     ylab="probability density")
hist((abc.postIM$unadj.values[,2] ), breaks=seq(0,4,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.postIM$adj.values[,2]), breaks=seq(0,4,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postIM$weights)
#bre=expression(italic("Bresle river"))
#text(8,6,bre,cex=1)
#legend(
# legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
box()

hist(M_IM[,4],
     breaks=seq(0,4,0.05),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="thetaAnc/thetaRef",
     xlab=expression((theta[Anc])),
     ylab="probability density")
hist(abc.postSC$unadj.values[,3], breaks=seq(0,4,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.postSC$adj.values[,3], breaks=seq(0,4,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postIM$weights/sum(abc.postIM$weights))
#bre=expression(italic("Bresle river"))
#text(8,6,bre,cex=1)
#legend( legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
box()

hist(M_IM[,5],
	breaks=seq(0,20,0.05),
	col="lightgrey",
	freq=F,
	ylim=c(0,2.5),
	main="Effective migration rate (from Lp to Lf)",
	xlab= "M1<-2",
	ylab="probability density")
hist(abc.postIM$unadj.values[,4], breaks=seq(0,20,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.postIM$adj.values[,4], breaks=seq(0,20,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postIM$weights/sum(abc.postIM$weights))
#legend( legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
box()

hist(M_IM[,6],
	breaks=seq(0,20,0.05),
	col="lightgrey",
	freq=F,
	ylim=c(0,2.5),
	main="Effective migration rate (from Lf to Lp)",
	xlab= "M1<-2",
	ylab="probab density")
hist(abc.postIM$unadj.values[,5], breaks=seq(0,20,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.postIM$adj.values[,5], breaks=seq(0,20,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postIM$weights/sum(abc.postIM$weights))
#legend( legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
box()
                  
hist(M_IM[,7],
	breaks=seq(0,25,0.05),
	col="lightgrey",
	freq=F,
	ylim=c(0,2.5),
	main="divergence time( 4N generations)",
	xlab= expression(tau),
	ylab="probab density")
hist(abc.postIM$unadj.values[,6], breaks=seq(0,25,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.postIM$adj.values[,6], breaks=seq(0,25,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.postIM$weights/sum(abc.postIM$weights))
#legend( legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
box()

dev.off()







