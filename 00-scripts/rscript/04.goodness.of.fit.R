

library(abc)
library(pegas)
library(weights)
library(laeken)
library(grid)

argv  <-commandArgs(TRUE)
model <-argv[1]
pop1  <- argv[2]
pop2  <- argv[3]
nsim  <- as.numeric(argv[4])
nloc  <- as.numeric(argv[5])
motif <- argv[6]

#iso
nlinesFul <- as.numeric(strsplit(system("wc -l ../si.simul.ABC.txt", intern=T), " ")[[1]][1])
#load simul
M_SI <- matrix(scan(model), byrow=T, nrow=nlinesFul)

means <- colMeans(M_SI, na.rm=TRUE)
for (j in 1:ncol(M_SI)){
     M_SI[is.na(M_SI[, j]), j] <- means[j]
 }

#exclude parameters (contained in the first columns) of the simulation table to keep only summary statistics
param.si=c(2:5)
target.stat=c(12:17,24:45)

obs <- read.table("../target_sumstats.txt", header=T)
#drop out Ho, NA,(NA=Number of allele not missing data!!):
obs2=obs[,-c(1:6,13:18)]

#estimate posterior
abc.posterior <- abc(target  = obs2, param = M_SI[,param.si], sumstat = M_SI[,target.stat], tol = 500/2e6, transf=c("logit", "logit", "logit", "logit"),
        logit.bounds = rbind(range(M_SI[, 2]), range(M_SI[, 3]), range(M_SI[, 4]), range(M_SI[, 5])),  hcorr=T, method  = "neuralnet", numnet=20, sizenet=5)


#sample from posterior
simpost<-nsim
newsamp<-sample(1:(dim(abc.posterior$adj)[1]),size=simpost,replace=T,prob=abc.posterior$weights) 
newsamp1<-abc.posterior$adj[newsamp,]
write.table(newsamp1,"ppc.si", quote=F,row.names=F,col.names=T)

theta1=newsamp1[,1]
theta2=newsamp1[,2]
thetaA=newsamp1[,3]
tau=newsamp1[,4]

#prepare simulation pipeline
nsim  <- nsim
add_simulations <- F
reftable_file   <- "reference_table_ppc.si.txt"
batch <- 1

# Uses RR Hudson coalescent simulator (ms)
source("/home/quentin/script/msdir/readms.output.R")

# Data input file consist on a text file for each population which contains allele data (size in base pairs)
# loci in columns
# gene copies (gene copies = 2 * number of individuals) in rows
data_pop1 <- read.table(pop1)
data_pop2 <- read.table(pop2)

#remove outlier alleles and missing data
for (i in 1:length(data_pop1)){
data_pop1[which(data_pop1[,i]==0),i]<-NA
  }

for (i in 1:length(data_pop2)){
data_pop2[which(data_pop2[,i]==0),i]<-NA
  }

# combines dataset
number_of_loci <- dim(data_pop1)[2]
sample_size_pop1 <- array(NA,number_of_loci)
for (i in 1:number_of_loci){
  sample_size_pop1[i]                  <- length(which(!is.na(data_pop1[,i])))
  data_pop1[which(data_pop1[,i]==0),i] <- NA
}

if (dim(data_pop2)[2]!=dim(data_pop1)[2]) print("ERROR IN NUMBER OF LOCI")
sample_size_pop2 <- array(NA,number_of_loci)
for (i in 1:number_of_loci){
  data_pop2[which(data_pop2[,i]==0),i] <- NA
  sample_size_pop2[i]                  <- length(which(!is.na(data_pop2[,i])))
}

data_total <- rbind(data_pop1,data_pop2)
# Setting dataset characteristics (13 loci)
#loci_names   <- c("LP-003","LP-006","LP-009","LP-018","LP-022","LP-027","LP-028","LP-30","LP-037","LP-039","LP-043","LP-045","LP-046")
loci_names  <- rep(paste("loc",seq(1,nloc),sep="."))

if (length(loci_names)!=dim(data_pop1)[2]) print("ERROR IN NUMBER OF LOCI")
#motif_size <- list( 2,2,2,1,3,3,3,3,2,3,2,3,3) #for GSM
motif_size <- as.list(as.matrix(read.table(motif)))


sample_size_total <- colSums(rbind(sample_size_pop1,sample_size_pop2))

min.n.pop1=min(sample_size_pop1) #for Ar computations
min.n.pop2=min(sample_size_pop2)
min.n.total=min(sample_size_total)


###################################
#Prepare ref.table for simulation #
##################################
#table were the output of the simulation will be written
if (!add_simulations){
  write.table( cbind( #### PARAMETERS
    "model",      # 0:isolation; 1:isolation with migration 2: AM- Ancient Migration 3: SC Secondary Contact
    "theta1",     # theta 1   (present)
    #### SUMMARY STATISTICS                      
    "mean_H_pop1",
    "var_H_pop1",
    "mean_H_pop2",
    "var_H_pop2",
    "mean_H_total",
    "var_H_total",
    "mean_He_pop1",
    "var_He_pop1",
    "mean_He_pop2",
    "var_He_pop2",
    "mean_He_total",
    "var_He_total",
    "mean_A_pop1",
    "var_A_pop1",
    "mean_A_pop2",
    "var_A_pop2",
    "mean_A_total", 
    "var_A_total",
    "mean_Ar_pop1",
    "var_Ar_pop1",
    "mean_Ar_pop2",
    "var_Ar_pop2",
    "mean_Ar_total",
    "var_Ar_total",
    "mean_V_pop1",  #variance allele size (in bp)
    "var_V_pop1",
    "mean_V_pop2",
    "var_V_pop2",
    "mean_V_total",
    "var_V_total",
    "mean_R_pop1", #allelic range
    "var_R_pop1",
    "mean_R_pop2",
    "var_R_pop2",
    "mean_R_total",
    "var_R_total",
    "mean_GW_pop1",
    "var_GW_pop1",
    "mean_GW_pop2",
    "var_GW_pop2",
    "mean_GW_total",
    "var_GW_total",
    "mean_GST",
    "var_GST",
    "mean_deltamu2",
    "var_deltamu2"#,  
  ),
  file=reftable_file,sep=" ",
  quote=F,col.names=F,row.names=F,append=F)
 
} #END if(!add_simulations)


for(i in 1:length(theta1)){
   
  cat(paste(i,"of",length(theta1)))
  thetaRef=1
  theta=theta1[i]
  thetaB=theta2[i]
  thetaC=thetaA[i]
  mig12a=0
  mig21a=0
  #tau_em1=tau_eM[i]
  tau1=tau[i]
  alpha <- runif( 1, min=0.5, max=1 )
  
 H_pop1   <- array(NA,number_of_loci)
 H_pop2   <- array(NA,number_of_loci)
 H_total  <- array(NA,number_of_loci)
 He_pop1  <- array(NA,number_of_loci)
 He_pop2  <- array(NA,number_of_loci)
 He_total <- array(NA,number_of_loci)
 A_pop1   <- array(NA,number_of_loci)
 A_pop2   <- array(NA,number_of_loci)
 A_total  <- array(NA,number_of_loci)
 Ar_pop1   <- array(NA,number_of_loci)
 Ar_pop2   <- array(NA,number_of_loci)
 Ar_total  <- array(NA,number_of_loci)
 R_pop1   <- array(NA,number_of_loci)
 R_pop2   <- array(NA,number_of_loci)
 R_total  <- array(NA,number_of_loci)
 GW_pop1  <- array(NA,number_of_loci)
 GW_pop2  <- array(NA,number_of_loci)
 GW_total <- array(NA,number_of_loci)
 GST      <- array(NA,number_of_loci)
 deltamu2 <- array(NA,number_of_loci)
  
  for (locus in 1:number_of_loci){
    #for (p in 1:length(theta1)){
    #prépare ms settings
    ms_out_file <- paste( "msout_locus_", locus,"_batch_",batch,".txt", sep="" )
    ms_run <- paste(sample_size_total[locus], "1",       # total sample size
                     "-t", thetaRef,                              # size population 1
                     "-I 2", sample_size_pop1[locus], 
                     sample_size_pop2[locus],           # sample sizes per population
                     "-ma X", mig12a, mig21a, "X",                # migration matrix
                     "-n 1", theta,                             # size population 2
                     "-n 2", thetaB,                       # size ancestral population (pop 1)
                     "-ej", tau1, "2 1",                       #migration sets to Zéro (ok for SC, but try for AM)
                     "-eN" , tau1, thetaC,                            # size population1+2
                    #"-seeds", theta1*100, theta1*5000, theta1*2500,
                     ">", ms_out_file)                          # output file
    #}
    if(.Platform$OS.type == "unix") {
      ms_run <- paste( "ms", ms_run, sep=" " ) 
      #ms_run <- paste( "/home/qrougemo/bin/ms_test/msdir/./ms",ms_run, sep=" ") #unil
      #ms_run <- paste( "/usr/local/bioinfo/src/ms/msdir/./ms", ms_run, sep=" ") #genotool
       
      system( ms_run )
    }else{
      shell( ms_run )
    }

 # add the mutationns to ms output (microsatellites transformations)   
    mutations <- NA
    msout <- read.ms.output(file=ms_out_file)
    mutations <- msout$segsites
    
    
    if  (length(motif_size[[locus]])==1){
      mutation_size <- (rgeom(mutations,alpha)+1) * motif_size[[locus]]
    }else{
      mutation_size <- (rgeom(mutations,alpha)+1) * sample(motif_size[[locus]],mutations,replace=T)
    }
  
    mutation_size <- mutation_size * sample(c(-1,1),mutations,replace=T)
     
    genotypes <- 500+colSums(t(msout$gametes[[1]])*mutation_size)
    
    #calculating summary statistics
    H_pop1[locus]  <- H(as.factor(genotypes[1:sample_size_pop1[locus]])) 
    H_pop2[locus]  <- H(as.factor(genotypes[(sample_size_pop1[locus]+1):sample_size_total[locus]]))
    H_total[locus] <- H(as.factor(genotypes))
    He_pop1[locus] <- 1-sum((table(as.factor(genotypes[1:sample_size_pop1[locus]]))/sum(table(as.factor(genotypes[1:sample_size_pop1[locus]])),na.rm=T))^2)
    He_pop2[locus] <- 1-sum((table(as.factor(genotypes[(sample_size_pop1[locus]+1):sample_size_total[locus]]))/sum(table(as.factor(genotypes[(sample_size_pop1[locus]+1):sample_size_total[locus]])),na.rm=T))^2)
    He_total[locus] <- 1-sum((table(as.factor(genotypes))/sum(table(as.factor(genotypes)),na.rm=T))^2)           
    A_pop1[locus]  <- length(levels(as.factor(genotypes[1:sample_size_pop1[locus]])))
    A_pop2[locus]  <- length(levels(as.factor(genotypes[(sample_size_pop1[locus]+1):sample_size_total[locus]])))
    A_total[locus] <- length(levels(as.factor(genotypes)))
    Ar_pop1[locus] <-sum(1-choose(sample_size_pop1[locus]-table(as.factor(genotypes[1:sample_size_pop1[locus]])),min.n.pop1)/choose(sample_size_pop1[locus],min.n.pop1))
    Ar_pop2[locus] <-sum(1-choose(sample_size_pop2[locus]-table(as.factor(genotypes[(sample_size_pop1[locus]+1):sample_size_total[locus]])),min.n.pop2)/choose(sample_size_pop2[locus],min.n.pop2))
    Ar_total[locus] <-sum(1-choose(sample_size_total[locus]-table(as.factor(genotypes)),min.n.total)/choose(sample_size_total[locus],min.n.total))
    R_pop1[locus]  <- max(genotypes[1:sample_size_pop1[locus]])-min(genotypes[1:sample_size_pop1[locus]])
    R_pop2[locus]  <- max(genotypes[(sample_size_pop1[locus]+1):sample_size_total[locus]])-min(genotypes[(sample_size_pop1[locus]+1):sample_size_total[locus]])
    R_total[locus] <- max(genotypes)-min(genotypes)
    GW_pop1[locus]  <- A_pop1[locus]/(R_pop1[locus]+1)
    GW_pop2[locus]  <- A_pop2[locus]/(R_pop2[locus]+1)
    GW_total[locus] <- A_total[locus]/(R_total[locus]+1)
    GST[locus] <- (H_total[locus] - mean(cbind(H_pop1[locus],H_pop2[locus])))/H_total[locus]
    deltamu2[locus] <- (mean(genotypes[1:sample_size_pop1[locus]])-mean(genotypes[(sample_size_pop1[locus]+1):sample_size_total[locus]]))^2
  }
  
  for (stat in c("H","He","A", "Ar","R","GW")){ #"HV","P" "V"
    for (subsample in c("pop1","pop2","total")){
      assign(paste("mean",stat,subsample,sep="_"),mean(get(paste(stat,subsample,sep="_")),na.rm=T))
      assign(paste( "var",stat,subsample,sep="_"), var(get(paste(stat,subsample,sep="_")),na.rm=T))
    }
  }
  mean_GST      <- mean(      GST, na.rm=T )
  var_GST       <-  var(      GST, na.rm=T )
  mean_deltamu2 <- mean( deltamu2, na.rm=T )
  var_deltamu2  <-  var( deltamu2, na.rm=T )
 
  write.table( cbind( ## PARAMETERS
    #model,	    # 0, 1, 2, 3, 
    #theta1,         # theta 1   (present)
    ## SUMMARY STATISTICS
    mean_H_pop1,
    var_H_pop1,
    mean_H_pop2,
    var_H_pop2,
    mean_H_total,
    var_H_total,
    mean_He_pop1,
    var_He_pop1,
    mean_He_pop2,
    var_He_pop2,
    mean_He_total,
    var_He_total,
    mean_A_pop1,
    var_A_pop1,
    mean_A_pop2,
    var_A_pop2,
    mean_A_total,
    var_A_total,
    mean_Ar_pop1,
    var_Ar_pop1,
    mean_Ar_pop2,
    var_Ar_pop2,
    mean_Ar_total,
    var_Ar_total,
    mean_R_pop1,
    var_R_pop1,
    mean_R_pop2,
    var_R_pop2,
    mean_R_total,
    var_R_total,
    mean_GW_pop1,
    var_GW_pop1,
    mean_GW_pop2,
    var_GW_pop2,
    mean_GW_total,
    var_GW_total,
    mean_GST,
    var_GST,
    mean_deltamu2,
    var_deltamu2#,  
  ),
  file=reftable_file,sep=" ",
  quote=F,col.names=F,row.names=F,append=T)
}


obs=read.table("../target_sumstats.txt", header=T)
#drop out Ho, NA,(NA=Number of allele not missing data!!):

ppc_im=read.table('reference_table_ppc.si.txt', skip=1)
mod_names=c("mean_H_pop1","var_H_pop1","mean_H_pop2","var_H_pop2","mean_H_total","var_H_total","mean_He_pop1","var_He_pop1","mean_He_pop2","var_He_pop2","mean_He_total","var_He_total","mean_A_pop1","var_A_pop1","mean_A_pop2","var_A_pop2","mean_A_total","var_A_total","mean_Ar_pop1","var_Ar_pop1","mean_Ar_pop2","var_Ar_pop2","mean_Ar_total","var_Ar_total","mean_R_pop1","var_R_pop1","mean_R_pop2","var_R_pop2","mean_R_total","var_R_total","mean_GW_pop1","var_GW_pop1","mean_GW_pop2","var_GW_pop2","mean_GW_total","var_GW_total","mean_GST","var_GST","mean_deltamu2","var_deltamu2")

vect=obs

#plot stat obs et stat sim:
pdf("PPC_I.pdf",width=40,height=40)
par(mfrow=c(6,7))
for (i in 1:ncol(ppc_im)){
	hist(ppc_im[,i], xlab=mod_names[i], main=mod_names[i])
	abline(v=obs[,i], col='red')
	pval=mean(abs(ppc_im[,i]) >= abs(vect[,i]))
	write.table(t(pval),file="ppc_I_Goodness.txt", append=T,quote=F,row.names=mod_names[i], col.names=F)
	}
dev.off()


