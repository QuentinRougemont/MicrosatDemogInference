#########################################################################
# full script to perform simulations based on ISOLATION  MODEL using ms#
########################################################################
#lines 1 to 365 of the 5 scripts (namely SimulI, SimulAM, SimimSC, SimulPAN, SimulIM) are almoste the same, 
#they differ only on wether a prior for ancient migration, secondary contact or ongoing migration is used.
#the main difference between the 5 scripts are around line 365-450 where the call to ms is done and the model is explicitly formulated.

### Remove last features
rm(list=ls())
ls()

#install.packages("pegas")
library(pegas)
#nsim <- 1000000 #if no parallelisation is done we set nsim = 1 milion
nsim <- 10000 #otherwise we will run 100*1e4=1simulation in parallel for greater efficiency
add_simulations <- F
reftable_file <- "reference_table.txt"
batch <- 1

# Uses RR Hudson coalescent simulator (ms)
source("/usr/local/bioinfo/src/ms/msdir/readms.output.R") #set the path to the readms.output which will read result from ms in R
#source("/home/ubuntu/R/ABC/msdir/readms.output.R")
 

# Data input file consist of a text file for each population which contains allele data (size in base pairs)
# loci are in columns
# gene copies (gene copies = 2 * number of individuals) in rows
# no header, no row nameas, no column names

data_pop1 <- read.table("pop1aa") #River lamprey (example for one river)
data_pop2 <- read.table("pop2aa") #Brook lamprey

#remove missing data
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


# Setting dataset characteristics (13 loci here)
loci_names   <- c("LP-003",
                    "LP-006","LP-009","LP-018",
                  "LP-022","LP-027","LP-028","LP-30",
                  "LP-037","LP-039","LP-043","LP-045","LP-046"  
)
if (length(loci_names)!=dim(data_pop1)[2]) print("ERROR IN NUMBER OF LOCI")
motif_size <- list( 2,2,2,1,3,3,3,3,2,3,2,3,3) #for GSM

sample_size_total <- colSums(rbind(sample_size_pop1,sample_size_pop2))

min.n.pop1=min(sample_size_pop1) #used for Ar computations
min.n.pop2=min(sample_size_pop2) #Ar computation
min.n.total=min(sample_size_total) #Ar computation
#########################################
# TARGET SUMMARY STATISTICS
#########################################
# Calculate summary statistics for observed data
#Prepare the table
write.table( cbind( "mean_H_pop1", #Observed Heterozygosity
                    "var_H_pop1", "mean_H_pop2", "var_H_pop2", "mean_H_total", "var_H_total",
                    "mean_He_pop1", #Expected Heterozygosity
                     "var_He_pop1","mean_He_pop2", "var_He_pop2", "mean_He_total", "var_He_total",
                    "mean_A_pop1", #number of alleles
                    "var_A_pop1", "mean_A_pop2", "var_A_pop2", "mean_A_total", "var_A_total",
                    "mean_Ar_pop1", #Allelic richness
                    "var_Ar_pop1", "mean_Ar_pop2",  "var_Ar_pop2",  "mean_Ar_total",  "var_Ar_total",
                    #"mean_V_pop1", #Variance in allele size (bp)
                    #"var_V_pop1",
                    #"mean_V_pop2",
                    #"var_V_pop2",
                    #"mean_V_total",
                    #"var_V_total",
                    "mean_R_pop1", #allelic range
                      "var_R_pop1", "mean_R_pop2",  "var_R_pop2",  "mean_R_total", "var_R_total",
                    "mean_GW_pop1", #Garza Williamson Index
                    "var_GW_pop1",  "mean_GW_pop2", "var_GW_pop2", "mean_GW_total","var_GW_total",
                    #"mean_P_pop1", #Polymorphism
                    #"mean_P_pop2",
                    #"mean_P_total",
                    #"mean_HV_pop1", #Heterozygotie/Variance
                    #"var_HV_pop1",
                    #"mean_HV_pop2",
                    #"var_HV_pop2",
                    #"mean_HV_total",
                    #"var_HV_total",
                    "mean_GST", #Gst
                    "var_GST",
                    "mean_deltamu2",#Goldstein delta µ2
                    "var_deltamu2"#,  
),
file="target_sumstats.txt",sep=" ",
quote=F,col.names=F,row.names=F,append=F)

#Ultimately, I decided not ot use the P, HV and V as the RandomFOorest indicated they were not very informative.
#They are commented and can be used by anyone but will be probably only noise...

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
#V_pop1   <- array(NA,number_of_loci)
#V_pop2   <- array(NA,number_of_loci)
#V_total  <- array(NA,number_of_loci)
R_pop1   <- array(NA,number_of_loci)
R_pop2   <- array(NA,number_of_loci)
R_total  <- array(NA,number_of_loci)
GW_pop1  <- array(NA,number_of_loci)
GW_pop2  <- array(NA,number_of_loci)
GW_total <- array(NA,number_of_loci)
#P_pop1   <- array(F ,number_of_loci)
#P_pop2   <- array(F ,number_of_loci)
#P_total  <- array(F ,number_of_loci)
#HV_pop1  <- array(NA,number_of_loci)
#HV_pop2  <- array(NA,number_of_loci)
#HV_total <- array(NA,number_of_loci)
GST      <- array(NA,number_of_loci)
deltamu2 <- array(NA,number_of_loci)

#Compute the statistics
for (locus in 1:number_of_loci){
 
  
  H_pop1[locus]  <- H(as.factor(data_pop1[,locus]))  
  H_pop2[locus]  <- H(as.factor(data_pop2[,locus])) 
  H_total[locus] <- H(as.factor(data_total[,locus]))
  He_pop1[locus] <- 1-sum((table(as.factor(data_pop1[,locus]))/sum(table(as.factor(data_pop1[,locus])),na.rm=T))^2) 
  He_pop2[locus] <- 1-sum((table(as.factor(data_pop2[,locus]))/sum(table(as.factor(data_pop2[,locus])),na.rm=T))^2) 
  He_total[locus] <- 1-sum((table(as.factor(data_total[,locus]))/sum(table(as.factor(data_total[,locus])),na.rm=T))^2)            
  A_pop1[locus]  <- length(levels(as.factor(data_pop1[,locus]))) 
  A_pop2[locus]  <- length(levels(as.factor(data_pop2[,locus]))) 
  A_total[locus] <- length(levels(as.factor(data_total[,locus]))) 
  Ar_pop1[locus] <-sum(1-choose(sample_size_pop1[locus]-table(as.factor(data_pop1[,locus])),min.n.pop1)/choose(sample_size_pop1[locus],min.n.pop1)) 
  Ar_pop2[locus] <-sum(1-choose(sample_size_pop2[locus]-table(as.factor(data_pop2[,locus])),min.n.pop2)/choose(sample_size_pop2[locus],min.n.pop2)) 
  Ar_total[locus] <-sum(1-choose(sample_size_total[locus]-table(as.factor(data_total[,locus])),min.n.total)/choose(sample_size_total[locus],min.n.total)) 
  #V_pop1[locus]  <- var(data_pop1[,locus]) 
  #V_pop2[locus]  <- var(data_pop2[,locus]) 
  #V_total[locus] <- var(data_total[,locus]) 
  R_pop1[locus]  <- max(data_pop1[,locus])-min(data_pop1[,locus]) 
  R_pop2[locus]  <- max(data_pop2[,locus])-min(data_pop2[,locus]) 
  R_total[locus] <- max(data_total[,locus])-min(data_total[,locus]) 
  GW_pop1[locus]  <- A_pop1[locus]/(R_pop1[locus]+1) 
  GW_pop2[locus]  <- A_pop2[locus]/(R_pop2[locus]+1) 
  GW_total[locus] <- A_total[locus]/(R_total[locus]+1) 
  #if (A_pop1[locus]>1)  P_pop1[locus]  <- T 
  #if (A_pop2[locus]>1)  P_pop2[locus]  <- T 
  #if (A_total[locus]>1) P_total[locus] <- T 
  #HV_pop1[locus]  <- H_pop1[locus]/V_pop1[locus]  
  #HV_pop2[locus]  <- H_pop2[locus]/V_pop2[locus]  
  #HV_total[locus] <- H_total[locus]/V_total[locus] 
  GST[locus] <- (H_total[locus] - mean(cbind(H_pop1[locus],H_pop2[locus])))/H_total[locus]
  deltamu2[locus] <- (mean(data_pop1[,locus])-mean(data_pop2[,locus]))^2
}

for (stat in c("H","He","A","Ar","R","GW")){   #"A","V","HV","P"
  for (subsample in c("pop1","pop2","total")){ #Compute mean and var of individual statistics
    assign(paste("mean",stat,subsample,sep="_"),mean(get(paste(stat,subsample,sep="_")),na.rm=T))
    assign(paste( "var",stat,subsample,sep="_"), var(get(paste(stat,subsample,sep="_")),na.rm=T))
  }
}
mean_GST      <- mean(      GST, na.rm=T ) #Compute mean and variance of pairwise statisitics
var_GST       <-  var(      GST, na.rm=T )
mean_deltamu2 <- mean( deltamu2, na.rm=T )
var_deltamu2  <-  var( deltamu2, na.rm=T )

#Write the summary statisticts from observed data
write.table( cbind( mean_H_pop1, var_H_pop1, mean_H_pop2, var_H_pop2, mean_H_total, var_H_total,
                    mean_He_pop1, var_He_pop1, mean_He_pop2, var_He_pop2,  mean_He_total, var_He_total,
                    mean_A_pop1, var_A_pop1, mean_A_pop2, var_A_pop2, mean_A_total, var_A_total, 
                    mean_Ar_pop1, var_Ar_pop1, mean_Ar_pop2, var_Ar_pop2, mean_Ar_total, var_Ar_total,
                    #mean_V_pop1,
                    #var_V_pop1,
                    #mean_V_pop2,
                    #var_V_pop2,
                    #mean_V_total,
                    #var_V_total,
                    mean_R_pop1, var_R_pop1, mean_R_pop2, var_R_pop2, mean_R_total, var_R_total,
                    mean_GW_pop1, var_GW_pop1, mean_GW_pop2, var_GW_pop2, mean_GW_total, var_GW_total,
                    #mean_P_pop1,
                    #mean_P_pop2,
                    #mean_P_total,
                    #mean_HV_pop1,
                    #var_HV_pop1,
                    #mean_HV_pop2,
                    #var_HV_pop2,
                    #mean_HV_total,
                    #var_HV_total,
                    mean_GST, var_GST,
                    mean_deltamu2,
                    var_deltamu2#,  
),
file="target_sumstats.txt",sep=" ",
quote=F,col.names=F,row.names=F,append=T)


###################################
#Prepare reference table for simulation #
##################################
#table where the output of the simulation will be written

if (!add_simulations){
  write.table( cbind( #### PARAMETERS
    "model",      # 0:isolation; 1:isolation with migration 2: AM- Ancient Migration 3: SC Secondary Contact
    "theta1",     # theta 1   (present)
    "theta2",     # theta 2   (present)
    "thetaA",     # theta ancestral
    #"mig12",      # migration (l. planeri(2) -> L. fluviatilis (1))
    #"mig21",      # migration (lf -> lp)
    #"tau_eM",     #time of secondary contact
    "tau",        # split time
    #### SUMMARY STATISTICS  (Same as those for observed data)
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
    #"mean_V_pop1",  
    #"var_V_pop1",
    #"mean_V_pop2",
    #"var_V_pop2",
    #"mean_V_total",
    #"var_V_total",
    "mean_R_pop1",
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
    #"mean_P_pop1",
    #"mean_P_pop2",
    #"mean_P_total",
    #"mean_HV_pop1",
    #"var_HV_pop1",
    #"mean_HV_pop2",
    #"var_HV_pop2",
    #"mean_HV_total",
    #"var_HV_total",
    "mean_GST",
    "var_GST",
    "mean_deltamu2",
    "var_deltamu2"#,  
  ),
  file=reftable_file,sep=" ",
  quote=F,col.names=F,row.names=F,append=F)
 
} #END if(!add_simulations)

#########################################
# SIMULATION
#########################################
#Briefly here are the main line for each of the 5 scenarios: each parameter in [] is a variable. 
#IM & I: ms [total sample size] 1 -t [theta ref] -I 2 [sample size pop 1] [sample size pop 2] 0 -m 1 2 [mig12] -m 2 1 [mig21] -n 1 [theta 1] -n 2 [theta 2] -ej [tau] 2 1 -eN [tau] [theta A]
#SC: ms [total sample size] 1 -t [theta ref ] -I 2 [sample size pop 1] [sample size pop 2] 0 -m 1 2 [mig12] -m 2 1 [mig21] -n 1 [theta 1] -n 2 [theta 2] -eM [tau_eM] 0 -ej [tau] 2 1 -eN [tau] [theta A]
#AM: ms [total sample size] 1 -t [theta ref ] -I 2 [sample size pop 1] [sample size pop 2] 0 -ema [tau_eM] 2 0 [mig 12] [mig 21] 0 -n 1 [theta 1] -n 2 [theta 2] -ej [tau] 2 1 -eN [tau] [theta A]
#PAN:ms [total sample size] 1 -t [theta 1]


for(i in 1:nsim){
   
  cat(paste(i,"of",nsim))
 
 
   thetaRef=1 #hypothesis (possible): Ne=1000 mu=2.5e-4   with theta = 4*Nref*µ
  theta1_min <- 0 
  theta1_max <- 3 
  
  model <- 0 #index of the model
  #0: I 
  #1: IM
  #2: AM
  #3: SC
  #4: PAN
 
  theta1=NULL #theta lf
  theta2=NULL #theta lp lp < lf !
  thetaA=NULL #theta anc

  # take parameter values from priors
  theta1 = runif(1, theta1_min, theta1_max) 
  theta2=runif(1,theta1_min, theta1)
  thetaA <- runif( 1, min=theta1_min, max=theta1_max)

  mig21<-mig12<-0

  tmin=0
  tmax=25
  #tau=tsplit/4nref
  tau = runif(1, tmin, tmax)
 
  
 #Psni  <- runif( 1, min=0, max=0.1 ) #Usefull for people that want to implement single nucleotide insertion-deletion as in DIY-ABC (Cornuel et al. 2008)
 #this allow to model allele size variation
 #(contact me for detailed modification of the computations that are not here)
 alpha <- runif( 1, min=0.5, max=1 ) #geomoetrical parameter controlling the GSM  (generalised stepwise mutation model in which changes of more than one repeat unit are allowed) 
  
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
 #V_pop1   <- array(NA,number_of_loci)
 #V_pop2   <- array(NA,number_of_loci)
 #V_total  <- array(NA,number_of_loci)
 R_pop1   <- array(NA,number_of_loci)
 R_pop2   <- array(NA,number_of_loci)
 R_total  <- array(NA,number_of_loci)
 GW_pop1  <- array(NA,number_of_loci)
 GW_pop2  <- array(NA,number_of_loci)
 GW_total <- array(NA,number_of_loci)
 #P_pop1   <- array(F ,number_of_loci)
 #P_pop2   <- array(F ,number_of_loci)
 #P_total  <- array(F ,number_of_loci)
 #HV_pop1  <- array(NA,number_of_loci)
 #HV_pop2  <- array(NA,number_of_loci)
 #HV_total <- array(NA,number_of_loci)
 GST      <- array(NA,number_of_loci)
 deltamu2 <- array(NA,number_of_loci)
  
  
  for (locus in 1:number_of_loci){
    
   
    ms_out_file <- paste( "msout_locus_", locus,"_batch_",batch,".txt", sep="" )
    ms_run <- paste(sample_size_total[locus], "1",              # total sample size
                     "-t", theta1,                              # theta of reference (arbitrarilly fixed)
                     "-I 2", sample_size_pop1[locus],           
                     sample_size_pop2[locus],                   # sample sizes per population
                     "-ma X", mig12, mig21, "X",                # migration matrix
                     "-n 1", theta1,                            # size population 1
                     "-n 2", theta2,                            # size population 2
                     "-ej", tau, "2 1",                         # split population 1 and 2
                     "-eN" , tau, thetaA,                       # size ancestral population at split time
                     "-seeds","${JOB_ID}","${SGE_TASK_ID}","${SGE_TASK_ID}",  #random seed for parallelisation (this is very important to perform different simuls)
                     ">", ms_out_file)                          # output file
    
    if(.Platform$OS.type == "unix") {
      #ms_run <- paste( "/home/ubuntu/R/ABC/msdir/./ms", ms_run, sep=" " ) 
     ms_run <- paste( "/usr/local/bioinfo/src/ms/msdir/./ms", ms_run, sep=" ") #path to ms on genotoul Cluster
      
      system( ms_run ) #run ms command line
    }else{
      shell( ms_run )
    }
    
    #add the mutationns to ms output (microsatellites transformations)   
    mutations <- NA
    msout <- read.ms.output(file=ms_out_file) # read ms output file (see ms users guide for more information)
    mutations <- msout$segsites #take the number of segregating sites
    
    #sni <- rbinom(1,mutations,Psni) #for insertion-deletion we modelled sni according to a binomial distribution
    
    if  (length(motif_size[[locus]])==1){ #if only 1 mutational 
      mutation_size <- (rgeom(mutations,alpha)+1) * motif_size[[locus]] #if one wants to add sni, this is done by inserting them in the mutation as: mutations-sni
    }else{ #otheriwse take into accont the motif size /modelled prob of change according to a GSM
      mutation_size <- (rgeom(mutations,alpha)+1) * sample(motif_size[[locus]],mutations,replace=T) #here also we have to replace mutations as mutations-sni for modelling sni
    }
    #mutation_size <- c(array(1,sni),mutation_size) #usefull only for inserting the sni
    mutation_size <- mutation_size * sample(c(-1,1),mutations,replace=T) #final number of mutations will take positive or negative values (add or substract mutations)
    
    genotypes <- 500+colSums(t(msout$gametes[[1]])*mutation_size) #create the final genotypic matrix, using the gametes output (haplotype array of the data) of ms, we add a randomly chosen value (here 500 to obtain positive length of microsat)
    
    
    #calculating summary statistics for simulated data
    #these are the same formulae, adaptated for the microsatellite format derived from the transformation
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
    #V_pop1[locus]  <- var(genotypes[1:sample_size_pop1[locus]])
    #V_pop2[locus]  <- var(genotypes[(sample_size_pop1[locus]+1):sample_size_total[locus]])
    #V_total[locus] <- var(genotypes)
    R_pop1[locus]  <- max(genotypes[1:sample_size_pop1[locus]])-min(genotypes[1:sample_size_pop1[locus]])
    R_pop2[locus]  <- max(genotypes[(sample_size_pop1[locus]+1):sample_size_total[locus]])-min(genotypes[(sample_size_pop1[locus]+1):sample_size_total[locus]])
    R_total[locus] <- max(genotypes)-min(genotypes)
    GW_pop1[locus]  <- A_pop1[locus]/(R_pop1[locus]+1)
    GW_pop2[locus]  <- A_pop2[locus]/(R_pop2[locus]+1)
    GW_total[locus] <- A_total[locus]/(R_total[locus]+1)
    #if (A_pop1[locus]>1)  P_pop1[locus]  <- T
    #if (A_pop2[locus]>1)  P_pop2[locus]  <- T
    #if (A_total[locus]>1) P_total[locus] <- T
    #HV_pop1[locus]  <- H_pop1[locus]/V_pop1[locus]
    #HV_pop2[locus]  <- H_pop2[locus]/V_pop2[locus]
    #HV_total[locus] <- H_total[locus]/V_total[locus]
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
    model,	    # 0, 1, 2, 3, 
    theta1,         # theta 1   (present)
    theta2,         # theta 2   (present)
    thetaA,         # theta ancestral
    #mig12,          # migration (lp(2) -> lf(1))
    #mig21,          # migration (lf -> lp)
    tau,            # split time
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
    #mean_V_pop1,
    #var_V_pop1,
    #mean_V_pop2,
    #var_V_pop2,
    #mean_V_total,
    #var_V_total,
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
    #mean_P_pop1,
    #mean_P_pop2,
    #mean_P_total,
    #mean_HV_pop1,
    #var_HV_pop1,
    #mean_HV_pop2,
    #var_HV_pop2,
    #mean_HV_total,
    #var_HV_total,
    mean_GST,
    var_GST,
    mean_deltamu2,
    var_deltamu2#,  
  ),
  file=reftable_file,sep=" ",
  quote=F,col.names=F,row.names=F,append=T)
}
