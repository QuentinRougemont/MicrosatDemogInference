##############################################################################
# full script to perform simulations based on SECONDARY CONTACT MODEL using ms#
##############################################################################
if("pegas" %in% rownames(installed.packages()) == FALSE) {install.packages("pegas") }

library(pegas)
argv  <- commandArgs(TRUE)

if (argv[1]=="-h" || length(argv)==0){
        cat("\n 8 parameter needed!! \n  
    (1) pop1 (genotype data for pop1) \n, 
    (2) pop2 (genotype data for pop2) \n, 
    (3) nsim (number of simulation), \n 
    (4) nloc (number of loci)\n 
    (5) pattern of repeat motif for each loci (a single colomn file with one motif on each ligne) \n 
    #(6) migration_model (for migration matrix) \n 
    #(7) model_type(either am, im, sc or si)  \n" )
}else{

pop1  <- argv[1]
pop2  <- argv[2]
nsim  <- as.numeric(argv[3])
nloc  <- as.numeric(argv[4])
motif <- argv[5]
}


#nsim <- 10000
nsim  <- nsim
add_simulations <- F
reftable_file   <- "reference_table.txt"
batch <- 1

# Uses RR Hudson coalescent simulator (ms)
source("../../00-scripts/msdir/readms.output.R")

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
motif_size <- as.list(as.matrix(read.table(motif)))

sample_size_total <- colSums(rbind(sample_size_pop1,sample_size_pop2))

min.n.pop1=min(sample_size_pop1) #for Ar computations
min.n.pop2=min(sample_size_pop2)
min.n.total=min(sample_size_total)


#########################################
# TARGET SUMMARY STATISTICS
#########################################
# Calculate summary statistics for observed data
write.table( cbind( "mean_H_pop1",
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
file="target_sumstats.txt",sep=" ",
quote=F,col.names=F,row.names=F,append=F)

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
 
  
  H_pop1[locus]  <- H(as.factor(data_pop1[,locus])) #ok 
  H_pop2[locus]  <- H(as.factor(data_pop2[,locus])) #ok
  H_total[locus] <- H(as.factor(data_total[,locus]))
  He_pop1[locus] <- 1-sum((table(as.factor(data_pop1[,locus]))/sum(table(as.factor(data_pop1[,locus])),na.rm=T))^2) #ok
  He_pop2[locus] <- 1-sum((table(as.factor(data_pop2[,locus]))/sum(table(as.factor(data_pop2[,locus])),na.rm=T))^2) #ok
  He_total[locus] <- 1-sum((table(as.factor(data_total[,locus]))/sum(table(as.factor(data_total[,locus])),na.rm=T))^2) #ok           
  A_pop1[locus]  <- length(levels(as.factor(data_pop1[,locus]))) #ok
  A_pop2[locus]  <- length(levels(as.factor(data_pop2[,locus]))) #ok
  A_total[locus] <- length(levels(as.factor(data_total[,locus]))) #ok
  Ar_pop1[locus] <-sum(1-choose(sample_size_pop1[locus]-table(as.factor(data_pop1[,locus])),min.n.pop1)/choose(sample_size_pop1[locus],min.n.pop1)) #ok
  Ar_pop2[locus] <-sum(1-choose(sample_size_pop2[locus]-table(as.factor(data_pop2[,locus])),min.n.pop2)/choose(sample_size_pop2[locus],min.n.pop2)) #ok
  Ar_total[locus] <-sum(1-choose(sample_size_total[locus]-table(as.factor(data_total[,locus])),min.n.total)/choose(sample_size_total[locus],min.n.total)) #ok
  #V_pop1[locus]  <- var(data_pop1[,locus]) #ok
  #V_pop2[locus]  <- var(data_pop2[,locus]) #ok
  #V_total[locus] <- var(data_total[,locus]) #ok
  R_pop1[locus]  <- max(data_pop1[,locus])-min(data_pop1[,locus]) #ok
  R_pop2[locus]  <- max(data_pop2[,locus])-min(data_pop2[,locus]) #ok
  R_total[locus] <- max(data_total[,locus])-min(data_total[,locus]) #ok
  GW_pop1[locus]  <- A_pop1[locus]/(R_pop1[locus]+1) #ok
  GW_pop2[locus]  <- A_pop2[locus]/(R_pop2[locus]+1) #ok
  GW_total[locus] <- A_total[locus]/(R_total[locus]+1) #ok
  #if (A_pop1[locus]>1)  P_pop1[locus]  <- T #ok
  #if (A_pop2[locus]>1)  P_pop2[locus]  <- T #ok
  #if (A_total[locus]>1) P_total[locus] <- T #ok
  #HV_pop1[locus]  <- H_pop1[locus]/V_pop1[locus] #ok 
  #HV_pop2[locus]  <- H_pop2[locus]/V_pop2[locus] #ok 
  #HV_total[locus] <- H_total[locus]/V_total[locus] #ok
  GST[locus] <- (H_total[locus] - mean(cbind(H_pop1[locus],H_pop2[locus])))/H_total[locus]
  deltamu2[locus] <- (mean(data_pop1[,locus])-mean(data_pop2[,locus]))^2
}

for (stat in c("H","He","A","Ar","R","GW")){   #"A","V","HV","P"
  for (subsample in c("pop1","pop2","total")){
    assign(paste("mean",stat,subsample,sep="_"),mean(get(paste(stat,subsample,sep="_")),na.rm=T))
    assign(paste( "var",stat,subsample,sep="_"), var(get(paste(stat,subsample,sep="_")),na.rm=T))
  }
}
mean_GST      <- mean(      GST, na.rm=T )
var_GST       <-  var(      GST, na.rm=T )
mean_deltamu2 <- mean( deltamu2, na.rm=T )
var_deltamu2  <-  var( deltamu2, na.rm=T )


write.table( cbind( mean_H_pop1,
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
file="target_sumstats.txt",sep=" ",
quote=F,col.names=F,row.names=F,append=T)


###################################
#Prepare ref.table for simulation #
##################################
#table were the output of the simulation will be written

if (!add_simulations){
  write.table( cbind( #### PARAMETERS
    "model",      # 0:isolation; 1:isolation with migration 2: AM- Ancient Migration 3: SC Secondary Contact
    "theta1",     # theta 1   (present)
    "theta2",     # theta 2   (present)
    "thetaA",     # theta ancestral
    "mig12",      # migration (l. planeri(2) -> L. fluviatilis (1))
    "mig21",      # migration (lf -> lp)
    "tau_eM",
    "tau",        # split time
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
    #"mean_V_pop1",  #variance allele size (in bp)
    #"var_V_pop1",
    #"mean_V_pop2",
    #"var_V_pop2",
    #"mean_V_total",
    #"var_V_total",
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
    #"mean_P_pop1",
    #"mean_P_pop2",
    #"mean_P_total",
    #"mean_HV_pop1", #Hetero/Variance (inutile)
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
#IM & I: ms [total sample size] 1 -t [theta ref] -I 2 [sample size pop 1] [sample size pop 2] 0 -m 1 2 [mig12] -m 2 1 [mig21] -n 1 [theta 1] -n 2 [theta 2] -ej [tau] 2 1 -eN [tau] [theta A]
#SC: ms [total sample size] 1 -t [theta ref ] -I 2 [sample size pop 1] [sample size pop 2] 0 -m 1 2 [mig12] -m 2 1 [mig21] -n 1 [theta 1] -n 2 [theta 2] -eM [tau_eM] 0 -ej [tau] 2 1 -eN [tau] [theta A]
#AM: ms [total sample size] 1 -t [theta ref ] -I 2 [sample size pop 1] [sample size pop 2] 0 -ema [tau_eM] 2 0 [mig 12] [mig 21] 0 -n 1 [theta 1] -n 2 [theta 2] -ej [tau] 2 1 -eN [tau] [theta A]
#PAN:ms [total sample size] 1 -t [theta 1]

for(i in 1:nsim){
   
  cat(paste(i,"of",nsim))
 
  thetaRef=1.2 #hypothèse: Ne=10000 mu=10e-5
  theta1_min <- 0 
  theta1_max <- 3
  thetaA_max <- 5

  model <- 3
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
  theta2=runif(1,theta1_min, theta1_max)
  thetaA <- runif( 1, min=theta1_min, max=thetaA_max)

 mig12   <- runif( 1, min=0, max=15 )
 mig21   <- runif( 1, min=0, max=15 )

 tmin = 0
 tmax = 25

 tau = runif(1, tmin, tmax)
 tau_eM = runif(1, tmin, tau)
 
  #Psni  <- runif( 1, min=0, max=0.1 )
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
    
    #prépare ms settings
    ms_out_file <- paste( "msout_locus_", locus,"_batch_",batch,".txt", sep="" )
    ms_run <- paste(sample_size_total[locus], "1",       # total sample size
                     "-t", thetaRef,                              # size population 1
                     "-I 2", sample_size_pop1[locus], 
                     sample_size_pop2[locus],           # sample sizes per population
                     "-ma X", mig12, mig21, "X",                # migration matrix
                     "-n 1", theta1,                             # size population 2
                     "-n 2", theta2,                       # size ancestral population (pop 1)
                     "-eM" , tau_eM, "0", 
                     "-ej", tau, "2 1",                       #migration sets to Zéro (ok for SC, but try for AM)
                     "-eN" , tau, thetaA,
                    # "-seeds","${JOB_ID}","${SGE_TASK_ID}","${SGE_TASK_ID}",
                     ">", ms_out_file)                          # output file
    
    if(.Platform$OS.type == "unix") {
     # ms_run <- paste( "/home/ubuntu/R/ABC/msdir/./ms", ms_run, sep=" " ) 
      #ms_run <- paste( "/home/qrougemo/bin/ms_test/msdir/./ms", ms_run, sep=" ") #unil
      ms_run <- paste( "ms", ms_run, sep=" ") #genotool
      
      system( ms_run )
    }else{
      shell( ms_run )
    }
    

   # add the mutationns to ms output (microsatellites transformations)   
    mutations <- NA
    msout <- read.ms.output(file=ms_out_file) #ou read.ms.output
    mutations <- msout$segsites
    
    #sni <- rbinom(1,mutations,Psni)
    
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
    mig12,          # migration (lp(2) -> lf(1))
    mig21,          # migration (lf -> lp)
    tau_eM,        # end of ancien migration
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
