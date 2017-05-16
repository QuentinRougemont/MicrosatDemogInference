
#purpose = "run coalescent simulation using ms for microsatellite data for a 4population model"
#date    = "05/05/2017"
#authors =  "Quentin Rougemont"
#modification from the original 2 pop scripts developped by M. Navascues
#Here allows to compare different rates and direction of migration (im model) and to compare against a model without migration
#output  = "file of empirically observed statistics (target_sumstats.txt ) and reference table of coalescent sims.(reference_table.txt)"

if("pegas" %in% rownames(installed.packages()) == FALSE) {install.packages("pegas") }

library(pegas)

argv  <- commandArgs(TRUE)

if (argv[1]=="-h" || length(argv)==0){
        cat("\n 8 parameter needed!! \n  
    (1) pop1 (genotype data for pop1) \n, 
    (2) pop2 (genotype data for pop2) \n, 
    (3) pop3 (genotype data for pop3) \n 
    (4) pop3 (genotype data for pop4) \n 
    (5) nsim (number of simulation), \n 
    (6) nloc (number of loci)\n 
    (7) pattern of repeat motif for each loci (a single colomn file with one motif on each ligne) \n 
    (8) migration_model (for migration matrix) \n 
    (9) model_type(either im or si)  \n" )
}else{

pop1  <- argv[1]
pop2  <- argv[2]
pop3  <- argv[3]
pop4  <- argv[4]
nsim  <- as.numeric(argv[5])
nloc  <- as.numeric(argv[6])
motif <- argv[7]
mig_model  <- argv[8]
model_type <- argv[9]
}

nsim  <- nsim
add_simulations <- F
reftable_file   <- "reference_table.txt"
batch <- 1

# Uses RR Hudson coalescent simulator (ms)
source("../../00-scripts/msdir/readms.output.R")
# source("/home/ubuntu/programmation/R/ABC/msdir/readms.output.R")
# READ DATA
# Data input file consist on a text file for each population which contains allele data (size in base pairs)
# loci in columns
# gene copies (gene copies = 2 * number of individuals) in rows

data_pop1 <- read.table(pop1)
data_pop2 <- read.table(pop2)
data_pop3 <- read.table(pop3)
data_pop4 <- read.table(pop4)


#remove outlier alleles and missing data
for (i in 1:length(data_pop1)){
data_pop1[which(data_pop1[,i]==0),i]<-NA
  }

for (i in 1:length(data_pop2)){
data_pop2[which(data_pop2[,i]==0),i]<-NA
  }
  
for (i in 1:length(data_pop3)){
data_pop3[which(data_pop3[,i]==0),i]<-NA
  }

for (i in 1:length(data_pop4)){
data_pop4[which(data_pop4[,i]==0),i]<-NA
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

if (dim(data_pop3)[2]!=dim(data_pop2)[2]) print("ERROR IN NUMBER OF LOCI")
sample_size_pop3 <- array(NA,number_of_loci)
for (i in 1:number_of_loci){
  data_pop3[which(data_pop3[,i]==0),i] <- NA
  sample_size_pop3[i]                  <- length(which(!is.na(data_pop3[,i])))
}

if (dim(data_pop4)[2]!=dim(data_pop3)[2]) print("ERROR IN NUMBER OF LOCI")
sample_size_pop4 <- array(NA,number_of_loci)
for (i in 1:number_of_loci){
  data_pop4[which(data_pop4[,i]==0),i] <- NA
  sample_size_pop4[i]                  <- length(which(!is.na(data_pop4[,i])))
}

data_total <- rbind(data_pop1,data_pop2,data_pop3,data_pop4)
# Setting dataset characteristics 
loci_names  <- rep(paste("loc",seq(1,nloc),sep="."))

if (length(loci_names)!=dim(data_pop1)[2]) print("ERROR IN NUMBER OF LOCI")
motif_size <- as.list(as.matrix(read.table(motif)))

sample_size_total <- colSums(rbind(sample_size_pop1,sample_size_pop2, sample_size_pop3, sample_size_pop4))

min.n.pop1 = min(sample_size_pop1) 
min.n.pop2 = min(sample_size_pop2)
min.n.pop3 = min(sample_size_pop3)
min.n.pop4 = min(sample_size_pop4)
min.n.total= min(sample_size_total)

#########################################
# TARGET SUMMARY STATISTICS
#########################################
# Calculate summary statistics for observed data
write.table( cbind( "mean_H_pop1",     "var_H_pop1",
                    "mean_H_pop2",     "var_H_pop2",
                    "mean_H_pop3",     "var_H_pop3",
                    "mean_H_pop4",     "var_H_pop4",
                    "mean_H_total",    "var_H_total",
                    "mean_He_pop1",    "var_He_pop1",
                    "mean_He_pop2",    "var_He_pop2",
                    "mean_He_pop3",    "var_He_pop3",
                    "mean_He_pop4",    "var_He_pop4",
                    "mean_He_total",   "var_He_total",
                    "mean_A_pop1",     "var_A_pop1",
                    "mean_A_pop2",     "var_A_pop2",
                    "mean_A_pop3",     "var_A_pop3",
                    "mean_A_pop4",     "var_A_pop4",
                    "mean_A_total",    "var_A_total",
                    "mean_Ar_pop1",    "var_Ar_pop1",
                    "mean_Ar_pop2",    "var_Ar_pop2",
                    "mean_Ar_pop3",    "var_Ar_pop3",
                    "mean_Ar_pop4",    "var_Ar_pop4",
                    "mean_Ar_total",   "var_Ar_total",
                    "mean_V_pop1",     "var_V_pop1",
                    "mean_V_pop2",     "var_V_pop2",
                    "mean_V_pop3",     "var_V_pop3",
                    "mean_V_pop4",     "var_V_pop4",
                    "mean_V_total",    "var_V_total",
                    "mean_R_pop1",     "var_R_pop1",
                    "mean_R_pop2",     "var_R_pop2",
                    "mean_R_pop3",     "var_R_pop3",
                    "mean_R_pop4",     "var_R_pop4",
                    "mean_R_total",    "var_R_total",
                    "mean_GW_pop1",    "var_GW_pop1",
                    "mean_GW_pop2",    "var_GW_pop2",
                    "mean_GW_pop3",    "var_GW_pop3",
                    "mean_GW_pop4",    "var_GW_pop4",
                    "mean_GW_total",   "var_GW_total",
                    "mean_P_pop1",     "var_P_pop1", 
                    "mean_P_pop2",     "var_P_pop2", 
                    "mean_P_pop3",     "var_P_pop3", 
                    "mean_P_pop4",     "var_P_pop4", 
                    "mean_P_total",    "var_P_total", 
                    "mean_HV_pop1",    "var_HV_pop1",
                    "mean_HV_pop2",    "var_HV_pop2",
                    "mean_HV_pop3",    "var_HV_pop3",
                    "mean_HV_pop4",    "var_HV_pop4",
                    "mean_HV_total",   "var_HV_total",
                    "mean_GST1_2",     "var_GST1_2",
                    "mean_GST1_3",     "var_GST1_3",
                    "mean_GST2_3",     "var_GST2_3",
                    "mean_GST1_4",     "var_GST1_4",
                    "mean_GST2_4",     "var_GST2_4",
                    "mean_GST3_4",     "var_GST3_4",
                    "mean_GST",        "var_GST",
                    "mean_deltamu1_2", "var_deltamu1_2",
                    "mean_deltamu1_3", "var_deltamu1_3",
                    "mean_deltamu2_3", "var_deltamu2_3",
                    "mean_deltamu1_4", "var_deltamu1_4",
                    "mean_deltamu2_4", "var_deltamu2_4",
                    "mean_deltamu3_4", "var_deltamu3_4"
),
file="target_sumstats.txt",sep=" ",
quote=F,col.names=F,row.names=F,append=F)

H_pop1   <- array(NA,number_of_loci)
H_pop2   <- array(NA,number_of_loci)
H_pop3   <- array(NA,number_of_loci)
H_pop4   <- array(NA,number_of_loci)
H_total  <- array(NA,number_of_loci)
He_pop1  <- array(NA,number_of_loci)
He_pop2  <- array(NA,number_of_loci)
He_pop3  <- array(NA,number_of_loci)
He_pop4  <- array(NA,number_of_loci)
He_total <- array(NA,number_of_loci)
A_pop1   <- array(NA,number_of_loci)
A_pop2   <- array(NA,number_of_loci)
A_pop3   <- array(NA,number_of_loci)
A_pop4   <- array(NA,number_of_loci)
A_total  <- array(NA,number_of_loci)
Ar_pop1   <- array(NA,number_of_loci)
Ar_pop2   <- array(NA,number_of_loci)
Ar_pop3   <- array(NA,number_of_loci)
Ar_pop4   <- array(NA,number_of_loci)
Ar_total  <- array(NA,number_of_loci)
V_pop1    <- array(NA,number_of_loci)
V_pop2    <- array(NA,number_of_loci)
V_pop3    <- array(NA,number_of_loci)
V_pop4    <- array(NA,number_of_loci)
V_total   <- array(NA,number_of_loci)
R_pop1    <- array(NA,number_of_loci)
R_pop2    <- array(NA,number_of_loci)
R_pop3    <- array(NA,number_of_loci)
R_pop4    <- array(NA,number_of_loci)
R_total   <- array(NA,number_of_loci)
GW_pop1   <- array(NA,number_of_loci)
GW_pop2   <- array(NA,number_of_loci)
GW_pop3   <- array(NA,number_of_loci)
GW_pop4   <- array(NA,number_of_loci)
GW_total  <- array(NA,number_of_loci)
P_pop1    <- array(F ,number_of_loci)
P_pop2    <- array(F ,number_of_loci)
P_pop3    <- array(F ,number_of_loci)
P_pop4    <- array(F ,number_of_loci)
P_total   <- array(F ,number_of_loci)
HV_pop1   <- array(NA,number_of_loci)
HV_pop2   <- array(NA,number_of_loci)
HV_pop3   <-array(NA,number_of_loci)
HV_pop4   <-array(NA,number_of_loci)
HV_total  <- array(NA,number_of_loci)
GST1_2    <-array(NA,number_of_loci)
GST1_3    <-array(NA,number_of_loci)
GST2_3    <-array(NA,number_of_loci)
GST1_4    <-array(NA,number_of_loci)
GST2_4    <-array(NA,number_of_loci)
GST3_4    <-array(NA,number_of_loci)
GST       <- array(NA,number_of_loci)
deltamu1_2  <- array(NA,number_of_loci)
deltamu1_3  <- array(NA,number_of_loci)
deltamu2_3  <- array(NA,number_of_loci)
deltamu1_4  <- array(NA,number_of_loci)
deltamu2_4  <- array(NA,number_of_loci)
deltamu3_4  <- array(NA,number_of_loci)

#declare all necessary function
He <- function(x)
{
    He <-1-sum( (table(x)/sum(table(x), na.rm=T))^2)
    return(He)
}
A <- function(x)
{
	A <- length(levels(x))
	return(A)
}
Ar <- function(x,y,z)
{
	Ar <- sum(1-choose(x-table(y),z)/choose(x,z))
}
V <- function(x) 
{
	V <- var(x)
	return(V)
}    
R <- function(x) #y = sample_size
{
	R <- max(x)-min(x)
	return(R)
}    
GW <- function(x) #y = sample_size
{
	GW <-A(x)/(R(x) +1)
	return(GW)
}    
HV <- function(x) #y = sample_size
{
	HV <- H(x)/V(x)
	return(HV)
} 
Gst <- function(x,y,z)
{
	Gst <-(H(z) - mean(cbind(H(x), H(y) ))) /H(z)
	return(Gst)
}
#TO DO: Ajouter Dest, Fst', Fis, etc...

n_pop1=array(NA,number_of_loci)
n_pop2=array(NA,number_of_loci)
n_pop3=array(NA,number_of_loci)
n_pop4=array(NA,number_of_loci)
Hs_pop1=array(NA,number_of_loci)
Hs_pop2=array(NA,number_of_loci)
Hs_pop3=array(NA,number_of_loci)
Hs_pop4=array(NA,number_of_loci)

for (locus in 1:number_of_loci){
  H_pop1[locus]  <- H(as.factor(data_pop1[,locus])) # 
  H_pop2[locus]  <- H(as.factor(data_pop2[,locus])) #
  H_pop3[locus]  <- H(as.factor(data_pop3[,locus])) #
  H_pop4[locus]  <- H(as.factor(data_pop4[,locus])) #
  H_total[locus] <- H(as.factor(data_total[,locus]))

  He_pop1[locus]  <- He(as.factor(data_pop1[,locus])) # 
  He_pop2[locus]  <- He(as.factor(data_pop2[,locus])) # 
  He_pop3[locus]  <- He(as.factor(data_pop3[,locus])) # 
  He_pop4[locus]  <- He(as.factor(data_pop4[,locus])) # 

  n_pop1[locus] <- length(which(!is.na(data_pop1[,locus])))
  Hs_pop1[locus] <- 1 - sum((table(as.factor(data_pop1[,locus]))/sum(table(as.factor(data_pop1[,locus])),na.rm=T))^2) - H_pop1[locus]/2/n_pop1[locus]
  Hs_pop1[locus] <- n_pop1[locus]/(n_pop1[locus]-1)*Hs_pop1[locus]
  n_pop2[locus] <- length(which(!is.na(data_pop2[,locus])))
  Hs_pop2[locus] <- 1 - sum((table(as.factor(data_pop2[,locus]))/sum(table(as.factor(data_pop2[,locus])),na.rm=T))^2) - H_pop2[locus]/2/n_pop2[locus]
  Hs_pop2[locus] <- n_pop2[locus]/(n_pop2[locus]-1)*Hs_pop2[locus] 
  n_pop3[locus] <- length(which(!is.na(data_pop3[,locus])))
  Hs_pop3[locus] <- 1 - sum((table(as.factor(data_pop3[,locus]))/sum(table(as.factor(data_pop3[,locus])),na.rm=T))^2) - H_pop3[locus]/2/n_pop3[locus]
  Hs_pop3[locus] <- n_pop3[locus]/(n_pop3[locus]-1)*Hs_pop3[locus]
  n_pop4[locus] <- length(which(!is.na(data_pop4[,locus])))
  Hs_pop4[locus] <- 1 - sum((table(as.factor(data_pop4[,locus]))/sum(table(as.factor(data_pop4[,locus])),na.rm=T))^2) - H_pop4[locus]/2/n_pop4[locus]
  Hs_pop4[locus] <- n_pop4[locus]/(n_pop4[locus]-1)*Hs_pop4[locus]

  A_pop1[locus]  <- A(as.factor(data_pop1[,locus]))
  A_pop2[locus]  <- A(as.factor(data_pop2[,locus]))
  A_pop3[locus]  <- A(as.factor(data_pop3[,locus]))
  A_pop4[locus]  <- A(as.factor(data_pop4[,locus]))
  A_total[locus]  <- A(as.factor(data_total[,locus]))

  Ar_pop1[locus]  <- Ar(sample_size_pop1[locus],as.factor(data_pop1[,locus]),min.n.pop1)
  Ar_pop2[locus]  <- Ar(sample_size_pop2[locus],as.factor(data_pop2[,locus]),min.n.pop2)
  Ar_pop3[locus]  <- Ar(sample_size_pop3[locus],as.factor(data_pop3[,locus]),min.n.pop3)
  Ar_pop4[locus]  <- Ar(sample_size_pop4[locus],as.factor(data_pop4[,locus]),min.n.pop4)

  V_pop1[locus]  <- var(data_pop1[,locus]) #
  V_pop2[locus]  <- var(data_pop2[,locus]) #
  V_pop3[locus]  <- var(data_pop3[,locus]) #
  V_pop4[locus]  <- var(data_pop4[,locus]) #
  V_total[locus] <- var(data_total[,locus]) #

  R_pop1[locus]  <- max(data_pop1[,locus])-min(data_pop1[,locus]) #
  R_pop2[locus]  <- max(data_pop2[,locus])-min(data_pop2[,locus]) #
  R_pop3[locus]  <- max(data_pop3[,locus])-min(data_pop3[,locus]) #
  R_pop4[locus]  <- max(data_pop4[,locus])-min(data_pop4[,locus]) #
  R_total[locus] <- max(data_total[,locus])-min(data_total[,locus]) #

  GW_pop1[locus]  <- A_pop1[locus]/(R_pop1[locus]+1) #
  GW_pop2[locus]  <- A_pop2[locus]/(R_pop2[locus]+1) #
  GW_pop3[locus]  <- A_pop3[locus]/(R_pop3[locus]+1) #
  GW_pop4[locus]  <- A_pop4[locus]/(R_pop4[locus]+1) #

  GW_total[locus] <- A_total[locus]/(R_total[locus]+1) #
  if (A_pop1[locus]>1)  P_pop1[locus]  <- T #
  if (A_pop2[locus]>1)  P_pop2[locus]  <- T #
  if (A_pop3[locus]>1)  P_pop3[locus]  <- T #
  if (A_pop4[locus]>1)  P_pop4[locus]  <- T #
  if (A_total[locus]>1) P_total[locus] <- T #
  HV_pop1[locus]  <- H_pop1[locus]/V_pop1[locus] # 
  HV_pop2[locus]  <- H_pop2[locus]/V_pop2[locus] # 
  HV_pop3[locus]  <- H_pop3[locus]/V_pop3[locus] # 
  HV_pop4[locus]  <- H_pop4[locus]/V_pop4[locus] # 
  HV_total[locus] <- H_total[locus]/V_total[locus] #

  GST1_2[locus] <- (H_total[locus] - mean(cbind(H_pop1[locus],H_pop2[locus])))/H_total[locus]
  GST1_3[locus] <- (H_total[locus] - mean(cbind(H_pop1[locus],H_pop3[locus])))/H_total[locus]
  GST2_3[locus] <- (H_total[locus] - mean(cbind(H_pop2[locus],H_pop3[locus])))/H_total[locus]
  GST1_4[locus] <- (H_total[locus] - mean(cbind(H_pop1[locus],H_pop4[locus])))/H_total[locus]
  GST2_4[locus] <- (H_total[locus] - mean(cbind(H_pop2[locus],H_pop4[locus])))/H_total[locus]
  GST3_4[locus] <- (H_total[locus] - mean(cbind(H_pop3[locus],H_pop4[locus])))/H_total[locus]
  GST[locus]    <- (H_total[locus] - mean(cbind(H_pop1[locus],H_pop2[locus],H_pop3[locus], H_pop4[locus])))/H_total[locus]

  deltamu1_2[locus] <- (mean(data_pop1[,locus])-mean(data_pop2[,locus]))^2
  deltamu1_3[locus] <- (mean(data_pop1[,locus])-mean(data_pop3[,locus]))^2
  deltamu2_3[locus] <- (mean(data_pop2[,locus])-mean(data_pop3[,locus]))^2
  deltamu1_4[locus] <- (mean(data_pop1[,locus])-mean(data_pop4[,locus]))^2
  deltamu2_4[locus] <- (mean(data_pop2[,locus])-mean(data_pop4[,locus]))^2
  deltamu3_4[locus] <- (mean(data_pop3[,locus])-mean(data_pop4[,locus]))^2
}

for (stat in c("H","He","A","Ar","R","GW","A","V","HV","P")){   #"A","V","HV","P"
  for (subsample in c("pop1","pop2","pop3","pop4","total")){
    assign(paste("mean",stat,subsample,sep="_"),mean(get(paste(stat,subsample,sep="_")),na.rm=T))
    assign(paste( "var",stat,subsample,sep="_"), var(get(paste(stat,subsample,sep="_")),na.rm=T))
  }
}
mean_GST1_2      <- mean(      GST1_2, na.rm=T )
mean_GST1_3      <- mean(      GST1_3, na.rm=T )
mean_GST2_3      <- mean(      GST2_3, na.rm=T )
mean_GST1_4      <- mean(      GST1_4, na.rm=T )
mean_GST2_4      <- mean(      GST2_4, na.rm=T )
mean_GST3_4      <- mean(      GST3_4, na.rm=T )
var_GST1_2       <-  var(      GST1_2, na.rm=T )
var_GST1_3       <-  var(      GST1_3, na.rm=T )
var_GST2_3       <-  var(      GST2_3, na.rm=T )
var_GST1_4       <-  var(      GST1_4, na.rm=T )
var_GST2_4       <-  var(      GST2_4, na.rm=T )
var_GST3_4       <-  var(      GST3_4, na.rm=T )
mean_GST         <- mean(      GST,    na.rm=T )
var_GST          <-  var(      GST,    na.rm=T )
mean_deltamu1_2  <- mean( deltamu1_2,  na.rm=T )
mean_deltamu1_3  <- mean( deltamu1_3,  na.rm=T )
mean_deltamu2_3  <- mean( deltamu2_3,  na.rm=T )
mean_deltamu1_4  <- mean( deltamu1_4,  na.rm=T )
mean_deltamu2_4  <- mean( deltamu2_4,  na.rm=T )
mean_deltamu3_4  <- mean( deltamu3_4,  na.rm=T )
var_deltamu2_3   <-  var( deltamu2_3,  na.rm=T )
var_deltamu1_2   <-  var( deltamu1_2,  na.rm=T )
var_deltamu1_3   <-  var( deltamu1_3,  na.rm=T )
var_deltamu1_4   <-  var( deltamu1_4,  na.rm=T )
var_deltamu2_4   <-  var( deltamu2_4,  na.rm=T )
var_deltamu3_4   <-  var( deltamu3_4,  na.rm=T )

write.table( cbind( 
                    mean_H_pop1,     var_H_pop1,
                    mean_H_pop2,     var_H_pop2,
                    mean_H_pop3,     var_H_pop3,
                    mean_H_pop4,     var_H_pop4,
                    mean_H_total,    var_H_total,
                    mean_He_pop1,    var_He_pop1,
                    mean_He_pop2,    var_He_pop2,
                    mean_He_pop3,    var_He_pop3,
                    mean_He_pop4,    var_He_pop4,
                    mean_He_total,   var_He_total,
                    mean_A_pop1,     var_A_pop1,
                    mean_A_pop2,     var_A_pop2,
                    mean_A_pop3,     var_A_pop3,
                    mean_A_pop4,     var_A_pop4,
                    mean_A_total,    var_A_total,
                    mean_Ar_pop1,    var_Ar_pop1,
                    mean_Ar_pop2,    var_Ar_pop2,
                    mean_Ar_pop3,    var_Ar_pop3,
                    mean_Ar_pop4,    var_Ar_pop4,
                    mean_Ar_total,   var_Ar_total,
                    mean_V_pop1,     var_V_pop1,
                    mean_V_pop2,     var_V_pop2,
                    mean_V_pop3,     var_V_pop3,
                    mean_V_pop4,     var_V_pop4,
                    mean_V_total,    var_V_total,
                    mean_R_pop1,     var_R_pop1,
                    mean_R_pop2,     var_R_pop2,
                    mean_R_pop3,     var_R_pop3,
                    mean_R_pop4,     var_R_pop4,
                    mean_R_total,    var_R_total,
                    mean_GW_pop1,    var_GW_pop1,
                    mean_GW_pop2,    var_GW_pop2,
                    mean_GW_pop3,    var_GW_pop3,
                    mean_GW_pop4,    var_GW_pop4,
                    mean_GW_total,   var_GW_total,
                    mean_P_pop1,     var_P_pop1, 
                    mean_P_pop2,     var_P_pop2, 
                    mean_P_pop3,     var_P_pop3, 
                    mean_P_pop4,     var_P_pop4, 
                    mean_P_total,    var_P_total, 
                    mean_HV_pop1,    var_HV_pop1,
                    mean_HV_pop2,    var_HV_pop2,
                    mean_HV_pop3,    var_HV_pop3,
                    mean_HV_pop4,    var_HV_pop4,
                    mean_HV_total,   var_HV_total,
                    mean_GST1_2,     var_GST1_2,
                    mean_GST1_3,     var_GST1_3,
                    mean_GST2_3,     var_GST2_3,
                    mean_GST1_4,     var_GST1_4,
                    mean_GST2_4,     var_GST2_4,
                    mean_GST3_4,     var_GST3_4,
                    mean_GST,        var_GST,
                    mean_deltamu1_2, var_deltamu1_2,
                    mean_deltamu1_3, var_deltamu1_3,
                    mean_deltamu2_3, var_deltamu2_3,
                    mean_deltamu1_4, var_deltamu1_4,
                    mean_deltamu2_4, var_deltamu2_4,
                    mean_deltamu3_4, var_deltamu3_4
),
file="target_sumstats.txt",sep=" ",
quote=F,col.names=F,row.names=F,append=T)

cat("empirical statistics succesfully computed" )
###################################
#Prepare ref.table for simulation #
##################################
#table were the output of the simulation will be written

if (!add_simulations){
  write.table( cbind( #### PARAMETERS
    "model","theta1",  "theta2", "theta3", "theta4", "theta5","theta6", "thetaA",
    if(model_type=="im") { 
    cbind("mig12", "mig13", "mig14", "mig21", "mig31", "mig34","mig23", "mig32","mig24","mig42","mig34", "mig43")
    } else {},   
    "tau_4_3","tau_3_2", "tau_2_1",
    # SUMMARY STATISTICS                      
    "mean_H_pop1",     "var_H_pop1",
    "mean_H_pop2",     "var_H_pop2",
    "mean_H_pop3",     "var_H_pop3",
    "mean_H_pop4",     "var_H_pop4",
    "mean_H_total",    "var_H_total",
    "mean_He_pop1",    "var_He_pop1",
    "mean_He_pop2",    "var_He_pop2",
    "mean_He_pop3",    "var_He_pop3",
    "mean_He_pop4",    "var_He_pop4",
    "mean_He_total",   "var_He_total",
    "mean_A_pop1",     "var_A_pop1",
    "mean_A_pop2",     "var_A_pop2",
    "mean_A_pop3",     "var_A_pop3",
    "mean_A_pop4",     "var_A_pop4",
    "mean_A_total",    "var_A_total",
    "mean_Ar_pop1",    "var_Ar_pop1",
    "mean_Ar_pop2",    "var_Ar_pop2",
    "mean_Ar_pop3",    "var_Ar_pop3",
    "mean_Ar_pop4",    "var_Ar_pop4",
    "mean_Ar_total",   "var_Ar_total",
    "mean_V_pop1",     "var_V_pop1",
    "mean_V_pop2",     "var_V_pop2",
    "mean_V_pop3",     "var_V_pop3",
    "mean_V_pop4",     "var_V_pop4",
    "mean_V_total",    "var_V_total",
    "mean_R_pop1",     "var_R_pop1",
    "mean_R_pop2",     "var_R_pop2",
    "mean_R_pop3",     "var_R_pop3",
    "mean_R_pop4",     "var_R_pop4",
    "mean_R_total",    "var_R_total",
    "mean_GW_pop1",    "var_GW_pop1",
    "mean_GW_pop2",    "var_GW_pop2",
    "mean_GW_pop3",    "var_GW_pop3",
    "mean_GW_pop4",    "var_GW_pop4",
    "mean_GW_total",   "var_GW_total",
    "mean_P_pop1",     "var_P_pop1", 
    "mean_P_pop2",     "var_P_pop2", 
    "mean_P_pop3",     "var_P_pop3", 
    "mean_P_pop4",     "var_P_pop4", 
    "mean_P_total",    "var_P_total", 
    "mean_HV_pop1",    "var_HV_pop1",
    "mean_HV_pop2",    "var_HV_pop2",
    "mean_HV_pop3",    "var_HV_pop3",
    "mean_HV_pop4",    "var_HV_pop4",
    "mean_HV_total",   "var_HV_total",
    "mean_GST1_2",     "var_GST1_2",
    "mean_GST1_3",     "var_GST1_3",
    "mean_GST2_3",     "var_GST2_3",
    "mean_GST1_4",     "var_GST1_4",
    "mean_GST2_4",     "var_GST2_4",
    "mean_GST3_4",     "var_GST3_4",
    "mean_GST",        "var_GST",
    "mean_deltamu1_2", "var_deltamu1_2",
    "mean_deltamu1_3", "var_deltamu1_3",
    "mean_deltamu2_3", "var_deltamu2_3",
    "mean_deltamu1_4", "var_deltamu1_4",
    "mean_deltamu2_4", "var_deltamu2_4",
    "mean_deltamu3_4", "var_deltamu3_4"
  ),
  file=reftable_file,sep=" ", quote=F, col.names=F,row.names=F, append=F)
} #END if(!add_simulations)

#########################################
# SIMULATION
#########################################
#SI: ms [total sample size] 1 -t [theta ref] -I 3 [sample size pop 1] [sample size pop 2] [sample size pop 3] [sample size pop 4] -n 1 [theta3] -n 2 [theta2] -n 4 [theta4] -ej [tau_3_2] 3 2 -eN [tau_3_2] [theta 4] -ej [tau_2_1] -eN [tau_2_1] [theta A]
#IM: ms [total sample size] 1 -t [theta ref] -I 3 [sample size pop 1] [sample size pop 2] [sample size pop 3] [sample size pop 4] 0 -m 1 2 [mig12] -m 2 1 [mig21] -n 1 [theta 1] -n 2 [theta 2] -ej [tau] 2 1 -eN [tau] [theta A]

for(i in 1:nsim){
   
  thetaRef=0.05 #hypothèse: Ne=1000 mu=1.25*10^-5
  theta1_min <- 0 
  theta1_max <- 30 

  if (model_type=="im") {
  model <- 1
  } else {
  model <- 0
  }

  theta1=theta2=theta3=theta4=theta5=theta6=thetaA=NULL 
  # take parameter values from priors
  theta1 <- runif(1, theta1_min, theta1_max) 
  theta2 <- runif(1, theta1_min, theta1_max)
  theta3 <- runif(1, theta1_min, theta1_max) 
  theta4 <- runif(1, theta1_min, theta1_max) 
  theta5 <- runif(1, theta1_min, theta1_max) 
  theta6 <- runif(1, theta1_min, theta1_max) 
  thetaA <- runif(1, theta1_min, theta1_max)

  if (mig_model=="full") {
    mig12   <- runif( 1, min=0, max=20 )
    mig13   <- runif( 1, min=0, max=20 )
    mig14   <- runif( 1, min=0, max=20 )
    mig21   <- runif( 1, min=0, max=20 )
    mig23   <- runif( 1, min=0, max=20 )
    mig24   <- runif( 1, min=0, max=20 )
    mig31   <- runif( 1, min=0, max=20 )
    mig32   <- runif( 1, min=0, max=20 )
    mig34   <- runif( 1, min=0, max=20 )
    mig41   <- runif( 1, min=0, max=20 )
    mig42   <- runif( 1, min=0, max=20 )
    mig43   <- runif( 1, min=0, max=20 )
   } else if (mig_model=="234_432") {
    mig12   <- 0 #runif( 1, min=0, max=20 )
    mig13   <- 0 #runif( 1, min=0, max=20 )
    mig14   <- 0 #runif( 1, min=0, max=20 )
    mig21   <- 0 #runif( 1, min=0, max=20 )
    mig23   <- runif( 1, min=0, max=20 )
    mig24   <- runif( 1, min=0, max=20 )
    mig31   <- 0 #runif( 1, min=0, max=20 )
    mig32   <- runif( 1, min=0, max=20 )
    mig34   <- runif( 1, min=0, max=20 )
    mig41   <- 0 #runif( 1, min=0, max=20 )
    mig42   <- runif( 1, min=0, max=20 )
    mig43   <- runif( 1, min=0, max=20 )
   } else if (mig_model=="23_32_34_43") {
    mig12   <- 0 #runif( 1, min=0, max=20 )
    mig13   <- 0 #runif( 1, min=0, max=20 )
    mig14   <- 0 #runif( 1, min=0, max=20 )
    mig21   <- 0 #runif( 1, min=0, max=20 )
    mig23   <- runif( 1, min=0, max=20 )
    mig24   <- 0 #runif( 1, min=0, max=20 )
    mig31   <- 0 #runif( 1, min=0, max=20 )
    mig32   <- runif( 1, min=0, max=20 )
    mig34   <- runif( 1, min=0, max=20 )
    mig41   <- 0 #runif( 1, min=0, max=20 )
    mig42   <- 0 #runif( 1, min=0, max=20 )
    mig43   <- runif( 1, min=0, max=20 )
  } else if (mig_model=="123_321") {
    mig12   <- runif( 1, min=0, max=20 )
    mig13   <- runif( 1, min=0, max=20 )
    mig14   <- 0 #runif( 1, min=0, max=20 )
    mig21   <- runif( 1, min=0, max=20 )
    mig23   <- runif( 1, min=0, max=20 )
    mig24   <- 0 #runif( 1, min=0, max=20 )
    mig31   <- runif( 1, min=0, max=20 )
    mig32   <- runif( 1, min=0, max=20 )
    mig34   <- 0 #runif( 1, min=0, max=20 )
    mig41   <- 0 #runif( 1, min=0, max=20 )
    mig42   <- 0 #runif( 1, min=0, max=20 )
    mig43   <- 0 #runif( 1, min=0, max=20 )
   } else if (mig_model=="12_21_32_23") {
    mig12   <- runif( 1, min=0, max=20 )
    mig13   <- 0 #runif( 1, min=0, max=20 )
    mig14   <- 0 #runif( 1, min=0, max=20 )
    mig21   <- runif( 1, min=0, max=20 )
    mig23   <- runif( 1, min=0, max=20 )
    mig24   <- 0 #runif( 1, min=0, max=20 )
    mig31   <- 0 #runif( 1, min=0, max=20 )
    mig32   <- runif( 1, min=0, max=20 )
    mig34   <- 0 #runif( 1, min=0, max=20 )
    mig41   <- 0 #runif( 1, min=0, max=20 )
    mig42   <- 0 #runif( 1, min=0, max=20 )
    mig43   <- 0 #runif( 1, min=0, max=20 )

#release your creativity here to add other migration models
 
  } else if (mig_model=="23_32") {
   mig12   <- mig21 <- mig13 <- mig31 <- 0 #runif( 1, min=0, max=20 )
   mig23   <- runif( 1, min=0, max=20 )
   mig32   <- runif( 1, min=0, max=20 )
   mig11   <- mig22   <- mig33   <- "X"
  } else if (mig_model=="12_21") {
    mig23   <- mig32 <- mig13 <- mig31 <- 0 #runif( 1, min=0, max=20 )
    mig12   <- runif( 1, min=0, max=20 )
    mig21   <- runif( 1, min=0, max=20 )
    mig11   <- mig22   <- mig33   <- "X"
  } else if (mig_model=="13_31") {
    mig12   <- mig21 <- mig23 <- mig32 <- 0  #runif( 1, min=0, max=20 )
    mig13   <- runif( 1, min=0, max=20 )
    mig31   <- runif( 1, min=0, max=20 )
    mig11   <- mig22   <- mig33   <- "X"
  } else if (mig_model=="23") {
    mig12   <- mig21 <- mig13 <- mig31 <- 0 #runif( 1, min=0, max=20 )
    mig23   <- runif( 1, min=0, max=20 )
    mig32   <- 0 #runif( 1, min=0, max=20 )
    mig11   <- mig22   <- mig33   <- "X"
  } else if (mig_model=="32") {
    mig12   <- mig21 <- mig13 <- mig31 <- 0 #runif( 1, min=0, max=20 )
    mig23   <- 0 #runif( 1, min=0, max=20 )
    mig32   <- runif( 1, min=0, max=20 )
    mig11   <- mig22   <- mig33   <- "X"
  } else if (mig_model=="12") {
    mig32   <- mig21 <- mig13 <- mig31 <- 0  #runif( 1, min=0, max=20 )
    mig23   <- 0 #runif( 1, min=0, max=20 )
    mig12   <- runif( 1, min=0, max=20 )
    mig11   <- mig22   <- mig33   <- "X"
  } else if (mig_model=="21") {
    mig12   <- mig32 <- mig13 <- mig31 <- 0 #runif( 1, min=0, max=20 )
    mig23   <- 0 #runif( 1, min=0, max=20 )
    mig21   <- runif( 1, min=0, max=20 )
    mig11   <- mig22   <- mig33   <- "X"
  } else if (mig_model=="13") {
    mig12   <- mig21 <- mig32 <- mig31 <- 0 #runif( 1, min=0, max=20 )
    mig23   <- 0 #runif( 1, min=0, max=20 )
    mig13   <- runif( 1, min=0, max=20 )
    mig11   <- mig22   <- mig33   <- "X"
  } else if (mig_model=="31") {
    mig12   <- mig21 <- mig13 <- mig32 <- 0 #runif( 1, min=0, max=20 )
    mig23   <- 0 #runif( 1, min=0, max=20 )
    mig31   <- runif( 1, min=0, max=20 )
    mig11   <- mig22   <- mig33   <- "X"
 } else if (mig_model=="no_mig") {
         mig21<-mig12<-mig13<-mig23<-mig31<-mig32<-0
 }
  #timing  
  tmin=0
  tmax=25
  tau_4_3 = runif(1, tmin, tmax)
  tau_3_2 = runif(1, tau_4_3, tmax)
  tau_2_1 = runif(1, tau_3_2, tmax) 
  
  #alpha parameter
  alpha <- runif( 1, min=0.5, max=1 )
  
  H_pop1   <- array(NA,number_of_loci)
  H_pop2   <- array(NA,number_of_loci)
  H_pop3   <- array(NA,number_of_loci)
  H_pop4   <- array(NA,number_of_loci)
  H_total  <- array(NA,number_of_loci)
  He_pop1  <- array(NA,number_of_loci)
  He_pop2  <- array(NA,number_of_loci)
  He_pop3  <- array(NA,number_of_loci)
  He_pop4  <- array(NA,number_of_loci)
  He_total <- array(NA,number_of_loci)
  A_pop1   <- array(NA,number_of_loci)
  A_pop2   <- array(NA,number_of_loci)
  A_pop3   <- array(NA,number_of_loci)
  A_pop4   <- array(NA,number_of_loci)
  A_total  <- array(NA,number_of_loci)
  Ar_pop1   <- array(NA,number_of_loci)
  Ar_pop2   <- array(NA,number_of_loci)
  Ar_pop3   <- array(NA,number_of_loci)
  Ar_pop4   <- array(NA,number_of_loci)
  Ar_total  <- array(NA,number_of_loci)
  V_pop1    <- array(NA,number_of_loci)
  V_pop2    <- array(NA,number_of_loci)
  V_pop3    <- array(NA,number_of_loci)
  V_pop4    <- array(NA,number_of_loci)
  V_total   <- array(NA,number_of_loci)
  R_pop1    <- array(NA,number_of_loci)
  R_pop2    <- array(NA,number_of_loci)
  R_pop3    <- array(NA,number_of_loci)
  R_pop4    <- array(NA,number_of_loci)
  R_total   <- array(NA,number_of_loci)
  GW_pop1   <- array(NA,number_of_loci)
  GW_pop2   <- array(NA,number_of_loci)
  GW_pop3   <- array(NA,number_of_loci)
  GW_pop4   <- array(NA,number_of_loci)
  GW_total  <- array(NA,number_of_loci)
  P_pop1    <- array(F ,number_of_loci)
  P_pop2    <- array(F ,number_of_loci)
  P_pop3    <- array(F ,number_of_loci)
  P_pop4    <- array(F ,number_of_loci)
  P_total   <- array(F ,number_of_loci)
  HV_pop1   <- array(NA,number_of_loci)
  HV_pop2   <- array(NA,number_of_loci)
  HV_pop3   <-array(NA,number_of_loci)
  HV_pop4   <-array(NA,number_of_loci)
  HV_total  <- array(NA,number_of_loci)
  GST1_2    <-array(NA,number_of_loci)
  GST1_3    <-array(NA,number_of_loci)
  GST2_3    <-array(NA,number_of_loci)
  GST1_4    <-array(NA,number_of_loci)
  GST2_4    <-array(NA,number_of_loci)
  GST3_4    <-array(NA,number_of_loci)
  GST       <- array(NA,number_of_loci)
  deltamu1_2  <- array(NA,number_of_loci)
  deltamu1_3  <- array(NA,number_of_loci)
  deltamu2_3  <- array(NA,number_of_loci)
  deltamu1_4  <- array(NA,number_of_loci)
  deltamu2_4  <- array(NA,number_of_loci)
  deltamu3_4  <- array(NA,number_of_loci)

   
  for (locus in 1:number_of_loci){

    ms_out_file <- paste( "msout_locus_", locus,"_batch_",batch,".txt", sep="" )
    if(model_type=="im") {
    ms_run <- paste(sample_size_total[locus], "1",      
                     "-t", thetaRef,                            
                     "-I 4", sample_size_pop1[locus], sample_size_pop2[locus],  sample_size_pop3[locus], sample_size_pop4[locus],   
                     "-ma X", mig12, mig13, mig14, mig21,"X", mig23, mig24, mig31, mig32 , "X", mig34, mig41, mig42, mig43, "X",     
                     "-n 1", theta1,                          
                     "-n 2", theta2,
                     "-n 3", theta3,
                     "-n 4", theta4,
                     "-n 5", theta5,
                     "-n 6", theta6,
                     "-ej" , tau_4_3, "4 3",
                     "-eN" , tau_4_3, theta5,                   
                     "-ej" , tau_3_2, "3 2",                         
                     "-eN" , tau_3_2, theta6,
                     "-ej" , tau_2_1, "2 1",
                     "-eN" , tau_2_1, thetaA,
                    #"-seeds","${JOB_ID}","${SGE_TASK_ID}","${SGE_TASK_ID}",
                     ">", ms_out_file)                          # output file
    } else { 
    #if(model_type=="si") {
      ms_run <- paste(sample_size_total[locus], "1",      
                     "-t", thetaRef,                              
                     "-I 4", sample_size_pop1[locus], sample_size_pop2[locus],  sample_size_pop3[locus], sample_size_pop4[locus],
                     #"-ma X", mig12, mig13, mig21,"X", mig23, mig31, mig32, "X",
                     "-n 1", theta1,
                     "-n 2", theta2,
                     "-n 3", theta3,
                     "-n 4", theta4,
                     "-n 5", theta5,
                     "-n 6", theta6,
                     "-ej" , tau_4_3, "4 3",
                     "-eN" , tau_4_3, theta5 ,                    
                     "-ej" , tau_3_2, "3 2",                         
                     "-eN" , tau_3_2, theta6,
                     "-ej", tau_2_1, "2 1",
                     "-eN", tau_2_1, thetaA,
                    #"-seeds","${JOB_ID}","${SGE_TASK_ID}","${SGE_TASK_ID}",
                     ">", ms_out_file)                          # output file
   }

    if(.Platform$OS.type == "unix") {
      #ms_run <- paste( "/home/ubuntu/R/ABC/msdir/./ms", ms_run, sep=" " ) 
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
    
    geno_pop1 <- genotypes[1:sample_size_pop1[locus]]
    geno_pop2 <- genotypes[(sample_size_pop1[locus]+1):(sample_size_pop1[locus]+sample_size_pop2[locus])]
    geno_pop3 <- genotypes[(sample_size_pop1[locus]+sample_size_pop2[locus]+1):(sample_size_pop1[locus]+sample_size_pop2[locus]+sample_size_pop3[locus])]
    geno_pop4 <- genotypes[(sample_size_pop1[locus]+sample_size_pop2[locus]+sample_size_pop3[locus]+1):sample_size_total[locus]]
    
     #calculating summary statistics
    H_pop1[locus]  <- H(as.factor(geno_pop1)) 
    H_pop2[locus]  <- H(as.factor(geno_pop2))
    H_pop3[locus]  <- H(as.factor(geno_pop3))
    H_pop4[locus]  <- H(as.factor(geno_pop4))
    H_total[locus] <- H(as.factor(genotypes))
    
    He_pop1[locus] <- He(geno_pop1)
    He_pop2[locus] <- He(geno_pop2)
    He_pop3[locus] <- He(geno_pop3)
    He_pop4[locus] <- He(geno_pop4)
    He_total[locus] <-He(genotypes)
    
    A_pop1[locus]  <- A(as.factor(geno_pop1))
    A_pop2[locus]  <- A(as.factor(geno_pop2))
    A_pop3[locus]  <- A(as.factor(geno_pop3))
    A_pop4[locus]  <- A(as.factor(geno_pop4))
    A_total[locus] <- A(as.factor(genotypes))

    Ar_pop1[locus] <-Ar(sample_size_pop1[locus],as.factor(geno_pop1),min.n.pop1)
    Ar_pop2[locus] <-Ar(sample_size_pop2[locus],as.factor(geno_pop2),min.n.pop2)
    Ar_pop3[locus] <-Ar(sample_size_pop3[locus],as.factor(geno_pop3),min.n.pop3)
    Ar_pop4[locus] <-Ar(sample_size_pop4[locus],as.factor(geno_pop4),min.n.pop4)
    Ar_total[locus]<-Ar(sample_size_total[locus],as.factor(genotypes),min.n.total)

    V_pop1[locus]  <- V(geno_pop1)
    V_pop2[locus]  <- V(geno_pop2)
    V_pop3[locus]  <- V(geno_pop3)
    V_pop4[locus]  <- V(geno_pop4)
    V_total[locus] <- V(genotypes)

    R_pop1[locus]  <- R(geno_pop1)
    R_pop2[locus]  <- R(geno_pop2)
    R_pop3[locus]  <- R(geno_pop3)
    R_pop4[locus]  <- R(geno_pop4)
    R_total[locus] <- R(geno_pop1)

    GW_pop1[locus]  <- A_pop1[locus]/(R_pop1[locus]+1)
    GW_pop2[locus]  <- A_pop2[locus]/(R_pop2[locus]+1)
    GW_pop3[locus]  <- A_pop3[locus]/(R_pop3[locus]+1)
    GW_pop4[locus]  <- A_pop4[locus]/(R_pop4[locus]+1)
    GW_total[locus] <- A_total[locus]/(R_total[locus]+1)

    if (A_pop1[locus]>1)  P_pop1[locus]  <- T
    if (A_pop2[locus]>1)  P_pop2[locus]  <- T
    if (A_pop3[locus]>1)  P_pop3[locus]  <- T
    if (A_pop4[locus]>1)  P_pop4[locus]  <- T
    if (A_total[locus]>1) P_total[locus] <- T
    HV_pop1[locus]  <- H_pop1[locus]/V_pop1[locus]
    HV_pop2[locus]  <- H_pop2[locus]/V_pop2[locus]
    HV_pop3[locus]  <- H_pop3[locus]/V_pop3[locus]
    HV_pop4[locus]  <- H_pop4[locus]/V_pop4[locus]
    HV_total[locus] <- H_total[locus]/V_total[locus]

    GST1_2[locus] <- Gst(as.factor(geno_pop1),as.factor(geno_pop2),as.factor(genotypes))
    GST1_3[locus] <- Gst(as.factor(geno_pop1),as.factor(geno_pop3),as.factor(genotypes))
    GST1_4[locus] <- Gst(as.factor(geno_pop1),as.factor(geno_pop4),as.factor(genotypes))
    GST2_3[locus] <- Gst(as.factor(geno_pop2),as.factor(geno_pop3),as.factor(genotypes))
    GST2_4[locus] <- Gst(as.factor(geno_pop2),as.factor(geno_pop4),as.factor(genotypes))
    GST3_4[locus] <- Gst(as.factor(geno_pop3),as.factor(geno_pop4),as.factor(genotypes))
    
    #GST1_2[locus] <- (H_total[locus] - mean(cbind(H_pop1[locus],H_pop2[locus])))/H_total[locus]
    #GST1_3[locus] <- (H_total[locus] - mean(cbind(H_pop1[locus],H_pop3[locus])))/H_total[locus]
    #GST2_3[locus] <- (H_total[locus] - mean(cbind(H_pop2[locus],H_pop2[locus])))/H_total[locus]
    #GST1_4[locus] <- (H_total[locus] - mean(cbind(H_pop1[locus],H_pop4[locus])))/H_total[locus]
    #GST2_4[locus] <- (H_total[locus] - mean(cbind(H_pop2[locus],H_pop4[locus])))/H_total[locus]
    #GST3_4[locus] <- (H_total[locus] - mean(cbind(H_pop3[locus],H_pop4[locus])))/H_total[locus]
    GST[locus] <- (H_total[locus] - mean(cbind(H_pop1[locus],H_pop2[locus],H_pop3[locus])))/H_total[locus]

    deltamu1_2[locus] <- (mean(geno_pop1) - mean(geno_pop2))^2
    deltamu1_3[locus] <- (mean(geno_pop1) - mean(geno_pop3))^2
    deltamu2_3[locus] <- (mean(geno_pop2) - mean(geno_pop3))^2
    deltamu1_4[locus] <- (mean(geno_pop1) - mean(geno_pop4))^2
    deltamu2_4[locus] <- (mean(geno_pop2) - mean(geno_pop4))^2 
    deltamu3_4[locus] <- (mean(geno_pop3) - mean(geno_pop4))^2 
  }
  
for (stat in c("H","He","A","Ar","R","GW","A","V","HV","P")){   #"A","V","HV","P"
  for (subsample in c("pop1","pop2","pop3","pop4", "total")){
    assign(paste("mean",stat,subsample,sep="_"),mean(get(paste(stat,subsample,sep="_")),na.rm=T))
    assign(paste( "var",stat,subsample,sep="_"), var(get(paste(stat,subsample,sep="_")),na.rm=T))
  }
}
 mean_GST1_2      <- mean(      GST1_2, na.rm=T )
 mean_GST1_3      <- mean(      GST1_3, na.rm=T )
 mean_GST2_3      <- mean(      GST2_3, na.rm=T )
 mean_GST1_4      <- mean(      GST1_4, na.rm=T )
 mean_GST2_4      <- mean(      GST2_4, na.rm=T )
 mean_GST3_4      <- mean(      GST3_4, na.rm=T )
 var_GST1_2       <-  var(      GST1_2, na.rm=T )
 var_GST1_3       <-  var(      GST1_3, na.rm=T )
 var_GST2_3       <-  var(      GST2_3, na.rm=T )
 var_GST1_4       <-  var(      GST1_4, na.rm=T )
 var_GST2_4       <-  var(      GST2_4, na.rm=T )
 var_GST3_4       <-  var(      GST3_4, na.rm=T )
 mean_GST         <- mean(      GST,    na.rm=T )
 var_GST          <-  var(      GST,    na.rm=T )
 mean_deltamu1_2  <- mean( deltamu1_2,  na.rm=T )
 mean_deltamu1_3  <- mean( deltamu1_3,  na.rm=T )
 mean_deltamu2_3  <- mean( deltamu2_3,  na.rm=T )
 mean_deltamu1_4  <- mean( deltamu1_4,  na.rm=T )
 mean_deltamu2_4  <- mean( deltamu2_4,  na.rm=T )
 mean_deltamu3_4  <- mean( deltamu3_4,  na.rm=T )
 var_deltamu2_3   <-  var( deltamu2_3,  na.rm=T )
 var_deltamu1_2   <-  var( deltamu1_2,  na.rm=T )
 var_deltamu1_3   <-  var( deltamu1_3,  na.rm=T )
 var_deltamu1_4   <-  var( deltamu1_4,  na.rm=T )
 var_deltamu2_4   <-  var( deltamu2_4,  na.rm=T )
 var_deltamu3_4   <-  var( deltamu3_4,  na.rm=T )
    
    write.table( cbind( ## PARAMETERS
    model,theta1, theta2, theta3, theta4, theta5, theta6, thetaA,  
    if(model_type=="im") { 
    cbind(mig12, mig13, mig14, mig21, mig31, mig34,mig23, mig32,mig24,mig42,mig34, mig43)
    } else {},   
    tau_4_3,tau_3_2, tau_2_1, 
    # SUMMARY STATISTICS
                    mean_H_pop1,     var_H_pop1,
                    mean_H_pop2,     var_H_pop2,
                    mean_H_pop3,     var_H_pop3,
                    mean_H_pop4,     var_H_pop4,
                    mean_H_total,    var_H_total,
                    mean_He_pop1,    var_He_pop1,
                    mean_He_pop2,    var_He_pop2,
                    mean_He_pop3,    var_He_pop3,
                    mean_He_pop4,    var_He_pop4,
                    mean_He_total,   var_He_total,
                    mean_A_pop1,     var_A_pop1,
                    mean_A_pop2,     var_A_pop2,
                    mean_A_pop3,     var_A_pop3,
                    mean_A_pop4,     var_A_pop4,
                    mean_A_total,    var_A_total,
                    mean_Ar_pop1,    var_Ar_pop1,
                    mean_Ar_pop2,    var_Ar_pop2,
                    mean_Ar_pop3,    var_Ar_pop3,
                    mean_Ar_pop4,    var_Ar_pop4,
                    mean_Ar_total,   var_Ar_total,
                    mean_V_pop1,     var_V_pop1,
                    mean_V_pop2,     var_V_pop2,
                    mean_V_pop3,     var_V_pop3,
                    mean_V_pop4,     var_V_pop4,
                    mean_V_total,    var_V_total,
                    mean_R_pop1,     var_R_pop1,
                    mean_R_pop2,     var_R_pop2,
                    mean_R_pop3,     var_R_pop3,
                    mean_R_pop4,     var_R_pop4,
                    mean_R_total,    var_R_total,
                    mean_GW_pop1,    var_GW_pop1,
                    mean_GW_pop2,    var_GW_pop2,
                    mean_GW_pop3,    var_GW_pop3,
                    mean_GW_pop4,    var_GW_pop4,
                    mean_GW_total,   var_GW_total,
                    mean_P_pop1,     var_P_pop1, 
                    mean_P_pop2,     var_P_pop2, 
                    mean_P_pop3,     var_P_pop3, 
                    mean_P_pop4,     var_P_pop4, 
                    mean_P_total,    var_P_total, 
                    mean_HV_pop1,    var_HV_pop1,
                    mean_HV_pop2,    var_HV_pop2,
                    mean_HV_pop3,    var_HV_pop3,
                    mean_HV_pop4,    var_HV_pop4,
                    mean_HV_total,   var_HV_total,
                    mean_GST1_2,     var_GST1_2,
                    mean_GST1_3,     var_GST1_3,
                    mean_GST2_3,     var_GST2_3,
                    mean_GST1_4,     var_GST1_4,
                    mean_GST2_4,     var_GST2_4,
                    mean_GST3_4,     var_GST3_4,
                    mean_GST,        var_GST,
                    mean_deltamu1_2, var_deltamu1_2,
                    mean_deltamu1_3, var_deltamu1_3,
                    mean_deltamu2_3, var_deltamu2_3,
                    mean_deltamu1_4, var_deltamu1_4,
                    mean_deltamu2_4, var_deltamu2_4,
                    mean_deltamu3_4, var_deltamu3_4
   ),
  file=reftable_file,sep=" ",
  quote=F,col.names=F,row.names=F,append=T)

}
