#!/00-scripts/bash
#PBS -A ihv-653-ab
#PBS -N gfit.si
#PBS -o gfit.si
#PBS -e gfit.si
#PBS -l walltime=03:30:00
#PBS -l nodes=1:ppn=8  -l mem=23gb
#PBS -M quentinrougemont@orange.fr
#PBS -m ea 

# Move to directory where job was submitted
cd $PBS_O_WORKDIR

#global variable
model="../si.simul.ABC.txt"
FOLDER=./02-results/ppc.si
NREPS=5000
POP1=../../01-data/pop5cov
POP2=../../01-data/pop6ber
NLOC=$( awk -F " " '{print NF }' $POP1  |head -1 )
repeat=../../01-data/repeat.motif

# Create folder and move into it
mkdir "$FOLDER"
cd "$FOLDER"

#load R
module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0

# Launch Rscript
Rscript ../../00-scripts/rscript/04.goodness.of.fit.R "$model" "$POP1" "$POP2" "$NREPS" "$NLOC" "$repeat"

wait
sleep 30
