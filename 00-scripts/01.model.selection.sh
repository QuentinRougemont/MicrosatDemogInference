#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N ModelSelection1
#PBS -o ModSel1.out
#PBS -e ModSel1.err
#PBS -l walltime=02:30:00
#PBS -l nodes=1:ppn=8 -l mem=23gb
#PBS -M quentinrougemont@orange.fr
#PBS -m ea 
##PBS -t [1-10]%10

# Move to directory where job was submitted
cd $PBS_O_WORKDIR

module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0

cd ./02-results
Rscript ../00-scripts/rscript/01.model.choice.R im.simul.ABC.txt si.simul.ABC.txt sc.simul.ABC.txt  am.simul.ABC.txt

wait
sleep 30
