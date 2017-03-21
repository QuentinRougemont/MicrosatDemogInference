#!/00-scripts/bash
#PBS -A ihv-653-ab
#PBS -N param.sc.het.het
#PBS -o Paramsc.het.het.out
#PBS -e Paramsc.het.het.err
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=8  -l mem=23gb
#PBS -M quentinrougemont@orange.fr
#PBS -m ea 

# Move to directory where job was submitted
cd $PBS_O_WORKDIR

module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0


cd ./02-results
Rscript ../00-scripts/rscript/03.param.estim.sc.R sc.simul.ABC.txt


wait
sleep 30