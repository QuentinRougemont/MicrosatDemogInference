#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N robustess
#PBS -o robust.out
#PBS -e robust.err
#PBS -l walltime=36:00:00
#PBS -l nodes=1:ppn=8 -l mem=23gb
##PBS -M quentinrougemont@orange.fr
#PBS -m ea
##PBS -t [1-10]%10

# Move to directory where job was submitted
cd $PBS_O_WORKDIR

module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0
module load apps/gnu-parallel/20131022
#parallel

# Get parameters
target=sc # $1
if [[ -z "$target" ]]
then
    echo "Error: need model name (eg: AM.1)"
    exit
fi

# Launching jobs in parallel
echo running model: "$target"
seq 1 100 1000 > lower_bounds

cat lower_bounds | parallel -j 8 00-scripts/03.robustess.sh {} "$target"
