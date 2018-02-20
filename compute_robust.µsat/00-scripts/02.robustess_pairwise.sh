#!/bin/bash
#SBATCH -J "ABCRobust"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p low-suspend
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=asdf
#SBATCH --time=10-00:00
#SBATCH --mem=200G

# Move to directory where job was submitted
#cd $PBS_O_WORKDIR

#module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0
#module load apps/gnu-parallel/20131022
#parallel

# Get parameters
target=$1 #si # $1
target2=$2 #sc
if [[ -z "$target" ]]
then
    echo "Error: need model name (eg: AM.1)"
    exit
fi

if [[ -z "$target2" ]]
then
    echo "Error: need model name (eg: AM.1)"
    exit
fi

# Launching jobs in parallel
echo running model: "$target" vs "$target2"
seq 1 100 1000 > lower_bounds
#cat lower_bounds | parallel -j 8 00-scripts/03.robustess.sh {} "$target"
cat lower_bounds | parallel -j 10 00-scripts/03.robustess_pairwise.sh {} "$target" "$target2"
