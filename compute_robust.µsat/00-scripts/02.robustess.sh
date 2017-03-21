#!/bin/bash
#SBATCH -J "rob.ABC"
#SBATCH -o log_%j
#SBATCH -c 4
#SBATCH -p ibismax
#SBATCH -A ibismax
##SBATCH --mail-type=FAIL
##SBATCH --mail-user=YOUREMAIL
#SBATCH --time=1-00:00
#SBATCH --mem=32G

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Get parameters
target=$1
if [[ -z "$target" ]]
then
    echo "Error: need model name (eg: AM.1)"
    exit
fi

# Launching jobs in parallel
echo running model: "$target"
seq 1 100 1000 > lower_bounds

cat lower_bounds | parallel -j 4 00-scripts/03.robustess.sh {} "$target"
