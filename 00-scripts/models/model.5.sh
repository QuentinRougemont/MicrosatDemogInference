#!/bin/bash

# Global variables
FOLDER=$1
NREPS=37500
POP1=../../01-data/pop1aa
POP2=../../01-data/pop2aa
repeat=../../01-data/repeat.motif
# Create folder and move into it
mkdir "$FOLDER"
cd "$FOLDER"
NLOC=$( awk -F " " '{print NF }' $POP1  |head -1 )

#load R
module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0

# Launch Rscript
Rscript ../../00-scripts/rscript/Simul_equilibrium.R "$POP1" "$POP2" "$NREPS" "$NLOC" "$repeat"

