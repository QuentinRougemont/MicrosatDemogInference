#!/bin/bash
# Verify if 100 lines are done, if not run 03_robustness.R

# Global variables
lowbnd=$1
target=$2
model_type="si"
model1="01-data/si.simul.den.ABC.txt"
model2="01-data/si.simul.ber.ABC.txt"
model3="01-data/si.simul.cov.ABC.txt" 
# Launch R script
if [[ ! -f 02.done/"$target"/"$lowbnd" ]]
then
        Rscript 00-scripts/rscript/02-robustess.R "$lowbnd" "$target" "$model1" "$model2" "$model3" "$model_type" && touch 02.done/"$target"/"$lowbnd"
fi

