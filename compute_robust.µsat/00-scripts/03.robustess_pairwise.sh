#!/bin/bash
# Verify if 100 lines are done, if not run 03_robustness.R

# Global variables
lowbnd=$1
target=$2
target2=$3
model1="02-results/im.simul.ABC.txt"
model2="02-results/si.simul.ABC.txt"
model3="02-results/sc.simul.ABC.txt"
model4="02-results/am.simul.ABC.txt"

# Launch R script
if [[ ! -f 02.done/"$target"_"$target2"/"$lowbnd" ]]
then
        Rscript 00-scripts/rscript/02.robustess_pairwise.R "$lowbnd" "$target" "$target2" "$model1" "$model2" "$model3" "$model4" && touch 02.done/"$target"_"$target2"/"$lowbnd"
fi

