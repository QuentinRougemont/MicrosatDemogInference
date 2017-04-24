#!/bin/bash

#script to submit robustness of MOAB architecture using msub
for model in $(cat 01-data/model_name)  ; do
      export model
      msub -v model 00-scripts/01.robustess.sh
done

