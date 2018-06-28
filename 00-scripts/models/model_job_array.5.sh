#!/bin/bash
#PBS -A ihv-653-aa
#PBS -N OutputTest
##PBS -o OutTest.out
##PBS -e OutTest.err
#PBS -l walltime=12:01:00
#PBS -l nodes=1:ppn=8
#PBS -M quentinrougemont@orange.fr
##PBS -m ea 
#PBS -t [1-15]%15

# Move to directory where job was submitted
cd "${PBS_O_WORKDIR}"

# Folder to run simulations
MODEL=./00-scripts/models/model.5.sh
FOLDER=./02-results/eq.simul.$MOAB_JOBARRAYINDEX

for i in $(seq 8)
do
    sleep 0 
    ./"$MODEL" "$FOLDER"_"$i"  &
done

#echo "sed -i  "s/\${JOB_ID}/\${JOB_ID}/g" SimulPAN_parallel.R "
#echo "sed -i  "s/\${SGE_TASK_ID}/\${SGE_TASK_ID}/g" SimulPAN_parallel.R " 
#echo "sed -i  "s/\${SGE_TASK_ID}/\${SGE_TASK_ID}/g" SimulPAN_parallel.R " 
#echo "R CMD BATCH SimulPAN_parallel.R"

# Wait for all simulations to finish
wait
sleep 30
