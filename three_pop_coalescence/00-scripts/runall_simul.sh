#!/bin/bash
# WARNING! All scripts must been edited prior to submitting and not changed until run!
# HINT: Give meaningful job names in submission scripts to distinguish them

# Global variables
SCRIPTPATH="./00-scripts"

echo "    DO NOT USE WITHOUT APPROPRIATE ABC BACKGROUND"

p01=$(msub  "$SCRIPTPATH"/models/model_job_array.1.sh  | tail -n 1 | awk '{print $1}') #si
p02=$(msub  "$SCRIPTPATH"/models/model_job_array.1b.sh  | tail -n 1 | awk '{print $1}')  #im 
p03=$(msub  "$SCRIPTPATH"/models/model_job_array.1c.sh  | tail -n 1 | awk '{print $1}')  #im 
#p04=$(msub  "$SCRIPTPATH"/models/model_job_array.2.sh  | tail -n 1 | awk '{print $1}')  #im 
#p05=$(msub  "$SCRIPTPATH"/models/model_job_array.2.a.sh  | tail -n 1 | awk '{print $1}')  #im 
#p06=$(msub  "$SCRIPTPATH"/models/model_job_array.2.b.sh  | tail -n 1 | awk '{print $1}')  #im 
#p07=$(msub  "$SCRIPTPATH"/models/model_job_array.2.c.sh  | tail -n 1 | awk '{print $1}')  #im 

#p03=$(msub  "$SCRIPTPATH"/models/model_job_array.3.sh  | tail -n 1 | awk '{print $1}')  #sc 
#p04=$(msub  "$SCRIPTPATH"/models/model_job_array.4.sh  | tail -n 1 | awk '{print $1}')  #am 

00-scripts/04.paramestim.si.3pop.sh si.simul.ber.ABC.txt  pop6ber pop1den pop5cov
#reshape the simul and process to model choice
s01=$(msub -l depend=$p01:$p02:$p03:$p04:$p05:$p06:$p07 "$SCRIPTPATH"/00.reshape.sh | tail -n 1 | awk '{print $1}')
s02=$(msub -l depend=$s01 "$SCRIPTPATH"/01.model.selection3pops.si.sh   | tail -n 1 | awk '{print $1}')
s03=$(msub -l depend=$s01 "$SCRIPTPATH"/04.paramestim.si.3pop.sh si.simul.ber.ABC.txt  pop6ber pop1den pop5cov | tail -n 1 | awk '{print $1}')
s04=$(msub -l depend=$s01 "$SCRIPTPATH"/04.paramestim.si.3pop.sh si.simul.cov.ABC.txt  pop5cov pop6ber pop1den | tail -n 1 | awk '{print $1}')
s04=$(msub -l depend=$s01 "$SCRIPTPATH"/04.paramestim.si.3pop.sh si.simul.den.ABC.txt  pop1den pop6ber pop5cov | tail -n 1 | awk '{print $1}')
#s06=$(msub -l depend=$s04 "$SCRIPTPATH"/02.robustess.sh         | tail -n 1 | awk '{print $1}')
#s07=$(msub -l depend=$s04 "$SCRIPTPATH"/02.robustess_sc.sh      | tail -n 1 | awk '{print $1}')
#s08=$(msub -l depend=$s04 "$SCRIPTPATH"/02.robustess_im.sh      | tail -n 1 | awk '{print $1}')
#s09=$(msub -l depend=$s04 "$SCRIPTPATH"/02.robustess_am.sh      | tail -n 1 | awk '{print $1}')

#s10=$(msub -l depend=$s04 "$SCRIPTPATH"/05.goodness.of.fit.im.sh  | tail -n 1 | awk '{print $1}')
#s11=$(msub -l depend=$s05 "$SCRIPTPATH"/05.goodness.of.fit.sc.sh  | tail -n 1 | awk '{print $1}')

#confirm job submissions
echo "All jobs submitted"
#end of script
