#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N reshape
##PBS -o reshape.out
##PBS -e reshape.err
#PBS -l walltime=00:25:00
#PBS -l nodes=1:ppn=8
##PBS -M quentinrougemont@orange.fr
##PBS -m bea

# Move to directory where job was submitted
cd $PBS_O_WORKDIR

target=1000000 #set the number of wanted simulations

#ici faire test pour vérifier qu'on est bien le bon nombre de simulations totale
#1 on list les modèles
#list=$( ls -d 02-results/*/ | sed -e 's/\([0-9]*\)//g' -e 's/._\///g'  -e 's/02-results\///g' |uniq |less )
ls -d 02-results/*/ | sed -e 's/\([0-9]*\)//g' -e 's/._\///g'  -e 's/-results\///g' |uniq > list

for j in $(cat list) ; do mkdir 02-results/$j.glob ; done 
for j in  $(cat list) ; do mv 02-results/$j.* 02-results/$j.glob ; done 
for i in $(cat list) ; do for k in $(find 02-results/$i.glob -name reference_table.txt) ; do cat $k |grep -v model >> 02-results/$i.ABC.txt ; done ; done  
sed -i '/^$/d' 02-results/*.ABC.txt 
for i in $(head -n 1 list) ; do for k in $(find 02-results/$i.glob -name  target_sumstats.txt ); do cp $k 02-results/ ; done ; done 

for j in  *.glob ; do rm 02-results/$j/*/msout_locus*  02-results/$j/*/seedms ; done
for i in 02-results ; do for j  in  $i/*.glob  ; do  tar -zcvf $j.tar.gz $j ; done ;done
for j in $(ls -d 02-results/*glob/ ) ; do  rm  -r $j ; done
rm *.err *.out


for i in $(cat list) ; do wc -l 02-results/$i.ABC.txt  |awk '{print $1}' >>  stat.tmp ; done

sum_stat=$(awk '{total += $1} END {print total/NR}' stat.tmp )
if [[ $sum_stat != $target ]]  ; 
then 
	echo "you loose" 
	echo "error! Number of simulation in one of the model is not = "$target" !!"
	echo "please check models"
	exit
fi

rm list *.tmp 
  
