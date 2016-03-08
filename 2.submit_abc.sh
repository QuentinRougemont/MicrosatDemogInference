#!/bin/bash

PROGNAME=${0##*/}
PROGVERSION=20-10-2015

echo ""
echo "Script to submit Coal Sims (script version $PROGVERSION)"

## parameters and their defualt value
nbJobs=100     		              # number of replciates

#exeRelative=0                     # exe is localy and not installed
queue=unlimitq		                  # Must be change according to the 

isTest=0                          # if 1 no job is submitted, but scripts are written 

## get the arguments
SHORTOPTS="i:n:v"
LONGOPTS="files: input: queue: version files: name:"
ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "$PROGNAME" -- "$@")
eval set -- "$ARGS"
while  [ ! -z "$1" ]; do
    case "$1" in
       -v|--version) echo "version $PROGVERSION"; exit;;
       -i|--input)   input=$2; shift;;
       -n)           nbJobs=$2; shift ;;
       --files)      copyFiles=$2; shift;;
       --name)       simName=$2; shift;;
       --queue)      queue=$2; shift;;       
       --)           shift; break;;
       -*)           break;;
       *)            echo "Unknown parameter '$1'!" ;;
    esac
    shift
done

## are the arguments passed without option
if [ "$#" == 2 -a "$input" == "" -a "$nbJobs" == "" ]; then
	input=$1
	nbJobs=$2
elif [ "$#" != 0 ]; then
	echo "ERROR: the number of parameters is wrong!"
	exit 1
fi

## pre-tests, strip the input to directory and file name
if [ "$input" == "" ]; then
	echo "ERROR: ini file (parameter -i) is not specified!"
	exit 1
fi	

ini=$(basename $input)
folder=$(dirname  $input)
if [ "$folder" == "" ]; then
	echo "ERROR: ini file is not located in a subfolder (parameter -i)!"
	exit 1
fi	

if [ "$simName" == "" ]; then
  if [ "$folder" == "." ]; then
	  simName=$exe
  else
	  simName=${folder##*/}
	fi
fi
origDir=`pwd`

#****************************************************************
# function to launch the job on a node
function write_submitJobs()
{ 
   # write the script
   
   (
      ## specification #############################
      echo "#!/bin/bash"
      echo "#$ -N "$simName
      echo "#$ -q "$queue
      echo "#$ -cwd"
      echo "#$ -o $origDir/sge_logs/"$simName"_\$JOB_ID_\$TASK_ID.out"
	  echo "#$ -e $origDir/sge_logs/"$simName"_\$JOB_ID_\$TASK_ID.err"
      
      # get the computer name ####################################
      echo "echo "
      echo "echo ==== Node name: \`hostname\` ===="
      echo "echo "
      
      ## working directory #############################    
      home=`pwd`/                                         # absolut current path
        echo "work=job_\${JOB_ID}_\${SGE_TASK_ID}"         # name of the working directory
      echo "mkdir -p \$work"                              # make the temp folder
      echo "cd \$work"          
      
      ## copy all the stuff #############################
      for i in $copyFiles                               # get the other parameters
			do
			  echo cp $home$i .
			done
      
      ## launching the programs #############################
     
     echo "sed -i  "s/\${JOB_ID}/\${JOB_ID}/g" SimulPAN_parallel.R "
     echo "sed -i  "s/\${SGE_TASK_ID}/\${SGE_TASK_ID}/g" SimulPAN_parallel.R " 
     echo "sed -i  "s/\${SGE_TASK_ID}/\${SGE_TASK_ID}/g" SimulPAN_parallel.R " 
     echo "R CMD BATCH SimulPAN_parallel.R"
       
      # get the size #############################################
      echo "echo "
      echo "echo ==== Output size: \`du -sh\` ===="
      echo "echo "

      # bring back output ########################################
      	echo "mv ../\$work "$home"job_\${SGE_TASK_ID}"
      
   ) > "submitJobs.sh"
}
#****************************************************************
# "the main function"

# create the folders for the output
mkdir -p sge_logs

# write scipts...
cd $folder

write_submitJobs 
chmod +x "submitJobs.sh"

# launch the jobs
	qsub -t 1-$nbJobs:1 "submitJobs.sh"

cd ../

