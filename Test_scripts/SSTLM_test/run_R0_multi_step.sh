#!/bin/bash



hpc_switch=1
date_time=$(date '+%d-%m-%Y %H:%M:%S')
###########__________Run script__________#############
    ################ Hpc machine ################
if [ "$hpc_switch" == 1 ]
 then

sim_name="-HPC"
module load python/3.6.5
module load python-libs/3.1.0

#$ -cwd -V
#$ -l h_rt=48:00:00
#$ -t 1-100

python3 R0-multi-step.py $SGE_TASK_ID $date_time $sim_name

###########__________Run script__________#############
    ################ My machine ################

elif [ "$hpc_switch" == 0 ]
 then
# Local machine

job_id=1
sim_name="-LOCAL-"

python3 R0-multi-step.py $job_id $date_time $sim_name

fi
echo "Simulations Finished"
