#!/bin/bash

###########__________Parameters and settings__________#############
# hpc_switch | choose between hpc or local usage
# hpc_mode | choose to find phase diagram or large individual probability distributions
# data | choose between : ["lattice", "rand_uk", "Qro"]
# hpc_mode | flips between local machine and HPC
# sim_name | input a string to append to the output file to identify simulation runs

hpc_switch=0
data_type="lattice" # ["lattice", "channel"]
sim_name="-vel-km-day-test"

###########__________Run script__________#############
if [ "$hpc_switch" == 1 ]
 then
################ Hpc machine ################
module load python/3.6.5
module load python-libs/3.1.0
date_time=$(date '+%d-%m-%Y %H:%M:%S')
#$ -cwd -V
#$ -l h_rt=48:00:00
#$ -t 1-100

python3 mkdir.py $date_time $data_type $sim_name
######### find epidemiological phase space diagram #########
python3 main_SSTLM_phase.py  $SGE_TASK_ID $date_time $data_type $sim_name

elif [ "$hpc_switch" == 0 ]
 then
# Local machine

job_id=25
date_time=$(date '+%d-%m-%Y %H:%M:%S')
python3 mkdir.py  $date_time $data_type $sim_name
######### find epidemiological phase space diagram #########
python3 main_SSTLM_phase.py $job_id $date_time $data_type $sim_name

fi
echo "Simulations Finished"
