#!/bin/bash

# Execute from termial to run the model, either to simulate indiviual realisations locally or ensemble simluations, either on the HPC or locally.

###########__________Parameters and settings__________#############
# hpc_switch | choose between hpc or local usage
# sim_name | input a string to append to the output file to identify simulation runs
# data_type | currently set to lattice i.e. simple square homogeneous lattice 

hpc_switch=1
data_type="lattice"

###########__________Run script__________#############
if [ "$hpc_switch" == 1 ]
 then
################ Hpc machine ################
module load python/3.6.5
module load python-libs/3.1.0
date_time=$(date '+%d-%m-%Y %H:%M:%S')

#$ -cwd -V
#$ -l h_rt=48:00:00
#$ -t 1-10

sim_name="HPC"
python3 mkdir.py $date_time $data_type $sim_name
SGE_TASK_ID=1

python3 main.py $SGE_TASK_ID $date_time $data_type $sim_name
elif [ "$hpc_switch" == 0 ]
 then

################ Local machine ###############
job_id=25
date_time=$(date '+%d-%m-%Y %H:%M:%S')
sim_name="LCL"
python3 mkdir.py  $date_time $data_type $sim_name
python3 main.py $job_id $date_time $data_type $sim_name

fi
echo "Simulations Finished"
