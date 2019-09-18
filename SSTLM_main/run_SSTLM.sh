#!/bin/bash

# Execute from termial to run the model, either to simulate indiviual realisations locally or ensemble simluations, either on the HPC or locally.

###########__________Parameters and settings__________#############
# hpc_switch | choose between hpc or local usage
# sim_name | input a string to append to the output file to identify simulation runs
# data_type | currently set to lattice i.e. simple square homogeneous lattice 

hpc_switch=1

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

mode="HPC"
sim_type="-high_res-"  # HPC --> two sim_types : ['-high_res-', '-full_param-']
sim_name="ell_25"

python3 mkdir.py $date_time $mode $sim_type $sim_name
python3 main.py $SGE_TASK_ID $date_time $data_type $mode $sim_type $sim_name
elif [ "$hpc_switch" == 0 ]
 then

################ Local machine ###############
job_id=25
date_time=$(date '+%d-%m-%Y %H:%M:%S')
mode="LCL"
sim_type="-anim"  # LCL --> two sim_types : ['-anim', '-ens']
sim_name="-channel_test"
python3 mkdir.py  $date_time $data_type $mode $sim_type $sim_name
python3 main.py $job_id $date_time $data_type $mode $sim_type $sim_name

fi
echo "Simulations Finished"
