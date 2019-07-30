#!/bin/bash

###########__________Parameters and settings__________#############
# hpc_switch | choose between hpc or local usage
# hpc_mode | choose to find phase diagram or large individual probability distributions
# data | choose between : ["lattice", "rand_uk", "Qro"]
# hpc_mode | flips between local machine and HPC
# sim_name | input a string to append to the output file to identify simulation runs

hpc_switch=0
sim_name="test"
L=3 # index of [50, 100, 150, 200, 250, 300]m
beta=40 #

###########__________Run script__________#############
if [ "$hpc_switch" == 1 ]
 then
################ Hpc machine ################
module load python/3.6.5
module load python-libs/3.1.0
date_time=$(date '+%d-%m-%Y %H:%M:%S')
#$ -cwd -V
#$ -l h_rt=24:00:00

# ######## Run simulation ######## #
python3 pde_model_main.py $sim_name $L $beta
elif [ "$hpc_switch" == 0 ]
 then

# ######## local machine ######### #
       #### Run simulation ####

python3 pde_model_main.py $sim_name $L $beta

fi
echo "Simulations Finished"
