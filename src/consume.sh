#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l walltime=05:00:00,mem=4096GB,ncpus=2048
#PBS -V
DB_HOST=r2452
# usage: consume.py [-h] -i INPUT_DIRECTORY -o OUTPUT_COLLECTION [-b DB_HOST]
#                   [-L LOG_LEVEL] [-l LOG_NAME] [-c CODON_POSITION]
mpiexec -np $PBS_NCPUS python consume.py -i $SHORT/data/mammals -o mammals.data -L DEBUG -b $DB_HOST -a 1500 -t ../config/mammals.nwk
