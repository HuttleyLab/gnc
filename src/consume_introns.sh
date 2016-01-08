#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l walltime=05:00:00,mem=4096GB,ncpus=2048
#PBS -V
# usage: consume.py [-h] -i INPUT_DIRECTORY -o OUTPUT_COLLECTION [-b DB_HOST]
#                   [-L LOG_LEVEL] [-l LOG_NAME] [-c CODON_POSITION]
mpiexec -np $PBS_NCPUS python consume.py -i $SHORT/data/human_macaque_marmoset_introns/data/introns -o introns.data -L DEBUG -b $DB_HOST -a 1500 -t ../config/introns.nwk
