#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l walltime=05:00:00,mem=1096GB,ncpus=512
#PBS -V
DB_HOST=r1822
# usage: consume.py [-h] -i INPUT_DIRECTORY -o OUTPUT_COLLECTION [-b DB_HOST]
#                   [-L LOG_LEVEL] [-l LOG_NAME] [-c CODON_POSITION]
SHORT=/short/xe9/caj248
mpiexec -np $PBS_NCPUS python consume.py -i $SHORT/../data/hum_xen_fug -o hum_xen_fug.data -L DEBUG -b ${DB_HOST} -a 1500 -t ../config/hum_xen_fug.nwk

