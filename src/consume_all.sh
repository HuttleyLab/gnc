#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l walltime=05:00:00,mem=1096GB,ncpus=512
#PBS -V
# usage: consume.py [-h] -i INPUT_DIRECTORY -o OUTPUT_COLLECTION [-b DB_HOST]
#                   [-L LOG_LEVEL] [-l LOG_NAME] [-c CODON_POSITION]
mpiexec -np $PBS_NCPUS python consume.py -i $SHORT/data/mammals -o mammals.data -L DEBUG -b $DB_HOST -a 1500 -t ../config/mammals.nwk
mpiexec -np $PBS_NCPUS python consume.py -i $SHORT/../data/ants -o ants.data -L DEBUG -b ${DB_HOST} -a 1500 -t ../config/ants.nwk
mpiexec -np $PBS_NCPUS python consume.py -i $SHORT/../data/hum_xen_fug -o hum_xen_fug.data -L DEBUG -b ${DB_HOST} -a 1500 -t ../config/hum_xen_fug.nwk
mpiexec -np $PBS_NCPUS python consume.py -i $SHORT/data/human_macaque_marmoset_introns/data/introns -o introns.data -L DEBUG -b $DB_HOST -a 1500 -t ../config/introns.nwk
