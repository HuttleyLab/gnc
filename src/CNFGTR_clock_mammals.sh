#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l walltime=05:00:00,mem=4096GB,ncpus=2048
#PBS -V
#usage: map_collection.py [-h]
#                         (-i INPUT_COLLECTION | -I INPUT_COLLECTIONS_FILE)
#                         [-b DB_HOST]
#                         (-o OUTPUT_COLLECTION | -O OUTPUT_COLLECTIONS_FILE)
#                         -f FUNCTION [-L LOG_LEVEL] [-l LOG_NAME] [-s]
#                         [-B BATCH_SIZE] [-k KWARGS_FILE] [-m]
mpiexec -np $PBS_NCPUS python map_collection.py -i mammals.data -b $DB_HOST -L DEBUG -f clock.ml -o mammals.CNFGTR_clock -k ../config/CNFGTR_clock_mammals.json -s