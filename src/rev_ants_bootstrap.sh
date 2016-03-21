#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l walltime=5:00:00,mem=4064GB,ncpus=2032
#PBS -V
#usage: map_collection.py [-h]
#                         (-i INPUT_COLLECTION | -I INPUT_COLLECTIONS_FILE)
#                         [-b DB_HOST]
#                         (-o OUTPUT_COLLECTION | -O OUTPUT_COLLECTIONS_FILE)
#                         -f FUNCTION [-L LOG_LEVEL] [-l LOG_NAME] [-s]
#                         [-B BATCH_SIZE] [-k KWARGS_FILE] [-m]
mpiexec -np $PBS_NCPUS python map_collection.py -i ants.CNFGTR -b $DB_HOST -L DEBUG -f ml.ml_bootstraps -o ants.CNFGTR_bootstraps3 -k ../config/dont_use_mpi.json
mpiexec -np $PBS_NCPUS python map_collection.py -i ants.Y98 -b $DB_HOST -L DEBUG -f ml.ml_bootstraps -o ants.Y98_bootstraps3 -k ../config/dont_use_mpi.json