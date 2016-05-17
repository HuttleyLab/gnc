#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l walltime=24:00:00,mem=996GB,ncpus=496
#PBS -V
#usage: map_collection.py [-h]
#                         (-i INPUT_COLLECTION | -I INPUT_COLLECTIONS_FILE)
#                         [-b DB_HOST]
#                         (-o OUTPUT_COLLECTION | -O OUTPUT_COLLECTIONS_FILE)
#                         -f FUNCTION [-L LOG_LEVEL] [-l LOG_NAME] [-s]
#                         [-B BATCH_SIZE] [-k KWARGS_FILE] [-m]
mpiexec -np $PBS_NCPUS $VIRTUALENVWRAPPER_HOOK_DIR/normal/bin/python map_collection.py -i mammals.GNC -b $DB_HOST -L DEBUG -f ml.ml_bootstraps -o mammals.GNC_bootstraps2 -k ../config/dont_use_mpi.json
