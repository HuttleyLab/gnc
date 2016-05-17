#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l walltime=10:00:00,mem=2016GB,ncpus=1008
#PBS -v DB_HOST
#PBS -P x36
#usage: map_collection.py [-h]
#                         (-i INPUT_COLLECTION | -I INPUT_COLLECTIONS_FILE)
#                         [-b DB_HOST]
#                         (-o OUTPUT_COLLECTION | -O OUTPUT_COLLECTIONS_FILE)
#                         -f FUNCTION [-L LOG_LEVEL] [-l LOG_NAME] [-s]
#                         [-B BATCH_SIZE] [-k KWARGS_FILE] [-m]
source ~/.profile
workon normal
mpirun -np $PBS_NCPUS $VIRTUAL_ENV/bin/python map_collection_slowly.py -i hum_xen_fug.GNC -b $DB_HOST -L DEBUG -f ml.ml_bootstraps -o hum_xen_fug.GNC_bootstraps2 -k ../config/dont_use_mpi.json
