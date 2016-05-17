#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l walltime=5:00:00,mem=2208GB,ncpus=1104
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
mpiexec -np $PBS_NCPUS python map_collection.py -i mammals.GNC_unrooted_sims -b $DB_HOST -L DEBUG -f ml.ml -o mammals.GNC_unrooted_sims_fits -k ../config/GNC.json -s
