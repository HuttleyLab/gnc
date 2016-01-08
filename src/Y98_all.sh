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
#                         [-k KWARGS_FILE] [-m]
mpiexec -np $PBS_NCPUS python map_collection.py -i ants.data -b $DB_HOST -L DEBUG -f ml.ml -o ants.Y98 -k ../config/Y98.json
mpiexec -np $PBS_NCPUS python map_collection.py -i mammals.data -b $DB_HOST -L DEBUG -f clock.ml -o mammals.Y98_clock -k ../config/Y98_clock_mammals.json -s
mpiexec -np $PBS_NCPUS python map_collection.py -i introns.data -b $DB_HOST -L DEBUG -f ml.ml -o introns.Y98 -k ../config/Y98_no_stop.json
mpiexec -np $PBS_NCPUS python map_collection.py -i mammals.data -b $DB_HOST -L DEBUG -f ml.ml -o mammals.Y98 -k ../config/Y98.json -s
mpiexec -np $PBS_NCPUS python map_collection.py -i hum_xen_fug.data -b $DB_HOST -L DEBUG -f ml.ml -o hum_xen_fug.Y98 -k ../config/Y98.json -s
