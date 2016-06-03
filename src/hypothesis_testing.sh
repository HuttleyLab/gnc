#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l walltime=5:00:00,mem=4096GB,ncpus=2048
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
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i mammals.data -b $DB_HOST -L DEBUG -f omega.ml -o mammals.GNC_null -k ../config/GNC_neutral.json
echo 'done mammals GNC'
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i mammals.data -b $DB_HOST -L DEBUG -f omega.ml -o mammals.Y98_null -k ../config/Y98_neutral.json
echo 'done mammals Y98'
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i mammals.data -b $DB_HOST -L DEBUG -f omega.ml -o mammals.CNFGTR_null -k ../config/CNFGTR_neutral.json
echo 'done mammals CNFGTR'
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i ants.data -b $DB_HOST -L DEBUG -f omega.ml -o ants.GNC_null -k ../config/GNC_neutral.json
echo 'done ants GNC'
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i ants.data -b $DB_HOST -L DEBUG -f omega.ml -o ants.Y98_null -k ../config/Y98_neutral.json
echo 'done ants Y98'
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i ants.data -b $DB_HOST -L DEBUG -f omega.ml -o ants.CNFGTR_null -k ../config/CNFGTR_neutral.json
echo 'done ants CNFGTR'
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i hum_xen_fug.data -b $DB_HOST -L DEBUG -f omega.ml -o hum_xen_fug.GNC_null -k ../config/GNC_neutral.json
echo 'done hxf GNC'
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i hum_xen_fug.data -b $DB_HOST -L DEBUG -f omega.ml -o hum_xen_fug.Y98_null -k ../config/Y98_neutral.json
echo 'done hxf Y98'
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i hum_xen_fug.data -b $DB_HOST -L DEBUG -f omega.ml -o hum_xen_fug.CNFGTR_null -k ../config/CNFGTR_neutral.json
echo 'done hxf CNFGTR'
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i introns.data -b $DB_HOST -L DEBUG -f omega.ml -o introns.GNC_null -k ../config/GNC_neutral_no_stop.json
echo 'done introns GNC'
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i introns.data -b $DB_HOST -L DEBUG -f omega.ml -o introns.Y98_null -k ../config/Y98_neutral_no_stop.json
echo 'done introns Y98'
mpiexec -np $PBS_NCPUS python map_collection_slowly.py -i introns.data -b $DB_HOST -L DEBUG -f omega.ml -o introns.CNFGTR_null -k ../config/CNFGTR_neutral_no_stop.json
echp 'done introns CNFGTR'
