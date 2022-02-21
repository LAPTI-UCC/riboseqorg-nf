#!usr/bin/env bash 

SUPERSET='/home/jack/projects/riboseq_data_processing/data/ribosome_profiling_superset.tsv'
cat $SUPERSET | while read line; do 
    ARRAY=${line}
    echo {$ARRAY[9]}

done 
