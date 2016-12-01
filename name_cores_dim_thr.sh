#!/bin/bash

jobname=$1
cores=$2
d=$3 # Dimensions
t=$4 # Threshold

jobname_date=$jobname-run-$(date +%Y%m%d).$(date +%H%M)

echo "---Submitted as---"
echo "Name: $jobname_date"
echo "Cores: $cores"
echo "Dimensions: $d"
echo "Threshold: $t"
echo "full call: qsub -N $jobname_date -q bio,pub64,free* -pe openmp $cores-$cores -m beas \
python3 ~/data/CC4.0/run_OSF.py -i $jobname.csv -o $jobname_date -d $d -p $cores -l 10000 -t $t \
> $jobname_date.log"

qsub -N $jobname_date -q bio,pub64,free* -pe openmp $cores-$cores -m beas \
python3 ~/data/CC4.0/run_OSF.py -i $jobname.csv -o $jobname_date -d $d -p $cores -l 10000 -t $t \
> $jobname_date.log