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
echo "full call: qsub -N $jobname_date -q bio,pub64,free* -pe openmp $cores-$cores -m beas"
echo "python3 ~/data/CC4.0/run_OSF.py -i $jobname.csv -o $jobname_date -d $d -p $cores -l 10000 -t $t \
> $jobname_date.log"

cat << _EOF_ > temp.sh
#!/bin/bash

python3 ~/data/CC4.0/run_OSF.py -i $jobname.csv -o $jobname_date -d $d -p $cores -l 10000 -t $t \
> $jobname_date.log
_EOF_

qsub -N $jobname_date -q bio,pub64,free* -pe openmp $cores-$cores -m beas temp.sh
sleep 2
rm temp.sh
echo "---Done---"
