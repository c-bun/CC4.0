#!/bin/bash

config=$1 # Path to the config file
cores=$2 # Number of cores to request
let processes=$cores-1

jobname_date="$config-c$cores-$(date +%Y%m%d).$(date +%H%M)"

echo "Using config file at $config:"
cat $config
echo "Requesting $cores cores"
echo "Full call:"
cat << _EOF_ > temp.sh
#!/bin/bash
python3 ~/data/CrossCompare/run_OSF.py -c $config -p $processes > ${jobname_date}.log
_EOF_

cat temp.sh
qsub -N $jobname_date -q bio,pub64,free* -pe openmp $cores-$cores -m bea temp.sh
sleep 5
rm temp.sh
echo "---Done---"
