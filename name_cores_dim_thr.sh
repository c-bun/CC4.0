#!/bin/bash

jobname=$1
cores=$2
let pcores=$cores-1
d=$3 # Dimensions
t=$4 # Threshold

jobname_date="$jobname-c$cores-d$d-t$t-$(date +%Y%m%d).$(date +%H%M)"

echo "---Submitted as---"
echo "Name: $jobname_date"
echo "Cores: $cores"
echo "Dimensions: $d"
echo "Threshold: $t"
echo "full call: "
cat << _EOF_ > temp.sh
#!/bin/bash

python3 ~/data/CC4.0/run_OSF.py -i ${jobname}.csv -o ${jobname_date}.csv -d $d -p $pcores -l 10000 -t $t \
> ${jobname_date}.log
_EOF_

cat temp.sh
qsub -N $jobname_date -q bio,pub64,free* -pe openmp $cores-$cores -m bea temp.sh
sleep 5
rm temp.sh
echo "---Done---"
