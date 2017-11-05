#!/bin/bash

jobname=$1
cores=$2
let pcores=$cores-1
dm=$3 # M Dimensions
dn=$4 # N Dimensions
t=$5 # Threshold

jobname_date="$jobname-c$cores-d$dm$dn-t$t-$(date +%Y%m%d).$(date +%H%M)"

echo "---Submitted as---"
echo "Name: $jobname_date"
echo "Cores: $cores"
echo "Dimensions: $dm x $dn"
echo "Threshold: $t"
echo "full call: "
cat << _EOF_ > temp.sh
#!/bin/bash

python3 ~/data/CrossCompare/run_OSF.py -i ${jobname}.csv -o ${jobname_date}.pkl -m $dm -n $dn -p $pcores -l 10000 -t $t -k -f \
> ${jobname_date}.log
_EOF_

cat temp.sh
qsub -N $jobname_date -q bio,pub64,free* -pe openmp $cores-$cores -m bea temp.sh
sleep 5
rm temp.sh
echo "---Done---"
