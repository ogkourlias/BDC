#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2

source ./venvs/bdc/bin/activate
num_cores=2

rm timings.csv
touch timings.csv
echo "runtime,start,end" >> timings.csv
parallel -j 2 python assignment6.py -f "data/{}/chr21.vcf.gz" -r data/ref_vcf/chr21.vcf.gz -ih data/ref_vcf/chr21_headers.txt -ch data/ref_vcf/chr21_headers.txt -chr chr21 -o output/{}.csv -n 1000 ::: subset1 subset2
