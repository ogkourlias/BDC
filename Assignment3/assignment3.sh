#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2

source /commons/conda/conda_load.sh

input_file="SRR22537909.fasta"
output_csv="output.csv"
num_cores=4

mkdir work
mkdir scores
python3 assignment3.py -s -i $input_file -ch 50 -n 4
find work/ -type f -print | parallel -j $num_cores python3 assignment3.py -c -seq {}
python3 assignment3.py -calc -co $output_csv -if "$(find scores/ -type f -exec echo {} \;)"
rm -rf scores
rm -rf work
