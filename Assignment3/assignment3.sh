#!/bin/bash -l

#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

source /commons/conda/conda_load.sh
export FILE="fastq/test1.fastq";

parallel --jobs 20 --pipepart --block -4 --regexp \
         --recstart '@.*(/1| 1:.*)\n[A-Za-z\n\.~]' --recend '\n' \
         -a $FILE python3 assignment3.py --chunker | \
         python3 assignment3.py --combine
