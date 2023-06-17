#!/bin/bash -l
#SBATCH --job-name=mpi_fastq   # create a name for your job
#SBATCH --ntasks=24               # total number of tasks
#SBATCH --mem-per-cpu=0        # memory per cpu-core
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)

source /commons/conda/conda_load.sh

mpiexec -n 24 python3 assignment4.py /commons/Themas/Thema12/HPC/rnaseq.fastq