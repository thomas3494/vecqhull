#!/bin/sh

#SBATCH --account=csmpi
#SBATCH --partition=csmpi_long
#SBATCH --mem=0
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH --output=stream_cn125.out

make clean
make stream_c
OMP_PLACES=cores ./stream_c
