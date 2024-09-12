#!/bin/sh

#SBATCH --account=csmpi
#SBATCH --partition=csmpi_long
#SBATCH --mem=0
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH --output=stream.out

make clean
make
OMP_PLACES=cores numactl --interleave all ./stream_c
