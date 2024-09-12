#!/bin/sh

#SBATCH --account=csmpi
#SBATCH --partition=csmpi_long
#SBATCH --mem=0
#SBATCH --nodelist=cn125
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH --output=quickhull_par.out

if [ "$#" -ne 1 ]; then
    printf 'Usage: %s N\n' "$0" >&2
           '\tN: number of points\n' >&2
    exit 1
fi

(
cd code || exit
make bin/test_quickhull_par
)

./code/bin/test_quickhull_par b s < $1
