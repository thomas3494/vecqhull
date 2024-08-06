#!/bin/sh

#SBATCH --account=csmpi
#SBATCH --partition=csmpi_long
#SBATCH --mem=0
#SBATCH --nodelist=cn125
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH --output=check_quickhull.out

if [ "$#" -ne 1 ]; then
    printf 'Usage: %s N\n' "$0" >&2
           '\tN: number of points\n' >&2
    exit 1
fi

n="$1"

(
cd code || exit
make bin/test_quickhull bin/index2points
)

(
cd pbbsbench/benchmarks/convexHull/serialHull || exit
make
)

check()
{
    name="$1"
    ./code/bin/test_quickhull p < data/"$name".in > data/"$name".quickhull.out
    ./pbbsbench/benchmarks/convexHull/serialHull/hull  \
                -o data/"$name".pbbs_serial_indices.out \
                data/"$name".in
    ./code/bin/index2points data/"$name".in \
            data/"$name".pbbs_serial_indices.out > \
            data/"$name".pbbs_serial.out
    if diff data/"$name".pbbs_serial.out data/"$name".quickhull.out; then
        printf 'Same result as pbbs on %s\n', "$name"
    fi
}

check disk_"$n"
check circle_"$n"
check kuzmin_"$n"
