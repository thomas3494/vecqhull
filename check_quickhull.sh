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
make clean
make uninstall
make install
)

(
cd code/examples || exit
make clean
make test_quickhull index2points
)

(
cd pbbsbench/benchmarks/convexHull/serialHull || exit
make
)

check()
{
    name="$1"

    echo ""
    echo "=============== $name ==============="

    ./code/examples/test_quickhull p b < data/"$name".bin > data/"$name".quickhull.out
    ./pbbsbench/benchmarks/convexHull/serialHull/hull  \
                -o data/"$name".pbbs_serial_indices.out \
                data/"$name".in
    ./code/examples/index2points data/"$name".in \
            data/"$name".pbbs_serial_indices.out > \
            data/"$name".pbbs_serial.out
    if diff data/"$name".pbbs_serial.out data/"$name".quickhull.out; then
        echo "Same result as pbbs on $name"
    fi

    echo "=============== $name ==============="
}

check disk_"$n"
check circle_"$n"
check kuzmin_"$n"
