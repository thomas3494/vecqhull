#!/bin/sh

#SBATCH --account=csmpi
#SBATCH --partition=csmpi_long
#SBATCH --mem=0
#SBATCH --nodelist=cn125
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --output=gen_input.out

if [ "$#" -ne 1 ]; then
    printf 'Usage: ./gen_input.sh N\n' >&2
           '\tN: number of points\n' >&2
    exit 1
fi

n="$1"

(cd pbbsbench || exit 1
make geometryData
./testData/geometryData/randPoints -s -d 2 "$n" ../data/disk_"$n".in
./testData/geometryData/randPoints -S -d 2 "$n" ../data/circle_"$n".in
./testData/geometryData/randPoints -p -d 2 "$n" ../data/kuzmin_"$n".in)

(cd code/examples
make text2bin)

code/examples/text2bin < data/disk_"$n".in > data/disk_"$n".bin
code/examples/text2bin < data/circle_"$n".in > data/circle_"$n".bin
code/examples/text2bin < data/kuzmin_"$n".in > data/kuzmin_"$n".bin
