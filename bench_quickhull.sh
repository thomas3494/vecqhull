#!/bin/sh

#SBATCH --account=csmpi
#SBATCH --partition=csmpi_long
#SBATCH --nodelist=cn125
#SBATCH --mem=0
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH --output=bench_quickhull.out

if [ "$#" -ne 3 ]; then
    printf 'Usage: %s N ITER OUTDIR\n' "$0" >&2
    printf '\tN:      number of points\n' >&2
    printf '\tITER:   number of times to repeat the experiment\n' >&2
    printf '\tOUTDIR: directory to store result\n' >&2
    exit 1
fi

n="$1"
iter="$2"
outdir="$3"

mkdir -p "$outdir"

(
cd code || exit
make bin/test_quickhull
)

bench()
{
    name="$1"
    {
        i=1
        while [ $i -le "$iter" ]
        do
            ./code/bin/test_quickhull b < data/"$name".bin
            i=$(( i + 1 ))
        done
    } | awk '{
               for (i = 1; i <= NF; i++) {
                   b[i] = a[i] + ($i - a[i]) / NR;
                   q[i] += ($i - a[i]) * ($i - b[i]);
                   a[i] = b[i];
               }
             } END {
               printf "%f, %f", a[1], sqrt(q[1] / NR);
               for (i = 2; i <= NF; i++) {
                   printf ", %f, %f", a[i], sqrt(q[i] / NR);
               }
               print "";
             }' > "${outdir}/${name}_quickhull.csv"
}

bench disk_"$n"
bench circle_"$n"
bench kuzmin_"$n"
