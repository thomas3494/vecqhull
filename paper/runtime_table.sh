#!/bin/sh

if [ "$#" -ne 2 ]; then
    printf 'Usage: %s N PREC\n' "$0" >&2
    printf '\tN:      number of points\n' >&2
    printf '\tPREC:   number of significant digits\n' >&2
    exit 1
fi

n="$1"
prec="$2"

print_backend()
{
    backend="$1"
    paste "results/kuzmin_${n}_${backend}.csv" \
          "results/circle_${n}_${backend}.csv" \
          "results/disk_${n}_${backend}.csv" | \
          sed 's/\t/,/g' | ./round.sh "$prec"
}

echo "Imp,Kuzmin,kuzminstd,Circle,circlestd,Disk,diskstd"
printf '%s,' "BlockQuickhull"
print_backend quickhull
printf '%s,' "PBBS"
print_backend pbbs_serial
printf '%s,' "BlockQuickhull par"
print_backend quickhull_par
printf '%s,' "PBBS par"
print_backend pbbs_multi
