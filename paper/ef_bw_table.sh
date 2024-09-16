#!/bin/sh
#
if [ "$#" -ne 3 ]; then
    printf 'Usage: %s N PREC RESULT\n' "$0" >&2
    printf '\tN:      number of points\n' >&2
    printf '\tPREC:   number of significant digits\n' >&2
    printf '\tRESULT: result directory\n' >&2
    exit 1
fi

n="$1"
prec="$2"
result="$3"

print_backend()
{
    backend="$1"
    tmp_kuzmin=$(mktemp)
    ./div.sh $(echo "scale=17; (1600000000 + 4) / 1000000000" | bc) < \
            "${result}/kuzmin_${n}_${backend}.csv" > "$tmp_kuzmin"
    tmp_circle=$(mktemp)
    ./div.sh $(echo "scale=17; (1600000000 + 10985422) / 1000000000" | bc) < \
            "${result}/circle_${n}_${backend}.csv" > "$tmp_circle"
    tmp_disk=$(mktemp)
    ./div.sh $(echo "scale=17; (1600000000 + 1593) / 1000000000" | bc) < \
            "${result}/disk_${n}_${backend}.csv" > "$tmp_disk"
    paste "$tmp_kuzmin" \
          "$tmp_circle" \
          "$tmp_disk" | \
          sed 's/\t/,/g' | ./round.sh "$prec"
}

echo "Imp,Kuzmin,kuzminstd,Circle,circlestd,Disk,diskstd"
printf '%s,' "PBBS"
print_backend pbbs_serial
printf '%s,' "VecQuickhull"
print_backend quickhull
printf '%s,' "PBBS par"
print_backend pbbs_multi
printf '%s,' "VecQuickhull par"
print_backend quickhull_par
