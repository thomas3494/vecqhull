#!/bin/sh

if [ "$#" -ne 1 ]; then
    printf 'Usage: %s <d>\n' "$0"
    printf '\t<d>: number of significant digits\n'
    printf '\tTakes a csv on stdin, and returns the rounded csv on stdout\n'
    exit 1
fi

prec="$(( $1 - 1 ))"

awk -F "," -v format="%.${prec}e" '{for (i = 1; i <= NF - 1; i++) { \
                                        printf(format ", ", $i)     \
                                    }                               \
                                    printf(format "\n", $NF)        \
                                }' < /dev/stdin
