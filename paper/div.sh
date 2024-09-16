#!/bin/sh

if [ "$#" -ne 1 ]; then
    printf 'Usage: %s <n>\n' "$0"
    printf '\tTakes a csv on stdin, transforms it by f(x) = x / n\n'
    printf '\tand returns it on stdout\n'
    exit 1
fi

n="$1"

awk -F "," -v var="$n" '{for (i = 1; i <= NF - 1; i++) { \
                                        printf("%e, ", var / $i)    \
                                    }                               \
                                    printf("%e\n", var / $NF)       \
                                }' < /dev/stdin
