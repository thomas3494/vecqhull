#!/bin/sh

if [ "$#" -ne 1 ]; then
    printf 'Usage %s <N>\n' "$0"
    exit 1
fi

n="$1"

hash qconvex

tempfile=$(mktemp)

printf '2 disk D2 %s\n%s\n' "$n" "$n" > "$tempfile"

tail -n +2 ./data/disk_"$n".in | cat "$tempfile" - | qconvex | grep CPU | sed 's/.*: //g'
