#!/usr/bin/env bash

for p in "${@}"; do
   printf '<%s>\n' "${p}" >> subparams.txt
done

cat ${19} >> subparams.txt

echo "---" >> subparams.txt

exec "${@}"
