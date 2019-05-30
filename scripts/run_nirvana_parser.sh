#!/usr/bin/env bash

END=16

for i in $(seq 1 ${END});
do
    # echo $i;
    str=$(printf "%02d" ${i})
    zcat ./KSTF_NIRVANA/KSTF_06${str}.json.gz | python ./nirvana_parser.py -f ./Census_all.tsv
done
