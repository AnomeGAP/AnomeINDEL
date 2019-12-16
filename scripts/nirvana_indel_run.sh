#!/usr/bin/env bash
END=16

echo "ID	Total	Indel	Insertion   Deletion    Tier-1  Tier-2"
for i in $(seq 1 ${END});
do
    str=$(printf "%02d" ${i})
    OUTPUT=`gzcat /Users/chungtsai_su/KSTF/KSTF_NIRVANA/KSTF_06${str}.json.gz | python ./nirvana_indel.py -f /Users/chungtsai_su/data/COSMIC/Census_all.tsv -o /Users/chungtsai_su/KSTF/KSTF_NIRVANA/indel/KSTF_06${str}`
    echo ${OUTPUT}
done
