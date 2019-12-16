#!/usr/bin/env bash

## Example command on cgbs
## bash ./ScanFolder.sh ../data/WGBS

END=10

ptn_r1="*_R1*"
ptn_r2="*_R2*"

for i in $(seq 1 ${END});
do
    pattern="$(printf "*NY%02d*" ${i})"
    fastq_r1="$(printf "$1/NY%02d_R1.fq.gz" ${i})"
    fastq_r2="$(printf "$1/NY%02d_R2.fq.gz" ${i})"
    rm -rf $fastq_r1
    rm -rf $fastq_r2
    for folder in `ls $1` ; do
        #echo $folder
        if [ -d "$1/$folder" ] ; then
            for subfolder in `ls $1/$folder` ; do
                if [ -d "$1/$folder/$subfolder" ] ; then
                    for entry in `ls $1/$folder/$subfolder` ; do
                        if [[ $entry == $pattern ]]
                        then
                            for items in `ls $1/$folder/$subfolder/$entry` ; do
                                if [[ $items == $ptn_r1 ]]
                                then
                                    echo -e "R1\t$1/$folder/$subfolder/$entry/$items"
                                    cat $1/$folder/$subfolder/$entry/$items >> $fastq_r1
                                elif [[ $items == $ptn_r2 ]]
                                then
                                    echo -e "R2\t$1/$folder/$subfolder/$entry/$items"
                                    cat $1/$folder/$subfolder/$entry/$items >> $fastq_r2
                                fi
                            done
                        fi
                    done
                fi
            done
        fi
    done
done
