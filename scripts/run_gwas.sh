#!/usr/bin/env bash

##Software package
SRC=/mnt/volume_1TB/atsai/src
PLINK=$SRC/plink_linux_x86_64_20200219/plink
SHAPEIT=$SRC/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit
IMPUTE2=$SRC/impute_v2.3.2_x86_64_static/impute2
GTOOL=$SRC/gtool_v0.7.5_x86_64/gtool
VCF2IMPUTE2=$SRC/vcf2impute_legend_haps.pl
PEDMERGER=$SRC/PEDMerger.py

## Dataset
DATA=/mnt/volume_1TB/atsai/data
THE1000GP_PHASE3=$DATA/1000GP_Phase3

## INPUT
INPUT=/mnt/volume_1TB/academia-sinica/data/baseline_grch37
TMP=/mnt/volume_1TB/academia-sinica/tmp
OUTPUT=/mnt/volume_1TB/academia-sinica/output

## SYSTEM
NUM_THREADS=32

##make vcf from bed+bim+fam
#mkdir -p $TMP
#for chr in $(seq 17 17); do
#  $PLINK \
#    --bfile $INPUT \
#    --recode vcf bgz \
#    --chr $chr \
#    --memory 2048 \
#    --out $TMP/chr$chr
#done

##phasing
for chr in $(seq 17 17); do
  REF_GENETICMAP=$THE1000GP_PHASE3/genetic_map_chr$chr\_combined_b37.txt
  REF_HAP=$THE1000GP_PHASE3/1000GP_Phase3_chr$chr.hap.gz
  REF_LEGEND=$THE1000GP_PHASE3/1000GP_Phase3_chr$chr.legend.gz
  REF_SAMPLE=$THE1000GP_PHASE3/1000GP_Phase3.sample
  echo $chr
#  REF_GENETICMAP=/mnt/volume_1TB/academia-sinica/data/genetic_map_chr$chr\_combined_b37.txt
#  REF_HAP=/mnt/volume_1TB/academia-sinica/data/1000GP_Phase3_chr$chr.hap.gz
#  REF_LEGEND=/mnt/volume_1TB/academia-sinica/data/1000GP_Phase3_chr$chr.legend.gz
#  REF_SAMPLE=

#  $SHAPEIT \
#    -check \
#    --input-vcf $TMP/chr$chr.vcf.gz \
#    -M $REF_GENETICMAP \
#    --seed 158277 \
#    --thread $NUM_THREADS \
#    --output-log $TMP/chr$chr.phased
##    --input-ref $REF_HAP $REF_LEGEND $REF_SAMPLE \

#  $SHAPEIT \
#    --input-vcf $TMP/chr$chr.vcf.gz \
#    -M $REF_GENETICMAP \
#    --seed 158277 \
#    --thread $NUM_THREADS \
#    -O $TMP/chr$chr.phased
##    --input-ref $REF_HAP $REF_LEGEND $REF_SAMPLE \
##    --exclude-snp $TMP/chr$chr.phased.snp.strand.exclude \

#  $SHAPEIT \
#    -convert \
#    --input-haps $TMP/chr$chr.phased \
#    --output-vcf $TMP/chr$chr.phased.vcf

#  bgzip $TMP/chr$chr.phased.vcf

done

START=(78 5000079 10000080 15000081 20000082 25000083 30000084 35000085 40000086 45000087 50000088 55000089 60000090 65000091 70000092 75000093 80000094)
END=(5000078 10000079 15000080 20000081 25000082 30000083 35000084 40000085 45000086 50000087 55000088 60000089 65000090 70000091 75000092 80000093 81189911)
for chr in $(seq 17 17); do
  for (( i=0; i<${#START[@]}; i++ )); do
    echo $chr:${START[$i]}-${END[$i]}
#    vcftools --gzvcf $TMP/chr$chr.phased.vcf.gz --chr $chr --from-bp ${START[$i]} --to-bp ${END[$i]} --recode --recode-INFO-all --out $TMP/chr$chr.phased.chunk$i
#    bgzip $TMP/chr$chr.phased.chunk$i.recode.vcf
  done
done

##impute2.3.2
for chr in $(seq 17 17); do
  for (( i=0; i<${#START[@]}; i++ )); do
    echo $chr:${START[$i]}-${END[$i]}
#    perl $VCF2IMPUTE2 -vcf $TMP/chr$chr.phased.chunk$i.recode.vcf.gz -leghap $TMP/chr$chr.phased.chunk$i.leghap -chr $chr
  done
done

for chr in $(seq 17 17); do
  for (( i=0; i<${#START[@]}; i++ )); do
    echo $chr:${START[$i]}-${END[$i]}
#    $IMPUTE2 \
#      -filt_rules_l  'EAS==0' 'TYPE==Multiallelic_SNP' 'TYPE==Multiallelic_INDEL' 'TYPE==Biallelic_INDEL' 'TYPE==Biallelic_DEL' 'TYPE==Biallelic_DUP' 'TYPE==Biallelic_INV' 'TYPE==Biallelic_MNP' 'TYPE==Multiallelic_CNV' 'TYPE==Biallelic_INS:ME:ALU' 'TYPE==Biallelic_INS:ME:LINE1' 'TYPE==Biallelic_INS:ME:SVA' 'TYPE==Biallelic_INS:MT'  -use_prephased_g \
#      -known_haps_g $TMP/chr$chr.phased.chunk$i.leghap.hap.gz \
#      -h $THE1000GP_PHASE3/1000GP_Phase3_chr$chr.hap.gz \
#      -l $THE1000GP_PHASE3/1000GP_Phase3_chr$chr.legend.gz \
#      -m $THE1000GP_PHASE3/genetic_map_chr$chr\_combined_b37.txt \
#      -int ${START[$i]} ${END[$i]}  \
#      -Ne 20000 \
#      -o $TMP/chr$chr.phased.chunk$i.impute2 > $TMP/chr$chr.phased.chunk$i.log
  done
done

mkdir $OUTPUT
#TODO:Please rename the corresponding files to skip this renaming process
for chr in $(seq 17 17); do
  for (( i=0; i<${#START[@]}; i++ )); do
    mv $OUTPUT/chr$chr.chunk$i.ped $OUTPUT/chr$chr.$(printf "chunk%02d.ped" ${i})
    mv $OUTPUT/chr$chr.chunk$i.map $OUTPUT/chr$chr.$(printf "chunk%02d.map" ${i})

  done
done

for chr in $(seq 17 17); do
  for (( i=0; i<${#START[@]}; i++ )); do
    ##gtools
    echo $i
#    $GTOOL -G \
#      --g $TMP/chr$chr.phased.chunk$i.impute2 \
#      --s $TMP/chr$chr.phased.chunk$i.leghap.sample_list \
#      --ped $OUTPUT/chr$chr.chunk$i.ped \
#      --map $OUTPUT/chr$chr.chunk$i.map \
#      --log $OUTPUT/chr$chr.chunk$i.log \
#      --chr $chr \
#      --phenotype plink_pheno  --threshold 0.9
  done
  # merge ped files
  rm -rf $OUTPUT/chr$chr.merged.ped
  python3 $PEDMERGER -i $OUTPUT -o $OUTPUT/chr$chr.merged.ped
  # merge map files
  rm -rf $OUTPUT/chr$chr.merged.map
  for (( i=0; i<${#START[@]}; i++ )); do
    cat $OUTPUT/chr$chr.$(printf "chunk%02d.map" ${i}) >> $OUTPUT/chr$chr.merged.map
  done
done

#TODO:concat all of chunk per chromosome (Deprecated)
#for chr in $(seq 17 17); do
#  echo $chr
#  for (( i=0; i<${#START[@]}; i++ )); do
#    if [ $i -eq 0 ]
#    then
#      cp $TMP/chr$chr.phased.chunk$i.impute2 $TMP/chr$chr.phased.impute2
#    else
#      cat $TMP/chr$chr.phased.chunk$i.impute2 >> $TMP/chr$chr.phased.impute2
#    fi
#  done
  ##gtools
#  $GTOOL -G \
#    --g $TMP/chr$chr.phased.impute2 \
#    --s $TMP/chr$chr.phased.sample \
#    --ped $OUTPUT/chr$chr.ped \
#    --map $OUTPUT/chr$chr.map \
#    --log $OUTPUT/chr$chr.log \
#    --chr $chr \
#    --phenotype plink_pheno  --threshold 0.9
#done

