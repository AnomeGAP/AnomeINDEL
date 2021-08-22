#!/usr/bin/env bash

LIST=('A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J')

for (( i=0; i<${#LIST[@]}; i++ ));
do
  echo ${LIST[$i]}
#  bgzip ./Pisces/${LIST[$i]}.vcf
#  bcftools index -f ./Pisces/${LIST[$i]}.vcf.gz
#
#  bgzip ./Mutect2/${LIST[$i]}.vcf
#  bcftools index -f ./Mutect2/${LIST[$i]}.vcf.gz
#
#  bgzip ./DeepVariant/${LIST[$i]}.vcf
#  bcftools sort ./DeepVariant/${LIST[$i]}.vcf.gz -o ./DeepVariant/${LIST[$i]}.sorted.vcf.gz -O z
#  bcftools index -f ./DeepVariant/${LIST[$i]}.sorted.vcf.gz
#
#  bgzip ./HaplotypeCaller/${LIST[$i]}.vcf
#  bcftools index -f ./HaplotypeCaller/${LIST[$i]}.vcf.gz
#
#  bcftools merge --force-samples ./Pisces/${LIST[$i]}.vcf.gz ./Mutect2/${LIST[$i]}.vcf.gz ./DeepVariant/${LIST[$i]}.sorted.vcf.gz ./HaplotypeCaller/${LIST[$i]}.vcf.gz > ./${LIST[$i]}.merged.vcf
  python ~/src/github/AnomeINDEL/scripts/PrecisionFDA-TMB.py -i ${LIST[$i]}.merged.vcf -o output/${LIST[$i]}.final.vcf
done

cd output
rm -rf all.germline_and_somatic.zip
zip -9 all.germline_and_somatic.zip *.final.vcf
