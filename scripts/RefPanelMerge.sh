#!/usr/bin/env bash

rm -rf /mnt/volume_1TB/atsai/script/merge
mkdir -p /mnt/volume_1TB/atsai/script/merge/log

while IFS= read -r line; do
  declare -a items
  items=(`echo "$line"`)

  if [ ${items[0]} != '#' ]
  then
    /mnt/volume_1TB/atsai/src/impute_v2.3.2_x86_64_static/impute2 \
      -allow_large_regions \
      -m /seqslab/mnt/gwas-reference/genetic_map_hg38_withX/genetic_map_${items[0]}_b38.txt \
      -h /seqslab/mnt/gwas-reference/haplotype_ref_for_merge/${items[0]}_1kg_b38_for_merge.hap \
         /seqslab/mnt/gwas-reference/haplotype_TWBref_for_merge/${items[0]}_TWB_b38_for_merge.hap \
      -l /seqslab/mnt/gwas-reference/haplotype_ref_for_merge/${items[0]}_1kg_b38_for_merge.legend \
         /seqslab/mnt/gwas-reference/haplotype_TWBref_for_merge/${items[0]}_TWB_b38_for_merge.legend \
      -merge_ref_panels \
      -merge_ref_panels_output_ref /mnt/volume_1TB/atsai/script/merge/${items[0]}_${items[3]} \
      -int ${items[1]} ${items[2]} \
      -Ne 20000 \
      -buffer 500 \
      -include_buffer_in_output \
      -o /mnt/volume_1TB/atsai/script/merge/log/${items[0]}_${items[3]}.log > /dev/null &

    if [ `echo "${items[3]} % 55" | bc` -eq 54 ]
    then
      echo ${items[0]} ${items[1]} ${items[2]} ${items[3]}
      while [ `ls -al /mnt/volume_1TB/atsai/script/merge/*.hap | wc -l` -le ${items[3]} ]
      do
        echo "sleep"
        sleep 300
      done
    fi
  fi
done < /seqslab/mnt/system/bed/38/contiguous_unmasked_regions_9406_parts

##ERROR: wrong order for two reference panel
#/mnt/volume_1TB/atsai/src/impute_v2.3.2_x86_64_static/impute2 -allow_large_regions -m /seqslab/mnt/gwas-reference/genetic_map_hg38_withX/genetic_map_chr1_b38.txt -h /seqslab/mnt/gwas-reference/haplotype_TWBref_for_merge/chr1_TWB_b38_for_merge.hap /seqslab/mnt/gwas-reference/haplotype_ref_for_merge/chr1_1kg_b38_for_merge.hap -l /seqslab/mnt/gwas-reference/haplotype_TWBref_for_merge/chr1_TWB_b38_for_merge.legend /seqslab/mnt/gwas-reference/haplotype_ref_for_merge/chr1_1kg_b38_for_merge.legend -merge_ref_panels -merge_ref_panels_output_ref chr1_merge -int 10416 1044017 -Ne 20000 -buffer 500 -include_buffer_in_output -o name_of_log_file

#/mnt/volume_1TB/atsai/src/impute_v2.3.2_x86_64_static/impute2 -allow_large_regions -m /seqslab/mnt/gwas-reference/genetic_map_hg38_withX/genetic_map_chr1_b38.txt -h /seqslab/mnt/gwas-reference/haplotype_ref_for_merge/chr1_1kg_b38_for_merge.hap /seqslab/mnt/gwas-reference/haplotype_TWBref_for_merge/chr1_TWB_b38_for_merge.hap -l /seqslab/mnt/gwas-reference/haplotype_ref_for_merge/chr1_1kg_b38_for_merge.legend /seqslab/mnt/gwas-reference/haplotype_TWBref_for_merge/chr1_TWB_b38_for_merge.legend -merge_ref_panels -merge_ref_panels_output_ref /mnt/volume_1TB/atsai/script/chr1_merge0 -int 10416 1044017 -Ne 20000 -buffer 500 -include_buffer_in_output -o name_of_log_file0

