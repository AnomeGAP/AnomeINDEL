LIST=('A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J')
TMB=('5.15' '4.28' '19.66' '8.14' '8.56' '6.09' '9.2' '31.38' '25.16' '7.12')

#for (( j=1; j<=3; j++ ));
#do
#  echo Panel${j}
#  for (( i=0; i<${#LIST[@]}; i++ ));
#  do
#    # echo ${LIST[$i]} ${j}
#    vcftools --vcf ~/data/PrecisionFDA/TMB/TMB2/${LIST[$i]}.vcf --bed ~/data/PrecisionFDA/TMB/TMB2/panel${j}.bed --out ~/data/PrecisionFDA/TMB/TMB2/${LIST[$i]}.p${j} --recode --recode-INFO-all 2>&1 | grep Sites
#    grep -v "#" ~/data/PrecisionFDA/TMB/TMB2/${LIST[$i]}.p${j}.recode.vcf | grep PASS | wc -l
#  done
#done

#for (( i=0; i<${#LIST[@]}; i++ ));
#do
#  echo "${TMB[$i]}"
#  python3 ~/src/github/AnomeINDEL/scripts/PrecisionFDA-TMB-VCF.py -i ~/data/PrecisionFDA/TMB/TMB2/${LIST[$i]}.vcf -a ~/data/PrecisionFDA/TMB/TMB2-VEP/${LIST[$i]}
#done
#
#for (( j=1; j<=3; j++ ));
#do
#  echo Panel${j}
#  touch ~/data/PrecisionFDA/TMB/TMB2/training${j}.tsv
#  echo "" > ~/data/PrecisionFDA/TMB/TMB2/training${j}.tsv
#  for (( i=0; i<${#LIST[@]}; i++ ));
#  do
#    # echo ${LIST[$i]} ${j}
#    echo "${TMB[$i]}" >> ~/data/PrecisionFDA/TMB/TMB2/training${j}.tsv
#    python3 ~/src/github/AnomeINDEL/scripts/PrecisionFDA-TMB-VCF.py -i ~/data/PrecisionFDA/TMB/TMB2/${LIST[$i]}.p${j}.recode.vcf -a ~/data/PrecisionFDA/TMB/TMB2-VEP/${LIST[$i]}.p${j}.recode >> ~/data/PrecisionFDA/TMB/TMB2/training${j}.tsv
#  done
#done
#
#exit

#NOTE: need to fix the new-line '\n' problems.

SVM=~/src/github/libsvm
#for (( j=1; j<=3; j++ ));
#do
#  echo Panel${j}
#  $SVM/svm-scale -s ~/data/PrecisionFDA/TMB/TMB2/training${j}.parameter ~/data/PrecisionFDA/TMB/TMB2/training${j}.tsv > ~/data/PrecisionFDA/TMB/TMB2/training${j}.scale
#  # python $SVM/tools/subset.py -s 1 ~/data/PrecisionFDA/TMB/TMB2/training${j}.scale 9 ~/data/PrecisionFDA/TMB/TMB2/training${j}.scale.t1 ~/data/PrecisionFDA/TMB/TMB2/training${j}.scale.v1
#  # python $SVM/tools/subset.py -s 1 ~/data/PrecisionFDA/TMB/TMB2/training${j}.scale 9 ~/data/PrecisionFDA/TMB/TMB2/training${j}.scale.t2 ~/data/PrecisionFDA/TMB/TMB2/training${j}.scale.v2
#  # python $SVM/tools/subset.py -s 1 ~/data/PrecisionFDA/TMB/TMB2/training${j}.scale 9 ~/data/PrecisionFDA/TMB/TMB2/training${j}.scale.t3 ~/data/PrecisionFDA/TMB/TMB2/training${j}.scale.v3
#  python $SVM/tools/grid.py -log2c -3,3,1 -log2g -3,0,1 -v 25 -svmtrain $SVM/svm-train -gnuplot gnuplot -s 3 ~/data/PrecisionFDA/TMB/TMB2/training${j}.scale
#  # panel1 - 0.03125 0.0625 182.41
#  # panel2 - 0.03125 0.0625 182.76
#  # panel3 - 0.03125 0.0625 182.25
#done

#${SVM}/svm-train -s 3 -c 8 -g 0.25 ~/data/PrecisionFDA/TMB/TMB2/training1.scale ~/data/PrecisionFDA/TMB/TMB2/training1.model
#${SVM}/svm-train -s 3 -c 8 -g 1 ~/data/PrecisionFDA/TMB/TMB2/training2.scale ~/data/PrecisionFDA/TMB/TMB2/training2.model
#${SVM}/svm-train -s 3 -c 8 -g 0.5 ~/data/PrecisionFDA/TMB/TMB2/training3.scale ~/data/PrecisionFDA/TMB/TMB2/training3.model

# ${SVM}/svm-predict ~/data/PrecisionFDA/TMB/TMB2/training1.scale  ~/data/PrecisionFDA/TMB/TMB2/training1.model ~/data/PrecisionFDA/TMB/TMB2/training1.out > /dev/null 2>&1
# ${SVM}/svm-predict ~/data/PrecisionFDA/TMB/TMB2/training2.scale  ~/data/PrecisionFDA/TMB/TMB2/training2.model ~/data/PrecisionFDA/TMB/TMB2/training2.out > /dev/null 2>&1
# ${SVM}/svm-predict ~/data/PrecisionFDA/TMB/TMB2/training3.scale  ~/data/PrecisionFDA/TMB/TMB2/training3.model ~/data/PrecisionFDA/TMB/TMB2/training3.out > /dev/null 2>&1


##Testing
SAMPLES=('01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29')
PANEL=(   '2'  '3'  '1'  '2'  '3'  '1'  '2'  '3'  '1'  '2'  '3'  '1'  '2'  '3'  '1'  '2'  '3'  '1'  '2'  '3'  '1'  '2'  '3'  '1'  '2'  '3'  '1'  '2'  '3')
for (( i=0; i<29; i++ ));
do
#  echo Sample${SAMPLES[$i]}
#  gunzip ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor${SAMPLES[$i]}.vcf.gz
#  echo 10 > ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor${SAMPLES[$i]}.tsv
#  python3 ~/src/github/AnomeINDEL/scripts/PrecisionFDA-TMB-VCF.py -i ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor${SAMPLES[$i]}.vcf -a ~/data/PrecisionFDA/TMB/TMB2-Testset-VEP/tumor${SAMPLES[$i]} >> ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor${SAMPLES[$i]}.tsv
  ##NOTE: need manual effort to remove the new-line
  ${SVM}/svm-scale -r ~/data/PrecisionFDA/TMB/TMB2/training${PANEL[${i}]}.parameter ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor${SAMPLES[$i]}.tsv > ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor${SAMPLES[$i]}.scale
  ${SVM}/svm-predict ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor${SAMPLES[$i]}.scale  ~/data/PrecisionFDA/TMB/TMB2/training${PANEL[$i]}.model ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor${SAMPLES[$i]}.out > /dev/null 2>&1
  cat ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor${SAMPLES[$i]}.out
done

#vcftools --vcf ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor${SAMPLES[$i]}.vcf --bed ~/data/PrecisionFDA/TMB/TMB2/panel${PANEL[${i}]}.bed --out ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor${SAMPLES[$i]}.p${j} --recode --recode-INFO-all 2>&1 | grep Sites
#vcftools --vcf ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor01.vcf --bed ~/data/PrecisionFDA/TMB/TMB2/panel2.bed --out ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor01.p2 --recode --recode-INFO-all 2>&1 | grep Sites
#vcftools --vcf ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor02.vcf --bed ~/data/PrecisionFDA/TMB/TMB2/panel3.bed --out ~/data/PrecisionFDA/TMB/TMB2-Testset/tumor01.p3 --recode --recode-INFO-all 2>&1 | grep Sites
