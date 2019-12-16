#!/usr/bin/env bash

#python ./LengthFilter.py -i ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_insertion.q1.fa
#python ./LengthFilter.py -i ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_deletion.q1.fa
#python ./LengthFilter.py -i ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_softclipping.q1.fa
#python ./LengthFilter.py -i ~/NA12878-novaseq/v1.0.2/U0.fa
#
#~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_insertion.q1.fa-s.fa
#~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_insertion.q1.fa-m.fa
#~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_insertion.q1.fa-l.fa
#
#~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_deletion.q1.fa-s.fa
#~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_deletion.q1.fa-m.fa
#~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_deletion.q1.fa-l.fa
#
#~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_softclipping.q1.fa-s.fa
#~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_softclipping.q1.fa-m.fa
#~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_softclipping.q1.fa-l.fa

~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/U0.fa-l.fa
~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/U0.fa-m.fa
~/src/RepeatMasker/RepeatMasker ~/NA12878-novaseq/v1.0.2/U0.fa-s.fa
