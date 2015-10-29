#!/bin/bash
set -e
set -u
set -o pipefail

PROJECT=IDH2mut_v_NBM
classBED=./analysis/classification_comparison/IDH2mut_v_NBM_classification_comparison.bed
result=`basename ${classBED} .bed`_bp_per_class.txt
cpgResult=`basename ${classBED} .bed`_cpg_annot.bed

cd ~/latte/mint/${PROJECT}

################################################################################
# Distribution of classifications overall and in CpG annotations
# This will take a very long time

  bedtools intersect \
    -a <(awk -v OFS='\t' '{print $1, $2, $3, $4}' ${classBED}) \
    -b <(cat \
        <(awk -v OFS='\t' '{print $1, $2, $3, "island"}' ../data/annotation/cpg_islands_hg19_ucsc.bed) \
        <(awk -v OFS='\t' '{print $1, $2, $3, "shore"}' ../data/annotation/cpg_shores_hg19_ucsc.bed) \
        <(awk -v OFS='\t' '{print $1, $2, $3, "shelf"}' ../data/annotation/cpg_shelves_hg19_ucsc.bed) \
        <(awk -v OFS='\t' '{print $1, $2, $3, "inter"}' ../data/annotation/cpg_inter_hg19_ucsc.bed) \
      | sort -T . -k1,1 -k2,2n) \
    -wa -wb \
  > ${cpgResult}
