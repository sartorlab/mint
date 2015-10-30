#!/bin/bash
set -e
set -u
set -o pipefail

# All classification BED files are formatted in the same way so this is universal

PROJECT=IDH2mut_v_NBM
classBED=./analysis/classification_comparison/IDH2mut_v_NBM_classification_comparison.bed
cpgResult=`basename ${classBED} .bed`_annot.bed

cd ~/latte/mint/${PROJECT}

################################################################################
# Intersection of classification (comparison, sample, or simple) regions with annotations
# For barplots showing proportion of classification regions in each annotation
# This will take a very long time and require lots of memory

  bedtools intersect -wa -wb \
    -a <(awk -v OFS='\t' '{print $1, $2, $3, $4}' ${classBED}) \
    -b <(cat \
        ../data/annotation/cpg_islands_hg19_ucsc.bed \
        ../data/annotation/cpg_shores_hg19_ucsc.bed \
        ../data/annotation/cpg_shelves_hg19_ucsc.bed \
        ../data/annotation/cpg_inter_hg19_ucsc.bed \
        ../data/annotation/annot_promoters_hg19.bed \
        ../data/annotation/annot_enhancers_hg19.bed \
      | sort -T . -k1,1 -k2,2n) \
  > ./analysis/summary/tables/${cpgResult}
