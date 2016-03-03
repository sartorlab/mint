#!/bin/bash
set -e
set -u
set -o pipefail

# All classification BED files are formatted in the same way so this is universal

# Arguments
# -project The project name given to init_project.sh
# -classBED The path to the classification BED file NOTE: Need full path
PROJECT=$2
classBED=$4

# The resulting file from the intersection
annotResult=`basename ${classBED} .bed`_annot.bed

cd ~/latte/mint/${PROJECT}

################################################################################
# Intersection of classification (comparison, sample, or simple) regions with annotations
# For barplots showing proportion of classification regions in each annotation
# This will take a very long time and require lots of memory

  echo "Computing intersections over annotations for ${annotResult}"
  bedtools intersect -wa -wb \
    -a <(awk -v OFS='\t' '{print $1, $2, $3, $4}' ${classBED}) \
    -b <(cat \
        ../data/annotation/annot_cpg_islands_hg19.bed \
        ../data/annotation/annot_cpg_shores_hg19.bed \
        ../data/annotation/annot_cpg_shelves_hg19.bed \
        ../data/annotation/annot_cpg_inter_hg19.bed \
        ../data/annotation/annot_promoters_hg19.bed \
        ../data/annotation/annot_enhancers_hg19.bed \
      | sort -T . -k1,1 -k2,2n) \
  > ./analysis/summary/tables/${annotResult}
