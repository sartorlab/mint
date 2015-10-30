#!/bin/bash
set -e
set -u
set -o pipefail

PROJECT=IDH2mut_v_NBM
humanID=IDH2mut_1_hmc

cd ~/latte/mint/${PROJECT}

################################################################################
# Take the mean of pulldown coverage over each region in each annotation
# For histograms displaying distributions in different annotations
# From the perspective of the annotation regions
# Seems to require >=48GB of RAM and a long time

  for annot in cpg_islands cpg_shores cpg_shelves cpg_inter promoters enhancers windows_5kb
  do
    echo "Computing average coverage over ${annot}"
    bedtools intersect -wa -wb \
      -a ../data/annotation/annot_${annot}_hg19.bed \
      -b ./analysis/pulldown_coverages/${humanID}_pulldown_coverage.bdg \
    | awk -v OFS='\t' '{print $1, $2, $3, $4, $8}' \
    | bedtools groupby -g 1-4 -c 5 -o mean \
    > ./analysis/summary/tables/${humanID}_avg_coverage_${annot}.bed
  done
