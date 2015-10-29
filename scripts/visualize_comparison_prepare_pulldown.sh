#!/bin/bash
set -e
set -u
set -o pipefail

PROJECT=IDH2mut_v_NBM
COMPARISON=IDH2mut_v_NBM_hmc

cd ~/latte/mint/${PROJECT}

################################################################################
# Distribution of PePr peak widths



################################################################################
# Distribution of DMRs in CpG annotation features

  bedtools intersect \
  -a <(cat \
      <(awk -v OFS='\t' '{print $1, $2, $3, "up"}' ./analysis/pepr_peaks/${COMPARISON}_PePr_up_peaks.bed) \
      <(awk -v OFS='\t' '{print $1, $2, $3, "down"}' ./analysis/pepr_peaks/${COMPARISON}_PePr_down_peaks.bed) \
    | sort -k1,1 -k2,2n) \
  -b <(cat \
      <(awk -v OFS='\t' '{print $1, $2, $3, "island"}' ../data/annotation/cpg_islands_hg19_ucsc.bed) \
      <(awk -v OFS='\t' '{print $1, $2, $3, "shore"}' ../data/annotation/cpg_shores_hg19_ucsc.bed) \
      <(awk -v OFS='\t' '{print $1, $2, $3, "shelf"}' ../data/annotation/cpg_shelves_hg19_ucsc.bed) \
      <(awk -v OFS='\t' '{print $1, $2, $3, "inter"}' ../data/annotation/cpg_inter_hg19_ucsc.bed) \
    | sort -T . -k1,1 -k2,2n) \
  -wa -wb \
  > ${COMPARISON}_PePr_cpg_annot.bed
