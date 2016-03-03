#!/bin/bash
set -e
set -u
set -o pipefail

# Arguments
# -project The project name given to init_project.sh
# -humanID A human readble ID
PROJECT=$2
humanID=$4

cd ~/latte/mint/${PROJECT}

################################################################################
# Distribution of MACS2 peak widths

  echo "Computing MACS2 peak widths for ${humanID}"
  awk -v OFS='\t' '{print $1, $2, $3, $3 - $2 + 1}' ./analysis/macs_peaks/${humanID}_pulldown_macs2_peaks.narrowPeak \
  > ./analysis/summary/tables/${humanID}_macs2_peak_widths.bed

################################################################################
# Intersection of peaks from MACS2 with annotations
# For barplots showing proportion of peaks in each annotation

  echo "Computing intersection of ${humanID} MACS2 peaks with annotations"
  bedtools intersect -wa -wb \
  -a ./analysis/macs_peaks/${humanID}_pulldown_macs2_peaks.narrowPeak \
  -b <(cat \
      ../data/annotation/annot_cpg_islands_hg19.bed \
      ../data/annotation/annot_cpg_shores_hg19.bed \
      ../data/annotation/annot_cpg_shelves_hg19.bed \
      ../data/annotation/annot_cpg_inter_hg19.bed \
      ../data/annotation/annot_promoters_hg19.bed \
      ../data/annotation/annot_enhancers_hg19.bed \
    | sort -T . -k1,1 -k2,2n) \
  > ./analysis/summary/tables/${humanID}_macs2_annot.bed
