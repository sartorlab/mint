#!/bin/bash
set -e
set -u
set -o pipefail

# Arguments
# -project The project name given to init_project.sh
# -comparison The comparison name used in project_create_runs.R
PROJECT=$2
COMPARISON=$4

cd ~/latte/mint/${PROJECT}

################################################################################
# Distribution of PePr peak widths

  echo "Computing PePr peak widths for ${COMPARISON}"
  cat \
    <(awk -v OFS='\t' '{print $1, $2, $3, "up", $3 - $2 + 1}' ./analysis/pepr_peaks/${COMPARISON}_PePr_up_peaks.bed) \
    <(awk -v OFS='\t' '{print $1, $2, $3, "down", $3 - $2 + 1}' ./analysis/pepr_peaks/${COMPARISON}_PePr_down_peaks.bed) \
  | sort -T . -k1,1 -k2,2n \
  > ./analysis/summary/tables/${COMPARISON}_PePr_peak_widths.bed

################################################################################
# Intersection of differentially methylation regions from PePr with annotations
# For barplots showing proportion of DMRs tested in each annotation
# Can slice by DM in annotations, DM up/down in annotations

  echo "Computing intersections of ${COMPARISON} PePr peaks with annotations"
  bedtools intersect -wa -wb \
  -a <(cat \
      <(awk -v OFS='\t' '{print $1, $2, $3, "up"}' ./analysis/pepr_peaks/${COMPARISON}_PePr_up_peaks.bed) \
      <(awk -v OFS='\t' '{print $1, $2, $3, "down"}' ./analysis/pepr_peaks/${COMPARISON}_PePr_down_peaks.bed) \
    | sort -T . -k1,1 -k2,2n) \
  -b <(cat \
      ../data/annotation/annot_cpg_islands_hg19.bed \
      ../data/annotation/annot_cpg_shores_hg19.bed \
      ../data/annotation/annot_cpg_shelves_hg19.bed \
      ../data/annotation/annot_cpg_inter_hg19.bed \
      ../data/annotation/annot_promoters_hg19.bed \
      ../data/annotation/annot_enhancers_hg19.bed \
    | sort -T . -k1,1 -k2,2n) \
  > ./analysis/summary/tables/${COMPARISON}_PePr_annot.bed
