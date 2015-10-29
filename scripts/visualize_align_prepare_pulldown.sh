#!/bin/bash
set -e
set -u
set -o pipefail

PROJECT=IDH2mut_v_NBM
humanID=IDH2mut_1_hmc

cd ~/latte/mint/${PROJECT}

################################################################################
# Average pulldown coverage in 5kb tiling regions
# This takes a lot of time, consider whether it is actually useful

  bedtools intersect \
    -a ../data/annotation/windows_5kb_hg19.bed \
    -b ./analysis/pulldown_coverages/${humanID}_pulldown_coverage.bdg \
    -wa -wb \
  | bedtools groupby -g 1-3 -c 7 -o mean \
  > ${humanID}_avg_coverage_5kb_windows.bed

################################################################################
# Average pulldown coverage in annotated regions

  # Find average pulldown coverage in CpG islands
  bedtools intersect \
    -a ../data/annotation/cpg_islands_hg19_ucsc.bed \
    -b ./analysis/pulldown_coverages/${humanID}_pulldown_coverage.bdg \
    -wa -wb \
  | bedtools groupby -g 1-3 -c 7 -o mean \
  > ${humanID}_avg_coverage_annot_cpg_islands.bed

  # Find average pulldown coverage in CpG shores
  bedtools intersect \
    -a ../data/annotation/cpg_shores_hg19_ucsc.bed \
    -b ./analysis/pulldown_coverages/${humanID}_pulldown_coverage.bdg \
    -wa -wb \
  | bedtools groupby -g 1-3 -c 7 -o mean \
  > ${humanID}_avg_coverage_annot_cpg_shores.bed

  # Find average pulldown coverage in promoters
  # On unique chrom, start, end, geneid because we want to preserve all of that information
  # For promoters, swap the geneid and quantity we're interested in columns so
  # that read-in to R consistently has the quantity in the 4th column
  bedtools intersect \
    -a ../data/annotation/ldef_5kb_hg19_reduced.bed \
    -b ./analysis/pulldown_coverages/${humanID}_pulldown_coverage.bdg \
    -wa -wb \
  | bedtools groupby -g 1-4 -c 8 -o mean \
  | awk -v OFS='\t' '{print $1, $2, $3, $5, $4}' \
  > ${humanID}_avg_coverage_annot_promoters.bed

  # Find average pulldown coverage in enhancers
  # Use the Andersson enhancers to start
  bedtools intersect \
    -a ../data/annotation/andersson_permissive_enhancers.bed \
    -b ./analysis/pulldown_coverages/${humanID}_pulldown_coverage.bdg \
    -wa -wb \
  | bedtools groupby -g 1-3 -c 16 -o mean \
  > ${humanID}_avg_coverage_annot_andersson_permissive_enhancers.bed
