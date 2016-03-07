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
# Take the mean of group methylation differences over each region in each annotation
# For histograms displaying distributions in different annotations
# From the perspective of the annotation regions rather than DMC/DMRs
# NOTE: The difference here is the average taken over each annotation region

  # methylSig columns
  # $5 = pvalue, $6 = qvalue, $7 = meth diff, $11 = meth in group 1, $12 = meth in group 0 (methylSig)

  # for annot in cpg_islands cpg_shores cpg_shelves cpg_inter promoters enhancers windows_5kb
  # do
  #   echo "Computing average differential methylation in DMs over ${annot} for ${COMPARISON}"
  #   bedtools intersect -wa -wb \
  #     -a ../data/annotation/annot_${annot}_hg19.bed \
  #     -b <(awk -v OFS='\t' 'NR > 1 && $5 < 0.05 { print $1, $2, $3, $5, $6, $7, $11, $12 }' ./analysis/methylsig_calls/${COMPARISON}.txt) \
  #   | bedtools groupby -g 1-4 -c 10 -o mean \
  #   > ./analysis/summary/tables/${COMPARISON}_methylSig_avg_diff_meth_in_DM_${annot}.bed
  # done

################################################################################
# Intersection of sites/regions tested in methylSig with annotations
# For barplots showing proportion of sites/regions tested in each annotation
# Can slice by all_tested in annotations, DM in annotations, DM up/down in annotations
# NOTE: CpG annotations are a partition of the genome. Promoter/enhancer is not.

  echo "Computing intersections with ${COMPARISON} methylSig results with annotations"
  bedtools intersect -wa -wb \
    -a <(awk -v OFS='\t' 'NR > 1 { print $1, $2, $3, $5, $6, $7, $11, $12 }' ./analysis/methylsig_calls/${COMPARISON}.txt) \
    -b <(cat \
        ../data/annotation/annot_cpg_islands_hg19.bed \
        ../data/annotation/annot_cpg_shores_hg19.bed \
        ../data/annotation/annot_cpg_shelves_hg19.bed \
        ../data/annotation/annot_cpg_inter_hg19.bed \
        ../data/annotation/annot_promoters_hg19.bed \
        ../data/annotation/annot_enhancers_hg19.bed \
      | sort -T . -k1,1 -k2,2n) \
  > ./analysis/summary/tables/${COMPARISON}_methylSig_annot.bed
