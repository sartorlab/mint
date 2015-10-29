#!/bin/bash
set -e
set -u
set -o pipefail

PROJECT=IDH2mut_v_NBM
COMPARISON=IDH2mut_v_NBM

cd ~/latte/mint/${PROJECT}

################################################################################
# Average change in methylation in DMRs in annotated regions

  # For methylSig output
  # $5 = pvalue, $6 = qvalue, $7 = meth diff, $11 = meth in group 1, $12 = meth in group 0 (methylSig)
  # $8 = pvalue, $9 = qvalue, $10 = meth diff, $14 = meth in group 1, $15 = meth in group 0 (intersection + CpG island/shore)
  # $9 = pvalue, $10 = qvalue, $11 = meth diff, $15 = meth in group 1, $16 = meth in group 0 (intersection + promoter)
  # $17 = pvalue, $18 = qvalue, $19 = meth diff, $23 = meth in group 1, $24 = meth in group 0 (intersection + enhancers)

  # Find average % methylation change in CpG islands
  bedtools intersect \
    -a ../data/annotation/cpg_islands_hg19_ucsc.bed
    -b <(awk -v OFS='\t' 'NR > 1 && $5 < 0.05 { print $0 }' ./analysis/methylsig_calls/${COMPARISON}.txt) \
    -wa -wb \
  | bedtools groupby -g 1-3 -c 10 -o mean \
  > ${COMPARISON}_mc_hmc_avg_diff_meth_in_DM_annot_cpg_islands.bed

  # Find average % methylation change in CpG islands
  bedtools intersect \
    -a ../data/annotation/cpg_shores_hg19_ucsc.bed \
    -b <(awk -v OFS='\t' 'NR > 1 && $5 < 0.05 { print $0 }' ./analysis/methylsig_calls/${COMPARISON}.txt) \
    -wa -wb \
  | bedtools groupby -g 1-3 -c 10 -o mean \
  > ${COMPARISON}_mc_hmc_avg_diff_meth_in_DM_annot_cpg_shores.bed

  # Find average % methylation change in promoters
  bedtools intersect \
    -a ../data/annotation/ldef_5kb_hg19_reduced.bed \
    -b <(awk -v OFS='\t' 'NR > 1 && $5 < 0.05 { print $0 }' ./analysis/methylsig_calls/${COMPARISON}.txt) \
    -wa -wb \
  | bedtools groupby -g 1-4 -c 11 -o mean \
  | awk -v OFS='\t' '{print $1, $2, $3, $5, $4}' \
  > ${COMPARISON}_mc_hmc_avg_diff_meth_in_DM_annot_promoters.bed

  # Find average % methylation change in enhancers
  bedtools intersect \
    -a ../data/annotation/andersson_permissive_enhancers.bed \
    -b <(awk -v OFS='\t' 'NR > 1 && $5 < 0.05 { print $0 }' ./analysis/methylsig_calls/${COMPARISON}.txt) \
    -wa -wb \
  | bedtools groupby -g 1-3 -c 19 -o mean \
  > ${COMPARISON}_mc_hmc_avg_diff_meth_in_DM_annot_andersson_permissive_enhanceres.bed

################################################################################
# Distribution of DMRs in CpG annotation features

  bedtools intersect \
    -a <(awk -v OFS='\t' 'NR > 1 { print $1, $2, $3, $5, $6, $7, $11, $12 }' ./analysis/methylsig_calls/${COMPARISON}.txt) \
    -b <(cat \
        <(awk -v OFS='\t' '{print $1, $2, $3, "island"}' ../data/annotation/cpg_islands_hg19_ucsc.bed) \
        <(awk -v OFS='\t' '{print $1, $2, $3, "shore"}' ../data/annotation/cpg_shores_hg19_ucsc.bed) \
        <(awk -v OFS='\t' '{print $1, $2, $3, "shelf"}' ../data/annotation/cpg_shelves_hg19_ucsc.bed) \
        <(awk -v OFS='\t' '{print $1, $2, $3, "inter"}' ../data/annotation/cpg_inter_hg19_ucsc.bed) \
      | sort -T . -k1,1 -k2,2n) \
    -wa -wb \
  > ${COMPARISON}_mc_hmc_methylSig_cpg_annot.bed
