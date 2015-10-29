#!/bin/bash
set -e
set -u
set -o pipefail

PROJECT=IDH2mut_v_NBM
humanID=IDH2mut_1_mc_hmc

cd ~/latte/mint/${PROJECT}

################################################################################
# Average % methylation and coverage in 5kb tiling regions

  # % methylation
  bedtools intersect \
    -a ../data/annotation/windows_5kb_hg19.bed \
    -b <(zcat ./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz) \
    -wa -wb \
  | bedtools groupby -g 1-3 -c 7 -o mean \
  > ${humanID}_avg_meth_5kb_windows.bed

  # Coverage
  bedtools intersect \
    -a ../data/annotation/windows_5kb_hg19.bed \
    -b <(zcat ./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz) \
    -wa -wb \
  | awk -v OFS='\t' '{print $1, $2, $3, $8 + $9}' \
  | bedtools groupby -g 1-3 -c 4 -o mean \
  > ${humanID}_avg_coverage_5kb_windows.bed

################################################################################
# Average % methylation in annotated regions

  # Find average % methylation in CpG islands
  bedtools intersect \
    -a ../data/annotation/cpg_islands_hg19_ucsc.bed \
    -b <(zcat ./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz) \
    -wa -wb \
  | bedtools groupby -g 1-3 -c 7 -o mean \
  > ${humanID}_avg_meth_annot_cpg_islands.bed

  # Find average % methylation in CpG shores
  bedtools intersect \
    -a ../data/annotation/cpg_shores_hg19_ucsc.bed \
    -b <(zcat ./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz) \
    -wa -wb \
  | bedtools groupby -g 1-3 -c 7 -o mean \
  > ${humanID}_avg_meth_annot_cpg_shores.bed

  # Find average % methylation in promoters
  # On unique chrom, start, end, geneid because we want to preserve all of that information
  # For promoters, swap the geneid and quantity we're interested in columns so
  # that read-in to R consistently has the quantity in the 4th column
  bedtools intersect \
    -a ../data/annotation/ldef_5kb_hg19_reduced.bed \
    -b <(zcat ./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz) \
    -wa -wb \
  | bedtools groupby -g 1-4 -c 8 -o mean \
  | awk -v OFS='\t' '{print $1, $2, $3, $5, $4}' \
  > ${humanID}_avg_meth_annot_promoters.bed

  # Find average % methylation in enhancers
  # Use the Andersson enhancers to start
  bedtools intersect \
    -a ../data/annotation/andersson_permissive_enhancers.bed \
    -b <(zcat ./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz) \
    -wa -wb \
  | awk -v OFS='\t' '{print $1, $2, $3, $13, $14, $15, $16}' \
  | bedtools groupby -g 1-3 -c 7 -o mean \
  > ${humanID}_avg_meth_annot_andersson_permissive_enhancers.bed

################################################################################
# Average bisulfite coverage in annotated regions

  # Find average bisulfite coverage in CpG islands
  bedtools intersect \
    -a ../data/annotation/cpg_islands_hg19_ucsc.bed \
    -b <(zcat ./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz) \
    -wa -wb \
  | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $8 + $9}' \
  | bedtools groupby -g 1-3 -c 7 -o mean \
  > ${humanID}_avg_coverage_annot_cpg_islands.bed

  # Find average bisulfite coverage in CpG shores
  bedtools intersect \
    -a ../data/annotation/cpg_shores_hg19_ucsc.bed \
    -b <(zcat ./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz) \
    -wa -wb | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $8 + $9}' \
  | bedtools groupby -g 1-3 -c 7 -o mean \
  > ${humanID}_avg_coverage_annot_cpg_shores.bed

  # Find average bisulfite coverage in promoters
  # On unique chrom, start, end, geneid because we want to preserve all of that information
  # For promoters, swap the geneid and quantity we're interested in columns so
  # that read-in to R consistently has the quantity in the 4th column
  bedtools intersect \
    -a ../data/annotation/ldef_5kb_hg19_reduced.bed \
    -b <(zcat ./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz) \
    -wa -wb \
  | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $9 + $10}' \
  | bedtools groupby -g 1-4 -c 8 -o mean \
  | awk -v OFS='\t' '{print $1, $2, $3, $5, $4}' \
  > ${humanID}_avg_coverage_annot_promoters.bed

  # Find average bisulfite coverage in enhancers
  # Use the Andersson enhancers to start
  bedtools intersect \
    -a ../data/annotation/andersson_permissive_enhancers.bed \
    -b <(zcat ./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz) \
    -wa -wb \
  | awk -v OFS='\t' '{print $1, $2, $3, $13, $14, $15, $17 + $18}' \
  | bedtools groupby -g 1-3 -c 7 -o mean \
  > ${humanID}_avg_coverage_annot_andersson_permissive_enhancers.bed
