cd ~/latte/mint/IDH2mut_v_NBM

################################################################################
# Average % methylation and coverage in 5kb tiling regions

  # % methylation
  bedtools intersect -a ../data/annotation/windows_5kb_hg19.bed -b <(zcat ./analysis/bismark_extractor_calls/IDH2mut_1_mc_hmc_trim.fastq.gz_bismark.bismark.cov.gz) -wa -wb | bedtools groupby -g 1-3 -c 7 -o mean > IDH2mut_1_mc_hmc_avg_meth_5kb_windows.bed

  # Coverage
  bedtools intersect -a ../data/annotation/windows_5kb_hg19.bed -b <(zcat ./analysis/bismark_extractor_calls/IDH2mut_1_mc_hmc_trim.fastq.gz_bismark.bismark.cov.gz) -wa -wb | awk -v OFS='\t' '{print $1, $2, $3, $8 + $9}' | bedtools groupby -g 1-3 -c 4 -o mean > IDH2mut_1_mc_hmc_avg_coverage_5kb_windows.bed

################################################################################
# Average % methylation in annotated regions

  # Find average % methylation in CpG islands
  bedtools intersect -a ../data/annotation/cpg_islands_hg19_ucsc.bed -b <(zcat ./analysis/bismark_extractor_calls/IDH2mut_1_mc_hmc_trim.fastq.gz_bismark.bismark.cov.gz) -wa -wb | bedtools groupby -g 1-3 -c 7 -o mean > IDH2mut_1_mc_hmc_avg_meth_annot_cpg_islands.bed

  # Find average % methylation in CpG shores
  bedtools intersect -a ../data/annotation/cpg_shores_hg19_ucsc.bed -b <(zcat ./analysis/bismark_extractor_calls/IDH2mut_1_mc_hmc_trim.fastq.gz_bismark.bismark.cov.gz) -wa -wb | bedtools groupby -g 1-3 -c 7 -o mean > IDH2mut_1_mc_hmc_avg_meth_annot_cpg_shores.bed

  # Find average % methylation in promoters
  # On unique chrom, start, end, geneid because we want to preserve all of that information
  # For promoters, swap the geneid and quantity we're interested in columns so
  # that read-in to R consistently has the quantity in the 4th column
  bedtools intersect -a ../data/annotation/ldef_5kb_hg19_reduced.bed -b <(zcat ./analysis/bismark_extractor_calls/IDH2mut_1_mc_hmc_trim.fastq.gz_bismark.bismark.cov.gz) -wa -wb | bedtools groupby -g 1-4 -c 8 -o mean | awk -v OFS='\t' '{print $1, $2, $3, $5, $4}' > IDH2mut_1_mc_hmc_avg_meth_annot_promoters.bed

  # Find average % methylation in enhancers
  # Use the Andersson enhancers to start
  bedtools intersect -a ../data/annotation/andersson_permissive_enhancers.bed -b <(zcat ./analysis/bismark_extractor_calls/IDH2mut_1_mc_hmc_trim.fastq.gz_bismark.bismark.cov.gz) -wa -wb | awk -v OFS='\t' '{print $1, $2, $3, $13, $14, $15, $16}' | bedtools groupby -g 1-3 -c 7 -o mean > IDH2mut_1_mc_hmc_avg_meth_annot_andersson_permissive_enhancers.bed

################################################################################
# Average bisulfite coverage in annotated regions

  # Find average bisulfite coverage in CpG islands
  bedtools intersect -a ../data/annotation/cpg_islands_hg19_ucsc.bed -b <(zcat ./analysis/bismark_extractor_calls/IDH2mut_1_mc_hmc_trim.fastq.gz_bismark.bismark.cov.gz) -wa -wb | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $8 + $9}' | bedtools groupby -g 1-3 -c 7 -o mean > IDH2mut_1_mc_hmc_avg_coverage_annot_cpg_islands.bed

  # Find average bisulfite coverage in CpG shores
  bedtools intersect -a ../data/annotation/cpg_shores_hg19_ucsc.bed -b <(zcat ./analysis/bismark_extractor_calls/IDH2mut_1_mc_hmc_trim.fastq.gz_bismark.bismark.cov.gz) -wa -wb | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $8 + $9}' | bedtools groupby -g 1-3 -c 7 -o mean > IDH2mut_1_mc_hmc_avg_coverage_annot_cpg_shores.bed

  # Find average bisulfite coverage in promoters
  # On unique chrom, start, end, geneid because we want to preserve all of that information
  # For promoters, swap the geneid and quantity we're interested in columns so
  # that read-in to R consistently has the quantity in the 4th column
  bedtools intersect -a ../data/annotation/ldef_5kb_hg19_reduced.bed -b <(zcat ./analysis/bismark_extractor_calls/IDH2mut_1_mc_hmc_trim.fastq.gz_bismark.bismark.cov.gz) -wa -wb | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $9 + $10}' | bedtools groupby -g 1-4 -c 8 -o mean | awk -v OFS='\t' '{print $1, $2, $3, $5, $4}' > IDH2mut_1_mc_hmc_avg_coverage_annot_promoters.bed

  # Find average bisulfite coverage in enhancers
  # Use the Andersson enhancers to start
  bedtools intersect -a ../data/annotation/andersson_permissive_enhancers.bed -b <(zcat ./analysis/bismark_extractor_calls/IDH2mut_1_mc_hmc_trim.fastq.gz_bismark.bismark.cov.gz) -wa -wb | awk -v OFS='\t' '{print $1, $2, $3, $13, $14, $15, $17 + $18}' | bedtools groupby -g 1-3 -c 7 -o mean > IDH2mut_1_mc_hmc_avg_coverage_annot_andersson_permissive_enhancers.bed

################################################################################
# Average change in methylation in DMRs in annotated regions

  # For methylSig output
  # $5 = pvalue, $6 = qvalue, $7 = meth diff, $11 = meth in group 1, $12 = meth in group 0 (methylSig)
  # $8 = pvalue, $9 = qvalue, $10 = meth diff, $14 = meth in group 1, $15 = meth in group 0 (intersection + CpG island/shore)
  # $9 = pvalue, $10 = qvalue, $11 = meth diff, $15 = meth in group 1, $16 = meth in group 0 (intersection + promoter)
  # $17 = pvalue, $18 = qvalue, $19 = meth diff, $23 = meth in group 1, $24 = meth in group 0 (intersection + enhancers)

  # Find average % methylation change in CpG islands
  bedtools intersect -a ../data/annotation/cpg_islands_hg19_ucsc.bed -b <(awk -v OFS='\t' 'NR > 1 && $5 < 0.05 { print $0 }' ./analysis/methylsig_calls/IDH2mut_v_NBM.txt) -wa -wb | bedtools groupby -g 1-3 -c 10 -o mean > IDH2mut_v_NBM_mc_hmc_avg_diff_meth_in_DM_annot_cpg_islands.bed

  # Find average % methylation change in CpG islands
  bedtools intersect -a ../data/annotation/cpg_shores_hg19_ucsc.bed -b <(awk -v OFS='\t' 'NR > 1 && $5 < 0.05 { print $0 }' ./analysis/methylsig_calls/IDH2mut_v_NBM.txt) -wa -wb | bedtools groupby -g 1-3 -c 10 -o mean > IDH2mut_v_NBM_mc_hmc_avg_diff_meth_in_DM_annot_cpg_shores.bed

  # Find average % methylation change in promoters
  bedtools intersect -a ../data/annotation/ldef_5kb_hg19_reduced.bed -b <(awk -v OFS='\t' 'NR > 1 && $5 < 0.05 { print $0 }' ./analysis/methylsig_calls/IDH2mut_v_NBM.txt) -wa -wb | bedtools groupby -g 1-4 -c 11 -o mean | awk -v OFS='\t' '{print $1, $2, $3, $5, $4}' > IDH2mut_v_NBM_mc_hmc_avg_diff_meth_in_DM_annot_promoters.bed

  # Find average % methylation change in enhancers
  bedtools intersect -a ../data/annotation/andersson_permissive_enhancers.bed -b <(awk -v OFS='\t' 'NR > 1 && $5 < 0.05 { print $0 }' ./analysis/methylsig_calls/IDH2mut_v_NBM.txt) -wa -wb | bedtools groupby -g 1-3 -c 19 -o mean > IDH2mut_v_NBM_mc_hmc_avg_diff_meth_in_DM_annot_andersson_permissive_enhanceres.bed
