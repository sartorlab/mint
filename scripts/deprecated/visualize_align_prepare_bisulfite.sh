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
# Take means of % methylation and bisulfite coverage over each region in each annotation
# For histograms displaying distributions in different annotations
# From the perspective of the annotation regions rather than CpGs

  for annot in cpg_islands cpg_shores cpg_shelves cpg_inter promoters enhancers windows_5kb
  do
    echo "Computing average methylation over ${annot} for ${humanID}"
    bedtools intersect -wa -wb \
      -a ../data/annotation/annot_${annot}_hg19.bed \
      -b <(zcat ./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz) \
    | awk -v OFS='\t' '{print $1, $2, $3, $4, $8, $9 + $10}' \
    | tee ./analysis/summary/tables/tmp_${humanID}_${annot}.tmp \
    | bedtools groupby -g 1-4 -c 5 -o mean \
    > ./analysis/summary/tables/${humanID}_avg_methylation_${annot}.bed

    echo "Computing average coverage over ${annot} for ${humanID}"
    bedtools groupby -g 1-4 -c 6 -o mean -i ./analysis/summary/tables/tmp_${humanID}_${annot}.tmp \
    > ./analysis/summary/tables/${humanID}_avg_coverage_${annot}.bed

    rm ./analysis/summary/tables/tmp_${humanID}_${annot}.tmp
  done
