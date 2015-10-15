#!/bin/bash
set -e
set -u
set -o pipefail

# Arguments
# -project The project name given to init_project.sh
# -humanID Some corresponding human readble ID
PROJECT=$2
humanID=$4

# Go to the project directory
cd ~/latte/mint/${PROJECT}

# Create appropriate file names
coverage=./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov.gz
noCoverage=./analysis/classification_simple/${humanID}_simple_tmp_no.bed
lowCoverage=./analysis/classification_simple/${humanID}_simple_tmp_low.bed
medCoverage=./analysis/classification_simple/${humanID}_simple_tmp_med.bed
highCoverage=./analysis/classification_simple/${humanID}_simple_tmp_high.bed
simpleTmp=./analysis/classification_simple/${humanID}_bisulfite_simple_tmp.bed
simpleSorted=./analysis/classification_simple/${humanID}_bisulfite_classification_simple.bed
simpleBigbed=./analysis/summary/${PROJECT}_hub/hg19/${humanID}_bisulfite_classification_simple.bb

# Create color variables
no=0,0,0
low=255,153,255
med=255,0,255
high=102,0,102

# Create individual low, med, high files
# Add in something for less than 3% - 4%
zcat $coverage | awk -v OFS="\t" '$4 < 5 { print $1, $2 - 1, $3, "5mC_5hmC_none", "1000", ".", $2, $3, "0,0,0"}' > $noCoverage
zcat $coverage | awk -v OFS="\t" '$4 < 33 && $4 >=5 { print $1, $2 - 1, $3, "5mC_5hmC_low", "1000", ".", $2, $3, "255,153,255"}' > $lowCoverage
zcat $coverage | awk -v OFS="\t" '$4 >= 33 && $4 < 66 { print $1, $2 - 1, $3, "5mC_5hmC_med", "1000", ".", $2, $3, "255,0,255"}' > $medCoverage
zcat $coverage | awk -v OFS="\t" '$4 >= 66 { print $1, $2 - 1, $3, "5mC_5hmC_high", "1000", ".", $2, $3, "102,0,102"}' > $highCoverage

# Concatenate and sort
cat $noCoverage $lowCoverage $medCoverage $highCoverage > $simpleTmp
sort -T . -k1,1 -k2,2n $simpleTmp > $simpleSorted

# Cleanup
rm $noCoverage
rm $lowCoverage
rm $medCoverage
rm $highCoverage
rm $simpleTmp

# Create bigBed
bedToBigBed $simpleSorted ~/latte/Homo_sapiens/chromInfo_hg19.txt $simpleBigbed
