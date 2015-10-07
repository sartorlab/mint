#!/bin/bash

# Arguments
# -project The project name given to init_project.sh
# -humanID Some corresponding human readble ID
PROJECT=$2
humanID=$4

# Go to the project directory
cd ~/latte/mint/${PROJECT}

# Create appropriate file names
coverage=./analysis/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bismark.cov
noCoverage=./analysis/classification_simple/${humanID}_simple_tmp_no.bed
lowCoverage=./analysis/classification_simple/${humanID}_simple_tmp_low.bed
medCoverage=./analysis/classification_simple/${humanID}_simple_tmp_med.bed
highCoverage=./analysis/classification_simple/${humanID}_simple_tmp_high.bed
simpleTmp=./analysis/classification_simple/${humanID}_errbs_simple_tmp.bed
simpleSorted=./analysis/classification_simple/${humanID}_bisulfite_classification_simple.bed
simpleBigbed=./analysis/summary/ucsc_trackhub/hg19/${humanID}_bisulfite_classification_simple.bb

# Create color variables
no=0,0,0
low=255,153,255
med=255,0,255
high=102,0,102

# Create individual low, med, high files
# Add in something for less than 3% - 4%
awk '$4 < 5 { print $1 "\t" $2 - 1 "\t" $3 "\t" "5mC_5hmC_none" "\t" "1000" "\t" "." "\t" $2 "\t" $3 "\t" "0,0,0"}' $coverage > $noCoverage
awk '$4 < 33 && $4 >=5 { print $1 "\t" $2 - 1 "\t" $3 "\t" "5mC_5hmC_low" "\t" "1000" "\t" "." "\t" $2 "\t" $3 "\t" "255,153,255"}' $coverage > $lowCoverage
awk '$4 >= 33 && $4 < 66 { print $1 "\t" $2 - 1 "\t" $3 "\t" "5mC_5hmC_med" "\t" "1000" "\t" "." "\t" $2 "\t" $3 "\t" "255,0,255"}' $coverage > $medCoverage
awk '$4 >= 66 { print $1 "\t" $2 - 1 "\t" $3 "\t" "5mC_5hmC_high" "\t" "1000" "\t" "." "\t" $2 "\t" $3 "\t" "102,0,102"}' $coverage > $highCoverage

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
