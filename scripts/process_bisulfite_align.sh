#!/bin/bash
set -e
set -u
set -o pipefail

# Arguments
# -project The project name given to init_project.sh
# -sampleID Some SRA, GEO, Core Sample ID, etc. identifier
# -humanID Some corresponding human readble ID
PROJECT=$2
sampleID=$4
humanID=$6

# Go to the project directory
cd ~/latte/mint/${PROJECT}

# Files
rawFastq=./data/raw_fastqs/${sampleID}.fastq.gz
sampleTrimFastq=./analysis/trim_fastqs/${sampleID}_trimmed.fq.gz
humanTrimFastq=./analysis/trim_fastqs/${humanID}_trimmed.fq.gz
bismarkBamPrefix=./bismark_bams/${humanID}_trimmed.fq.gz_bismark
bismarkBam=../bismark_bams/${humanID}_trimmed.fq.gz_bismark.bam
bismarkBedgraph=${humanID}_trimmed.fq.gz_bismark.bedGraph.gz
bismarkCytReport=${humanID}_trimmed.fq.gz_bismark.CpG_report.txt
methylSigCytReport=${humanID}_trimmed.fq.gz_bismark.CpG_report_for_methylSig.txt
bismarkBigwig=../summary/${PROJECT}_hub/hg19/${humanID}_trimmed.fq.gz_bismark.bw

# FastQC raw data
fastqc --format fastq --noextract --outdir ./analysis/raw_fastqcs $rawFastq

# Trim adapter sequence and by quality with trim_galore
trim_galore --fastqc --fastqc_args "--format fastq --noextract --outdir ./analysis/trim_fastqcs" --quality 20 --illumina --stringency 6 -e 0.2 --gzip --length 20 --rrbs --output_dir ./analysis/trim_fastqs $rawFastq

# Make the sampleID to humanID transition
# Required since trim_galore does not have an output name option
mv $sampleTrimFastq $humanTrimFastq

# Bismark on trimmed reads
bismark --bowtie2 -L 50 --output_dir ./analysis/bismark_bams --temp_dir ./analysis/bismark_bams ~/latte/Homo_sapiens/ $trimFastq

# Sort and index the .bam from bismark for more efficient storage and downstream use
samtools sort ${bismarkBamPrefix}.bam $bismarkBamPrefix
samtools index ${bismarkBamPrefix}.bam

# The --bedGraph module of the methylation extractor does not support path information
# in the reference to the extractor files, so we need to be in the bismark_extractor_calls folder
cd ./analysis/bismark_extractor_calls

# Methylation Extractor bismark BAM output
# Since we are in ./analysis/bismark_extractor calls, the input bams are in ../bismark_bams/
bismark_methylation_extractor --single-end --gzip --bedGraph --cutoff 5 --cytosine_report --genome_folder ~/latte/Homo_sapiens/ $bismarkBam

# Convert the CpG report into something useful for methylSigReadData
# We are requiring the coverage to be at least 5 reads
awk -v OFS="\t" '$4 + $5 > 4 { print $1 "." $2, $1, $2, $3, $4 + $5, ($4 / ($4 + $5))*100, ($5 / ($4 + $5))*100 }' $bismarkCytReport | sort -T . -k2,2 -k3,3n > $methylSigCytReport

# Visualize methylation rates in UCSC Genome Browser (sample-wise)
# v0.14.4 of Bismark automatically gz's bedGraph and coverage files
bedGraphToBigWig <(gunzip -c $bismarkBedgraphgz | awk 'NR > 1 {print $0}' | sort -T . -k1,1 -k2,2n) ~/latte/Homo_sapiens/chromInfo_hg19.txt $bismarkBigwig
