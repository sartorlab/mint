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
trimFastq=./analysis/trim_fastqs/${humanID}_trim.fastq.gz
bismarkBamPrefix=./bismark_bams/${humanID}_trim.fastq.gz_bismark
bismarkBam=../bismark_bams/${humanID}_trim.fastq.gz_bismark.bam
bismarkBedgraphgz=${humanID}_trim.fastq.gz_bismark.bedGraph.gz
bismarkBedgraph=${humanID}_trim.fastq.gz_bismark.bedGraph
bismarkSortedBedgraph=${humanID}_trim.fastq.gz_bismark_sorted.bedGraph
bismarkCytReport=${humanID}_trim.fastq.gz_bismark.CpG_report.txt
methylSigCytReport=${humanID}_trim.fastq.gz_bismark.CpG_report_for_methylSig.txt
bismarkBigwig=../summary/${PROJECT}_hub/hg19/${humanID}_trim.fastq.gz_bismark.bw

# FastQC raw data
fastqc --format fastq --noextract --outdir ./analysis/raw_fastqcs $rawFastq

# Trim reads of adapter sequence (maybe by quality with trim_galore later)
# Based on https://cutadapt.readthedocs.org/en/stable/guide.html#bisulfite-sequencing-rrbs
# it seems reasonable to add the two wildcards NN to the beginning of the adapter
cutadapt --error-rate=0.2 --adapter=NNTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGC --minimum-length=21 --overlap=6 --output=$trimFastq $rawFastq

# FastQC trimmed reads
fastqc --format fastq --noextract --outdir ./analysis/trim_fastqcs $trimFastq

# Bismark on trimmed reads
bismark --bowtie1 --bam --seedlen 50 --output_dir ./analysis/bismark_bams --temp_dir ./analysis/bismark_bams ~/latte/Homo_sapiens/ $trimFastq

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
gunzip $bismarkBedgraphgz

    # Remove first line of bedGraph
    sed -i "1d" $bismarkBedgraph
    # Sort bedGraph output from Methylation Extractor
    sort -T . -k1,1 -k2,2n $bismarkBedgraph > $bismarkSortedBedgraph
    # Convert to bigWig and replace original bedGraph with sorted version
    bedGraphToBigWig $bismarkSortedBedgraph ~/latte/Homo_sapiens/chromInfo_hg19.txt $bismarkBigwig
    mv $bismarkSortedBedgraph $bismarkBedgraph

gzip $bismarkBedgraph

    # Create custom track file and add the relevant track
    # echo '' > $customTracks
    # sed -i "1i\track type=bigWig name=${ID}_trim description=${ID}_trim db=hg19 bigDataUrl=http://www-personal.umich.edu/~rcavalca/GSE52945/$bismarkBigwig" $customTracks
