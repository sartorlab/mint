#!/bin/bash

# Arguments
# -project The project name given to init_project.sh
# -sampleID Some SRA, GEO, Core Sample ID, etc. identifier
# -humanID Some corresponding human readble ID
PROJECT=$2
sampleID=$4
humanID=$6

# Directories
DATADIR=~/latte/mint/data/${PROJECT}
ANALYSISDIR=~/latte/mint/analysis/${PROJECT}
HUBDIR=~/latte/mint/analysis/${PROJECT}/summary/ucsc_trackhub/hg19

# Files
rawFastq=${DATADIR}/raw_fastqs/${sampleID}.fastq.gz
trimFastq=${ANALYSISDIR}/trim_fastqs/${humanID}_trim.fastq.gz
bismarkBam=${ANALYSISDIR}/bismark_bams/${humanID}_trim.fastq.gz_bismark.bam
bismarkBedgraph=${ANALYSISDIR}/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark.bedGraph
bismarkSortedBedgraph=${ANALYSISDIR}/bismark_extractor_calls/${humanID}_trim.fastq.gz_bismark_sorted.bedGraph
bismarkBigwig=${HUBDIR}/${humanID}_trim.fastq.gz_bismark.bw

# FastQC raw data
fastqc --format fastq --noextract --outdir ${ANALYSISDIR}/raw_fastqcs $rawFastq

# Trim reads of adapter sequence (maybe by quality with trim_galore later)
# Based on https://cutadapt.readthedocs.org/en/stable/guide.html#bisulfite-sequencing-rrbs
# it seems reasonable to add the two wildcards NN to the beginning of the adapter
cutadapt --error-rate=0.2 --adapter=NNTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGC --minimum-length=21 --overlap=6 --output=$trimFastq $rawFastq

# FastQC trimmed reads
fastqc --format fastq --noextract --outdir ${ANALYSISDIR}/trim_fastqcs $trimFastq

# Bismark on trimmed reads
bismark --bam --seedlen 50 --output_dir ${ANALYSISDIR}/bismark_bams ~/latte/Homo_sapiens/ $trimFastq

# Methylation Extractor on trimmed reads
bismark_methylation_extractor --output ${ANALYSISDIR}/bismark_extractor_calls --single-end --bedGraph --cutoff 5 --cytosine_report --genome_folder ~/latte/Homo_sapiens/ $bismarkBam

# Visualize methylation rates in UCSC Genome Browser (sample-wise)

    # Remove first line of bedGraph
    sed -i "1d" $bismarkBedgraph
    # Sort bedGraph output from Methylation Extractor
    sort -T . -k1,1 -k2,2n $bismarkBedgraph > $bismarkSortedBedgraph
    # Convert to bigWig and replace original bedGraph with sorted version
    bedGraphToBigWig $bismarkSortedBedgraph ~/latte/Homo_sapiens/chromInfo_hg19.txt $bismarkBigwig
    mv $bismarkSortedBedgraph $bismarkBedgraph

    # Create custom track file and add the relevant track
    # echo '' > $customTracks
    # sed -i "1i\track type=bigWig name=${ID}_trim description=${ID}_trim db=hg19 bigDataUrl=http://www-personal.umich.edu/~rcavalca/GSE52945/$bismarkBigwig" $customTracks
