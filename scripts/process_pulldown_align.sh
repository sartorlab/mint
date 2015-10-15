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
bowtie2Bam=./analysis/bowtie2_bams/${humanID}_pulldown_aligned.bam
bowtie2BamPrefix=./analysis/bowtie2_bams/${humanID}_pulldown_aligned
pulldownBedgraph=./analysis/pulldown_coverages/${humanID}_pulldown_coverage.bdg
pulldownBigwig=./analysis/summary/${PROJECT}_hub/hg19/${humanID}_pulldown_coverage.bw

# FastQC
fastqc --format fastq --noextract --outdir ./analysis/raw_fastqcs $rawFastq

# Bowtie2 to align to reference genome
bowtie2 -q -x ~/latte/Homo_sapiens/genome -U $rawFastq | samtools view -bS - > $bowtie2Bam

# Sort and index resulting bam alignment
samtools sort $bowtie2Bam $bowtie2BamPrefix
samtools index $bowtie2Bam

# Visualization in UCSC Genome Browser

    # More general version is found in GSE52945_bam_to_bw.q
    bedtools genomecov -bg -ibam $bowtie2Bam -g ~/latte/Homo_sapiens/chromInfo_hg19.txt > $pulldownBedgraph
    bedGraphToBigWig $pulldownBedgraph ~/latte/Homo_sapiens/chromInfo_hg19.txt $pulldownBigwig

    # Add new track to the custom track file
    # sed -i "1i\track type=bigWig name=${NAME}_pulldown description=${NAME}_pulldown db=hg19 bigDataUrl=http://www-personal.umich.edu/~rcavalca/GSE52945/$pulldownBigwig" $customTracks
