#!/bin/bash

# Arguments
# -project The project name given to init_project.sh
# -chipID ID for the pulldown file
# -inputID ID for input file
# -ID the ID common to chip and input
PROJECT=$2
chipID=$4
inputID=$6
humanID=$8

# Directories
DATADIR=~/latte/mint/data/${PROJECT}
ANALYSISDIR=~/latte/mint/analysis/${PROJECT}
HUBDIR=~/latte/mint/analysis/${PROJECT}/summary/ucsc_trackhub/hg19

# Create appropriate file names
bowtie2Bam=${ANALYSISDIR}/bowtie2_bams/${chipID}_pulldown_aligned.bam
bowtie2InputBam=${ANALYSISDIR}/bowtie2_bams/${inputID}_pulldown_aligned.bam
macsPrefix=${humanID}_pulldown_macs2
macsNarrowpeak=${ANALYSISDIR}/macs_peaks/${humanID}_pulldown_macs2_peaks.narrowPeak
macsNarrowpeakSorted=${ANALYSISDIR}/macs_peaks/${humanID}_sorted_tmp.narrowPeak
macsNarrowpeakCeiling=${ANALYSISDIR}/macs_peaks/${humanID}_ceiling_tmp.narrowPeak
bowtie2InputBedgraph=${ANALYSISDIR}/pulldown_coverages/${inputID}_pulldown_zero.bdg
macsBigbed=${HUBDIR}/${humanID}_pulldown_macs2_peaks_ucsc.bb

# MACS2 to call peaks
macs2 callpeak -t $bowtie2Bam -c $bowtie2InputBam -f BAM -g hs --outdir ${ANALYSISDIR}/macs_peaks -n $macsPrefix

# Remove extraneous MACS2 output to minimize footprint
# Excel output is truly unnecessary, and summits information can be recovered
# from the narrowPeak output
    rm *peaks.xls
    rm *summits.bed

# Determine region of zero input coverage for classification
bedtools genomecov -bga -ibam $bowtie2InputBam -g ~/latte/Homo_sapiens/chromInfo_hg19.txt | grep -w '0$' > $bowtie2InputBedgraph

# Visualization in UCSC Genome Browser

    # Make sure MACS2 output is sorted
    sort -T . -k1,1 -k2,2n $macsNarrowpeak > $macsNarrowpeakSorted
    # Make sure MACS2 output score field doesn't exceed 1000
    # Fifth column of narrowPeak is the score
    awk '$5 > 1000 { print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "1000" "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 } $5 <= 1000 { print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 }' $macsNarrowpeakSorted > $macsNarrowpeakCeiling
    # Convert to bigBed
    # The *.as file apparently needs to be in the same folder as the files being converted
    bedToBigBed -type=bed6+4 -as=narrowPeak.as $macsNarrowpeakCeiling ~/latte/Homo_sapiens/chromInfo_hg19.txt $macsBigbed

    # Remove tmp file and replace original narrowPeak by sorted, ceilinged version
    # NOTE: Underlying peaks do not change
    rm $macsNarrowpeakSorted
    mv $macsNarrowpeakCeiling $macsNarrowpeak

    # Add new track to the custom track file
    # sed -i "1i\track type=bigBed name=${NAME}_peaks description=${NAME}_peaks db=hg19 bigDataUrl=http://www-personal.umich.edu/~rcavalca/GSE52945/$macsBigbed" $customTracks
