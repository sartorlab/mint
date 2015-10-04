#!/bin/bash

# Argument order
# -project The project name given to init_project.sh
# -input1 input_for_condition_1
# -input2 input_for_condition_2
# -c1 chip_for_condition_1
# -c2 chip_for_condition_2
# -n results_name
PROJECT=$2
INPUT1=$4
INPUT2=$6
CHIP1=$8
CHIP2=${10}
COMPARISON=${12}

# Directories
DATADIR=~/latte/mint/data/${PROJECT}
ANALYSISDIR=~/latte/mint/analysis/${PROJECT}
HUBDIR=~/latte/mint/analysis/${PROJECT}/summary/ucsc_trackhub/hg19

# Create appropriate file names
peprUp=${ANALYSISDIR}/pepr_peaks/${COMPARISON}__PePr_up_peaks.bed
peprDown=${ANALYSISDIR}/pepr_peaks/${COMPARISON}__PePr_down_peaks.bed
peprUpTmp=${ANALYSISDIR}/pepr_peaks/${COMPARISON}_up_tmp.bed
peprDownTmp=${ANALYSISDIR}/pepr_peaks/${COMPARISON}_down_tmp.bed
peprUcscTmp=${ANALYSISDIR}/pepr_peaks/${COMPARISON}_ucsc_tmp.bed
peprUcscSortedTmp=${ANALYSISDIR}/pepr_peaks/${COMPARISON}_ucsc_sorted_tmp.bed
peprBigbed=${HUBDIR}/${COMPARISON}_PePr_peaks_ucsc.bb

# PePr to call differential Methylation
python2.7 /home/rcavalca/.local/lib/python2.7/site-packages/PePr-1.0.8-py2.7.egg/PePr/PePr.py -input1=$INPUT1 --input2=$INPUT2 --chip1=$CHIP1 --chip2=$CHIP2 --name=$COMPARISON --file-format=bam --peaktype=sharp --diff --threshold 1e-03 --remove_artefacts --narrow_peak_width

# NOTE: Should include set differencing of PePr up/down peaks here

# Visualization in UCSC Genome Browser

    # This is BED format with the following \t separated columns
    # chrom, start, end, name, score, strand, thickStart, thickEnd, color
    awk '{ print $1 "\t" $2 "\t" $3 "\t" "up"$1":"$2 "\t" "1000" "\t" "." "\t" $2 "\t" $3 "\t" "0,0,255" }' $peprUp > $peprUpTmp
    awk '{ print $1 "\t" $2 "\t" $3 "\t" "down"$1":"$2 "\t" "1000" "\t" "." "\t" $2 "\t" $3 "\t" "102,102,255" }' $peprDown > $peprDownTmp
    cat $peprUpTmp $peprDownTmp > $peprUcscTmp

    rm $peprUpTmp
    rm $peprDownTmp

    sort -T . -k1,1 -k2,2n $peprUcscTmp > $peprUcscSortedTmp
    rm $peprUcscTmp
    bedToBigBed $peprUcscSortedTmp ~/latte/Homo_sapiens/chromInfo_hg19.txt $peprBigbed
    rm $peprUcscSortedTmp

    # Add new track to the custom track file
    # sed -i "1i\track type=narrowPeak name=${COMPARISON}_PePr description=${COMPARISON}_PePr db=hg19 bigDataUrl=http://www-personal.umich.edu/~rcavalca/GSE52945/$peprBigbed" $customTracks
