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

# Go to the project directory
cd ~/latte/mint/${PROJECT}

# Create appropriate file names
peprUp=${COMPARISON}__PePr_up_peaks.bed
peprDown=${COMPARISON}__PePr_down_peaks.bed
peprUpTmp=${COMPARISON}_up_tmp.bed
peprDownTmp=${COMPARISON}_down_tmp.bed
peprUcscTmp=${COMPARISON}_ucsc_tmp.bed
peprUcscSortedTmp=${COMPARISON}_ucsc_sorted_tmp.bed
peprBigbed=./analysis/summary/ucsc_trackhub/hg19/${COMPARISON}_PePr_peaks.bb

# Since PePr doesn't allow you to specify output directory, cd into the directory
# where we want the result files to live
cd ./analysis/pepr_peaks/

# PePr to call differential Methylation
# The input files should be preceded by ../bowtie2_bams/ since we are in ${PROJECT}/analysis/pepr_peaks/ from the cd above
python2.7 /home/rcavalca/.local/lib/python2.7/site-packages/PePr-1.0.8-py2.7.egg/PePr/PePr.py --input1=$INPUT1 --input2=$INPUT2 --chip1=$CHIP1 --chip2=$CHIP2 --name=$COMPARISON --file-format=bam --peaktype=sharp --diff --threshold 1e-03 --remove_artefacts

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
