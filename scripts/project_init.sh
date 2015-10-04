#!/bin/bash
# Script to initialize directory structure for the project whose name is passed as the first parameters

BASEDIR=~/latte/mint/
# The desired name of the project. Directory structure will be compartmentalized
PROJECT=$1

# Go to base dir
cd $BASEDIR

# Make project specific folders in data/, analysis/, and scripts/
mkdir data/${PROJECT}
mkdir analysis/${PROJECT}
mkdir scripts/${PROJECT}

# Make subfolders in data/ and analysis/
mkdir data/${PROJECT}/raw_fastqs

mkdir analysis/${PROJECT}/raw_fastqcs
mkdir analysis/${PROJECT}/trim_fastqs
mkdir analysis/${PROJECT}/trim_fastqcs
mkdir analysis/${PROJECT}/bowtie2_bams
mkdir analysis/${PROJECT}/pulldown_coverages
mkdir analysis/${PROJECT}/bismark_bams
mkdir analysis/${PROJECT}/bismark_extractor_calls
mkdir analysis/${PROJECT}/macs_peaks
mkdir analysis/${PROJECT}/pepr_peaks
mkdir analysis/${PROJECT}/methylsig_calls
mkdir analysis/${PROJECT}/classification_simple
mkdir analysis/${PROJECT}/classification_sample
mkdir analysis/${PROJECT}/classification_comparison

mkdir analysis/${PROJECT}/summary

mkdir analysis/${PROJECT}/summary/figures
mkdir analysis/${PROJECT}/summary/tables
mkdir analysis/${PROJECT}/summary/reports

mkdir analysis/${PROJECT}/summary/ucsc_trackhub
mkdir analysis/${PROJECT}/summary/ucsc_trackhub/hg19
