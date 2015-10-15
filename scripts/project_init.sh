#!/bin/bash
# Script to initialize directory structure for the project whose name is passed as the first parameters

BASEDIR=~/latte/mint/
# The desired name of the project. Directory structure will be compartmentalized
PROJECT=$1

# Go to base dir and create project directory
cd ${BASEDIR}

# Make project specific folders in data/, analysis/, and scripts/
mkdir ${PROJECT}
cd ${PROJECT}

mkdir data/
mkdir analysis/
mkdir scripts/

# Make subfolders in data/ and analysis/
mkdir data/raw_fastqs
mkdir analysis/raw_fastqcs
mkdir analysis/trim_fastqs
mkdir analysis/trim_fastqcs
mkdir analysis/bowtie2_bams
mkdir analysis/pulldown_coverages
mkdir analysis/bismark_bams
mkdir analysis/bismark_extractor_calls
mkdir analysis/macs_peaks
mkdir analysis/pepr_peaks
mkdir analysis/methylsig_calls
mkdir analysis/classification_simple
mkdir analysis/classification_sample
mkdir analysis/classification_comparison

mkdir analysis/summary

mkdir analysis/summary/figures
mkdir analysis/summary/tables
mkdir analysis/summary/reports

mkdir analysis/summary/${PROJECT}_hub
mkdir analysis/summary/${PROJECT}_hub/hg19
