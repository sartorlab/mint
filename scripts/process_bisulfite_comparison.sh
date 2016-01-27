#!/bin/bash
set -e
set -u
set -o pipefail

# Arguments
# -project name_of_project
# -cyt cytosine_report_files
# -samples sample_names
# -treatment treatment_vector
# -destrand T/F
# -maxCount max
# -minCount min
# -filterSNPs T/F
# -cores number_of_cores
# -comparison name_of_comparison
PROJECT=$2
CYT=$4
SAMPLES=$6
TREATMENT=${8}
DESTRAND=${10}
MAX=${12}
MIN=${14}
FILTER=${16}
TILE=${18}
GROUP=${20}
CORES=${22}
COMPARISON=${24}

# Go to the project directory
cd ~/latte/mint/${PROJECT}

# Directories
DATADIR=~/latte/mint/data/${PROJECT}
ANALYSISDIR=~/latte/mint/analysis/${PROJECT}
SCRIPTDIR=~/latte/mint/scripts
HUBDIR=~/latte/mint/analysis/${PROJECT}/summary/${PROJECT}_hub/hg19

# Files
mSigResults=./analysis/methylsig_calls/${COMPARISON}.txt
mSigTmpResults=./analysis/methylsig_calls/${COMPARISON}_tmp.txt
mSigBigwig=./analysis/summary/${PROJECT}_hub/hg19/${COMPARISON}_methylSig.bw

# Example call
# sh ~/latte/Methylation/Methylation_Code/process_bisulfite_comparison.sh -project GSE52945 -cyt ./analysis/bismark_extractor_calls/IDH2mut_1_errbs.fastq.gz_bismark.CpG_report.txt,./analysis/bismark_extractor_calls/IDH2mut_2_errbs.fastq.gz_bismark.CpG_report.txt,./analysis/bismark_extractor_calls/IDH2mut_3_errbs.fastq.gz_bismark.CpG_report.txt,./analysis/bismark_extractor_calls/NBM_1_errbs.fastq.gz_bismark.CpG_report.txt,./analysis/bismark_extractor_calls/NBM_2_errbs.fastq.gz_bismark.CpG_report.txt,./analysis/bismark_extractor_calls/NBM_3_errbs.fastq.gz_bismark.CpG_report.txt -samples IDH2mut_1,IDH2mut_2,IDH2mut_3,NBM_1,NBM_2,NBM_3 -treatment 1,1,1,0,0,0 -destrand TRUE -maxcount 500 -mincount 5 -filterSNPs TRUE -tile TRUE -ncores 6 -name IDH2vNBM_methylSig_regions

# This assumes process_errbs.sh has been run for all relevant samples

# Test for differentially methylated sites/regions with methylSig

    # CpG site resolution
    # The cytfiles should be preceded by ./analysis/bismark_extractor_calls/
    # getwd() = ~/latte/mint/${PROJECT} since that's where we are when we call Rscript
    Rscript ~/latte/mint/scripts/process_bisulfite_comparison_run_methylSig.R --project=$PROJECT --cytfiles=$CYT --sampleids=$SAMPLES --assembly=hg19 --pipeline=bismark --context=CpG --resolution=base --treatment=$TREATMENT --destranded=$DESTRAND --maxcount=$MAX --mincount=$MIN --filterSNPs=$FILTER --ncores=$CORES --quiet=FALSE --tile=$TILE --dispersion=both --minpergroup=$GROUP --comparison=$COMPARISON

# Visualize methylSig differential methylation rates in UCSC Genome Browser

    # Delete erroneous region on chr9 which appears for unknown reason
    # This is particular to GSE52945
    # sed '/249239388/d' IDH2vNBM_methylSig_regions.txt > IDH2vNBM_methylSig_regions_edit.txt
    # mv IDH2vNBM_methylSig_regions_edit.txt IDH2vNBM_methylSig_regions.txt

    # The fourth column of the bedGraph is the methylation difference
    # treatment (1) - control (0), which is the 7th column in methylSig output
    # NOTE: This filters the pvalue ($5) to be less than 0.05. May want qvalue ($6) later.
    # Convert to bigWig
    awk -v OFS='\t' '$5 < 0.05 {print $1, $2, $3, $7 }' $mSigResults | sort -T . -k1,1 -k2,2n > $mSigTmpResults
    bedGraphToBigWig $mSigTmpResults ~/latte/Homo_sapiens/chromInfo_hg19.txt.gz $mSigBigwig
    rm $mSigTmpResults
