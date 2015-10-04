#!/bin/bash

# Arguments
# -project name_of_project
# -cov coverage_files
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
COV=$4
CYT=$6
SAMPLES=$8
TREATMENT=${10}
DESTRAND=${12}
MAX=${14}
MIN=${16}
FILTER=${18}
TILE=${20}
CORES=${22}
COMPARISON=${24}

# Directories
DATADIR=~/latte/mint/data/${PROJECT}
ANALYSISDIR=~/latte/mint/analysis/${PROJECT}
HUBDIR=~/latte/mint/analysis/${PROJECT}/summary/ucsc_trackhub/hg19

# Files
mSigResults=${DATADIR}/methylsig_calls/${COMPARISON}.txt
mSigTmpBedgraph=${DATADIR}/methylsig_calls/${COMPARISON}_tmp.bdg
mSigTmpSortedBedgraph=${DATADIR}/methylsig_calls/${COMPARISON}_tmp_sorted.bdg
mSigBigwig=${HUBDIR}/${COMPARISON}.bw

# Example call
# sh ~/latte/Methylation/Methylation_Code/process_bisulfite_comparison.sh -project GSE52945 -cov IDH2mut_1_errbs.fastq.gz_bismark.bismark.cov,IDH2mut_2_errbs.fastq.gz_bismark.bismark.cov,IDH2mut_3_errbs.fastq.gz_bismark.bismark.cov,NBM_1_errbs.fastq.gz_bismark.bismark.cov,NBM_2_errbs.fastq.gz_bismark.bismark.cov,NBM_3_errbs.fastq.gz_bismark.bismark.cov -cyt IDH2mut_1_errbs.fastq.gz_bismark.CpG_report.txt,IDH2mut_2_errbs.fastq.gz_bismark.CpG_report.txt,IDH2mut_3_errbs.fastq.gz_bismark.CpG_report.txt,NBM_1_errbs.fastq.gz_bismark.CpG_report.txt,NBM_2_errbs.fastq.gz_bismark.CpG_report.txt,NBM_3_errbs.fastq.gz_bismark.CpG_report.txt -samples IDH2mut_1,IDH2mut_2,IDH2mut_3,NBM_1,NBM_2,NBM_3 -treatment 1,1,1,0,0,0 -destrand TRUE -maxcount 500 -mincount 5 -filterSNPs TRUE -tile TRUE -ncores 6 -name IDH2vNBM_methylSig_regions

# This assumes process_errbs.sh has been run for all relevant samples

# Test for differentially methylated sites/regions with methylSig

    # CpG site resolution
    Rscript ~/latte/Methylation/Methylation_Code/process_bisulfite_comparison_run_methylSig.R --project=$PROJECT --covfiles=$COV --cytfiles=$CYT --sampleids=$SAMPLES --assembly=hg19 --pipeline=bismark --context=CpG --resolution=base --treatment=$TREATMENT --destranded=$DESTRAND --maxcount=$MAX --mincount=$MIN --filterSNPs=$FILTER --ncores=$CORES --quiet=FALSE --tile=$TILE --dispersion=both --minpergroup=2,2 --comparison=$COMPARISON

# Visualize methylSig differential methylation rates in UCSC Genome Browser

    # Delete erroneous region on chr9 which appears for unknown reason
    # This is particular to GSE52945
    # sed '/249239388/d' IDH2vNBM_methylSig_regions.txt > IDH2vNBM_methylSig_regions_edit.txt
    # mv IDH2vNBM_methylSig_regions_edit.txt IDH2vNBM_methylSig_regions.txt

    # The fourth column of the bedGraph is the methylation difference
    # treatment (1) - control (0), which is the 7th column in methylSig output
    # NOTE: This filters the pvalue ($5) to be less than 0.05. May want qvalue ($6) later.
    awk '$5 < 0.05 { print $1 "\t" $2 "\t" $3 "\t" $7 }' $mSigResults > $mSigTmpBedgraph
    # Sort by first two columns
    sort -T . -k1,1 -k2,2n $mSigTmpBedgraph > $mSigTmpSortedBedgraph
    # Convert to bigWig
    bedGraphToBigWig $mSigTmpSortedBedgraph ~/latte/Homo_sapiens/chromInfo_hg19.txt $mSigBigwig

    rm $mSigTmpBedgraph
    rm $mSigTmpSortedBedgraph

    # Create custom track file and add the relevant track
    # echo '' > $customTracks
    # sed -i "1i\track type=bigWig name=${COMPARISON} description=${COMPARISON} db=hg19 bigDataUrl=http://www-personal.umich.edu/~rcavalca/GSE52945/$mSigBigwig" # $customTracks

    # track type=bigWig name=IDH2_v_NBM description=DM sites IDH2_v_NBM db=hg19 bigDataUrl=http://www-personal.umich.edu/~rcavalca/GSE52945/methylSig_IDH2vNBM_sites.bw
