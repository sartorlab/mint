# Configuration for mint pipeline analyses

################################################################################
# Project and experimental information

# Project name
PROJECT = project_name
# Genome
GENOME = genome

################################################################################
# Genome paths
GENOME_PATH := ~/latte/Homo_sapiens/
BOWTIE2_GENOME_PATH := ~/latte/Homo_sapiens/genome
CHROM_PATH := ~/latte/Homo_sapiens/chromInfo_hg19.txt

################################################################################
# Path to tools

PATH_TO_FASTQC := $(shell which fastqc)
PATH_TO_TRIMGALORE := $(shell which trim_galore)
PATH_TO_BISMARK := $(shell which bismark)
PATH_TO_EXTRACTOR := $(shell which bismark_methylation_extractor)
PATH_TO_R := $(shell which Rscript)
PATH_TO_MACS := $(shell which macs2)
PATH_TO_PEPR := python2.7 /home/rcavalca/.local/lib/python2.7/site-packages/PePr-1.0.8-py2.7.egg/PePr/PePr.py
PATH_TO_SAMTOOLS := $(shell which samtools)
PATH_TO_AWK := $(shell which awk)
PATH_TO_BDG2BW := $(shell which bedGraphToBigWig)
PATH_TO_BDG2BB := $(shell which bedToBigBed)

################################################################################
# Command line options for tools

# FastQC
OPTS_FASTQC = --format fastq --noextract
# trim_galore bisulfite
OPTS_TRIMGALORE_BISULFITE = --quality 20 --illumina --stringency 6 -e 0.2 --gzip --length 20 --rrbs
# trim_galore pulldown
OPTS_TRIMGALORE_PULLDOWN = --quality 20 --illumina --stringency 6 -e 0.2 --gzip --length 20
# bismark
OPTS_BISMARK = --bowtie2 $(GENOME_PATH)
# bismark_methylation_extractor
OPTS_EXTRACTOR = --single-end --gzip --bedGraph --cutoff 5 --cytosine_report --genome_folder $(GENOME_PATH) --multicore 5
# methylSig
OPTS_METHYLSIG = --context CpG --resolution base --destranded TRUE --maxcount 500 --mincount 5 --filterSNPs TRUE --ncores 4 --quiet FALSE --tile TRUE --dispersion both --minpergroup 2,2
# bowtie2
OPTS_BOWTIE2 = -q -x $(BOWTIE2_GENOME_PATH) -U
# macs2
OPTS_MACS = -t $bowtie2Bam -c $bowtie2InputBam -f BAM -g hs --outdir ./analysis/macs_peaks -n $macsPrefix
# PePr
OPTS_PEPR = --file-format=bam --peaktype=sharp --diff --threshold 1e-03 --remove_artefacts

################################################################################
# Directories

# track hub directory
DIR_TRACK := $(PROJECT)_hub/$(GENOME)

# bisulfite_align directories
DIR_BIS_BISMARK := bisulfite/bismark
DIR_BIS_TRIM_FASTQS := bisulfite/trim_fastqs
DIR_BIS_TRIM_FASTQCS := bisulfite/trim_fastqcs
DIR_BIS_RAW_FASTQCS := bisulfite/raw_fastqcs
DIR_BIS_RAW_FASTQS := bisulfite/raw_fastqs

# bisulfite_compare directories
DIR_BIS_MSIG := bisulfite/methylsig_calls

# pulldown_align directories
DIR_PULL_COVERAGES := pulldown/pulldown_coverages
DIR_PULL_BOWTIE2 := pulldown/bowtie2_bams
DIR_PULL_TRIM_FASTQS := pulldown/trim_fastqs
DIR_PULL_TRIM_FASTQCS := pulldown/trim_fastqcs
DIR_PULL_RAW_FASTQCS := pulldown/raw_fastqcs
DIR_PULL_RAW_FASTQS := pulldown/raw_fastqs

# pulldown_sample directories
DIR_PULL_MACS := pulldown/macs2_peaks

# pulldown_compare directories
DIR_PULL_PEPR := pulldown/pepr_peaks
