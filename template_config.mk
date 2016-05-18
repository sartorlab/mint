
################################################################################
# Genome paths
# Genomes, like as downloaded from iGenomes, will also contain bisulfite
# converted genomes for
GENOME_PATH := /path/to/genome
# Location of bowtie2 indexes
BOWTIE2_GENOME_PATH := /path/to/bowtie2/indices
# Location of chromosome length file
CHROM_PATH := /path/to/chromosome/length/file

################################################################################
# Path to tools

PATH_TO_FASTQC := $(shell which fastqc)
PATH_TO_TRIMGALORE := $(shell which trim_galore)
PATH_TO_BISMARK := $(shell which bismark)
PATH_TO_EXTRACTOR := $(shell which bismark_methylation_extractor)
PATH_TO_BOWTIE2 := $(shell which bowtie2)
PATH_TO_R := $(shell which Rscript)
PATH_TO_MACS := $(shell which macs2)
PATH_TO_PEPR := $(shell which PePr)
PATH_TO_SAMTOOLS := $(shell which samtools)
PATH_TO_BEDTOOLS := $(shell which bedtools)
PATH_TO_BEDOPS := $(shell which bedops)
PATH_TO_AWK := $(shell which awk)
PATH_TO_BDG2BW := $(shell which bedGraphToBigWig)
PATH_TO_BDG2BB := $(shell which bedToBigBed)

################################################################################
# Command line option for minimum coverage required for
# bismark_methylation_extractor and scripts/classify_prepare_bisulfite_sample.awk
# in the sample classification module
# NOTE: This is decoupled from the minCount parameter for methylSig runs
OPT_MIN_COV = 5

################################################################################
# Command line options for methylSig and compare classifications

# DMC for CpG resolution, and DMR for region resolution (window size parameter
# used in the methylSig options below).
OPT_DM_TYPE = DMR

# Thresholds to use for DMCs or DMRs (above) in methylSig
OPT_MSIG_DM_FDR_THRESHOLD = 0.05
OPT_MSIG_DM_DIFF_THRESHOLD = 10

################################################################################
# Command line options for tools

# FastQC
OPTS_FASTQC = --format fastq --noextract
# trim_galore bisulfite
OPTS_TRIMGALORE_BISULFITE = --quality 20 --illumina --stringency 6 -e 0.2 --gzip --length 25 --rrbs
# trim_galore pulldown
OPTS_TRIMGALORE_PULLDOWN = --quality 20 --illumina --stringency 6 -e 0.2 --gzip --length 25
# bismark
OPTS_BISMARK = --bowtie2 $(GENOME_PATH)
# bismark_methylation_extractor
OPTS_EXTRACTOR = --single-end --gzip --bedGraph --cutoff $(OPT_MIN_COV) --cytosine_report --genome_folder $(GENOME_PATH) --multicore 1
# bowtie2
OPTS_BOWTIE2 = -q -x $(BOWTIE2_GENOME_PATH) -U
# macs2
OPTS_MACS = -t $bowtie2Bam -c $bowtie2InputBam -f BAM -g hs --outdir ./analysis/macs_peaks -n $macsPrefix

################################################################################
# Comparison specific options
