
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
