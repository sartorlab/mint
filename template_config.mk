
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
