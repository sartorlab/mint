# Configuration for mint pipeline analyses

################################################################################
# Project and experimental information

PROJECT=test_hybrid_compare_make
# Sequencing setup
# hybrid, pulldown, or bisulfite (only support first two now)
SEQ_SETUP=hybrid
# Experimental design
# comparison name or 'none'
COMPARISON=IDH2mut_v_NBM
# Genome
GENOME=hg19
# Samples (comma-separated, no spaces)
SAMPLES=IDH2mut_1,IDH2mut_2,NBM_1,NBM_2
GROUP1=IDH2mut_1,IDH2mut_2
GROUP2=NBM_1,NBM_2

################################################################################
# File information

# Special variables for PePr pattern substitution
comma:= ,
empty:=
space:= $(empty) $(empty)

GENOME_PATH=~/latte/Homo_sapiens/
BOWTIE2_GENOME_PATH=~/latte/Homo_sapiens/genome
CHROM_PATH=~/latte/Homo_sapiens/chromInfo_hg19.txt

BIS_SAMPLES=IDH2mut_1_mc_hmc_bisulfite,IDH2mut_2_mc_hmc_bisulfite,NBM_1_mc_hmc_bisulfite,NBM_2_mc_hmc_bisulfite

BIS_MC_HMC_FILES := $(shell echo \
	$(PROJECT)/bis_mc_hmc/raw_fastqcs/{$(BIS_SAMPLES)}_fastqc.zip \
	$(PROJECT)/bis_mc_hmc/trim_fastqs/{$(BIS_SAMPLES)}_trimmed.fq.gz \
	$(PROJECT)/bis_mc_hmc/trim_fastqcs/{$(BIS_SAMPLES)}_trimmed.fq_fastqc.zip \
	$(PROJECT)/bis_mc_hmc/bismark/{$(BIS_SAMPLES)}_trimmed.fq.gz_bismark_bt2.bam \
	$(PROJECT)/bis_mc_hmc/bismark/{$(BIS_SAMPLES)}_trimmed.fq.gz_bismark_bt2.CpG_report.txt \
	$(PROJECT)/bis_mc_hmc/bismark/{$(BIS_SAMPLES)}_trimmed.fq.gz_bismark_bt2.bedGraph.gz \
	$(PROJECT)/bis_mc_hmc/bismark/{$(BIS_SAMPLES)}_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt \
	$(PROJECT)/bis_mc_hmc/bismark/{$(BIS_SAMPLES)}_trimmed.fq.gz_bismark_bt2.CpG_report_for_annotatr.txt \
	$(PROJECT)/$(PROJECT)_hub/$(GENOME)/{$(BIS_SAMPLES)}_trimmed.fq.gz_bismark_bt2.bw)

PULL_SAMPLES=IDH2mut_1_hmc_pulldown,IDH2mut_2_hmc_pulldown,NBM_1_hmc_pulldown,NBM_2_hmc_pulldown,IDH2mut_1_hmc_input_pulldown,IDH2mut_2_hmc_input_pulldown,NBM_1_hmc_input_pulldown,NBM_2_hmc_input_pulldown

PULL_HMC_FILES := $(shell echo \
	$(PROJECT)/pull_hmc/raw_fastqcs/{$(PULL_SAMPLES)}_fastqc.zip \
	$(PROJECT)/pull_hmc/trim_fastqs/{$(PULL_SAMPLES)}_trimmed.fq.gz \
	$(PROJECT)/pull_hmc/trim_fastqcs/{$(PULL_SAMPLES)}_trimmed.fq_fastqc.zip \
	$(PROJECT)/pull_hmc/bowtie2_bams/{$(PULL_SAMPLES)}_trimmed.fq.gz_aligned.bam \
	$(PROJECT)/pull_hmc/pulldown_coverages/{$(PULL_SAMPLES)}_coverage.bdg \
	$(PROJECT)/$(PROJECT)_hub/$(GENOME)/{$(PULL_SAMPLES)}_coverage.bw)

PULL_SAMPLES_MACS=IDH2mut_1_hmc,IDH2mut_2_hmc,NBM_1_hmc,NBM_2_hmc

PULL_SAMPLE_FILES := $(shell echo \
	$(PROJECT)/pull_hmc/macs_peaks/{$(PULL_SAMPLES_MACS)}_pulldown_macs2_peaks.narrowPeak \
	$(PROJECT)/pull_hmc/pulldown_coverages/{$(PULL_SAMPLES_MACS)}_pulldown_zero.bdg \
	$(PROJECT)/$(PROJECT)_hub/$(GENOME)/{$(PULL_SAMPLES_MACS)}_pulldown_macs2_peaks.bb)

PULL_COMPARE_GROUP1_CHIP := $(shell echo $(PROJECT)/pull_hmc/bowtie2_bams/{$(GROUP1)}_hmc_pulldown_trimmed.fq.gz_aligned.bam)
PULL_COMPARE_GROUP1_INPUT := $(shell echo $(PROJECT)/pull_hmc/bowtie2_bams/{$(GROUP1)}_hmc_input_pulldown_trimmed.fq.gz_aligned.bam)
PULL_COMPARE_GROUP2_CHIP := $(shell echo $(PROJECT)/pull_hmc/bowtie2_bams/{$(GROUP2)}_hmc_pulldown_trimmed.fq.gz_aligned.bam)
PULL_COMPARE_GROUP2_INPUT := $(shell echo $(PROJECT)/pull_hmc/bowtie2_bams/{$(GROUP2)}_hmc_input_pulldown_trimmed.fq.gz_aligned.bam)

PULL_COMPARE_FILES := $(shell echo \
	$(PROJECT)/pull_hmc/pepr_peaks/$(COMPARISON)__PePr_up_peaks.bed \
	$(PROJECT)/$(PROJECT)_hub/$(GENOME)/$(COMPARISON)_PePr_peaks.bb)

################################################################################
# Path to tools

PATH_TO_FASTQC :=$(shell which fastqc)
PATH_TO_TRIMGALORE :=$(shell which trim_galore)
PATH_TO_BISMARK :=$(shell which bismark)
PATH_TO_PEPR :=python2.7 /home/rcavalca/.local/lib/python2.7/site-packages/PePr-1.0.8-py2.7.egg/PePr/PePr.py

################################################################################
# Command line options for tools

# FastQC
OPTS_FASTQC = --format fastq --noextract
# trim_galore
OPTS_TRIMGALORE = --quality 20 --illumina --stringency 6 -e 0.2 --gzip --length 20 --rrbs
# bismark
OPTS_BISMARK = --bowtie2 $(GENOME_PATH)
# bismark_methylation_extractor
OPTS_EXTRACTOR = --single-end --gzip --bedGraph --cutoff 5 --cytosine_report --genome_folder $(GENOME_PATH) --multicore 5
# bowtie2
OPTS_BOWTIE2 = -q -x $(BOWTIE2_GENOME_PATH) -U
# macs2
OPTS_MACS = -t $bowtie2Bam -c $bowtie2InputBam -f BAM -g hs --outdir ./analysis/macs_peaks -n $macsPrefix
# PePr
OPTS_PEPR = --file-format=bam --peaktype=sharp --diff --threshold 1e-03 --remove_artefacts
