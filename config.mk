# Configuration for mint pipeline analyses

################################################################################
# Project and experimental information
PROJECT=test_hybrid_compare_make
# Sequencing setup
# hybrid, pulldown, or conversion (only support first two now)
SEQ_SETUP=hybrid
# Experimental design
# comparison name or 'none'
COMPARISON=IDH2mut_v_NBM
# Genome
GENOME=hg19

################################################################################
# File information
GENOME_PATH=~/latte/Homo_sapiens/
CHROM_PATH=~/latte/Homo_sapiens/chromInfo_hg19.txt

SAMPLES=IDH2mut_1_mc_hmc_bisulfite,IDH2mut_2_mc_hmc_bisulfite,NBM_1_mc_hmc_bisulfite,NBM_2_mc_hmc_bisulfite

BIS_MC_HMC_FILES := $(shell echo \
	$(PROJECT)/bis_mc_hmc/raw_fastqcs/{$(SAMPLES)}_fastqc.zip \
	$(PROJECT)/bis_mc_hmc/trim_fastqs/{$(SAMPLES)}_trimmed.fq.gz \
	$(PROJECT)/bis_mc_hmc/trim_fastqcs/{$(SAMPLES)}_trimmed.fq_fastqc.zip \
	$(PROJECT)/bis_mc_hmc/bismark/{$(SAMPLES)}_trimmed.fq.gz_bismark_bt2.bam \
	$(PROJECT)/bis_mc_hmc/bismark/{$(SAMPLES)}_trimmed.fq.gz_bismark_bt2.CpG_report.txt \
	$(PROJECT)/bis_mc_hmc/bismark/{$(SAMPLES)}_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt \
	$(PROJECT)/bis_mc_hmc/bismark/{$(SAMPLES)}_trimmed.fq.gz_bismark_bt2.CpG_report_for_annotatr.txt)

################################################################################
# Command configuration
# FastQC
OPTS_FASTQC=--format fastq --noextract
# trim_galore
OPTS_TRIMGALORE=--quality 20 --illumina --stringency 6 -e 0.2 --gzip --length 20 --rrbs
# bismark
OPTS_BISMARK=--bowtie2 $(GENOME_PATH)
# bismark_methylation_extractor
OPTS_EXTRACTOR=--single-end --gzip --bedGraph --cutoff 5 --cytosine_report --genome_folder $(GENOME_PATH) --multicore 5
