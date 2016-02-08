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
# bowtie2
OPTS_BOWTIE2=-q -x $(BOWTIE2_GENOME_PATH) -U
