SHELL=/bin/bash
include config.mk

################################################################################
# Bisulfite alignment
bisulfite_align : $(BIS_MC_HMC_FILES)

# Rule for FastQC on raw
$(PROJECT)/bis_mc_hmc/raw_fastqcs/%_fastqc.zip : $(PROJECT)/bis_mc_hmc/raw_fastqs/%.fastq.gz
	fastqc $(OPTS_FASTQC) --outdir $(@D) $^

# Rule for trim_galore
$(PROJECT)/bis_mc_hmc/trim_fastqs/%_trimmed.fq.gz : $(PROJECT)/bis_mc_hmc/raw_fastqs/%.fastq.gz
	trim_galore $(OPTS_TRIMGALORE) --output_dir $(@D) $^

# Rule for FastQC on trimmed
$(PROJECT)/bis_mc_hmc/trim_fastqcs/%_trimmed.fq_fastqc.zip : $(PROJECT)/bis_mc_hmc/trim_fastqs/%_trimmed.fq.gz
	fastqc $(OPTS_FASTQC) --outdir $(@D) $^

# Rule for bismark alignment
$(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.bam : $(PROJECT)/bis_mc_hmc/trim_fastqs/%_trimmed.fq.gz
	bismark $(OPTS_BISMARK) --output_dir $(@D) --temp_dir $(@D) $^
	samtools sort $@ $(patsubst %.bam,%,$@)
	samtools index $@

# Rule for bismark methylation extractor
$(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt : $(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.bam
	cd $(PROJECT)/bis_mc_hmc/bismark;\
	bismark_methylation_extractor $(OPTS_EXTRACTOR) $(^F)
$(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.bedGraph.gz : $(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt

# Rule for methylSig input
$(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt : $(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt
	awk -f scripts/to_methylSig.awk $^ | sort -T . -k2,2 -k3,3n > $@

# Rule for annotatr input
$(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_annotatr.txt : $(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt
	awk -f scripts/to_annotatr.awk $^ | sort -T . -k1,1 -k2,2n > $@

# Rules for UCSC bigWig track
.INTERMEDIATE : %_trimmed.fq.gz_bismark_bt2.bedGraph
$(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.bedGraph : $(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.bedGraph.gz
	gunzip -c $^ | awk 'NR > 1 {print $$0}' | sort -T . -k1,1 -k2,2n > $@

$(PROJECT)/$(PROJECT)_hub/$(GENOME)/%_trimmed.fq.gz_bismark_bt2.bw : $(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.bedGraph
	bedGraphToBigWig $^ $(CHROM_PATH) $@

################################################################################
# Pulldown alignment
pulldown_align : $(PULL_HMC_FILES)

# Rule for FastQC on raw
$(PROJECT)/pull_hmc/raw_fastqcs/%_fastqc.zip : $(PROJECT)/pull_hmc/raw_fastqs/%.fastq.gz
	fastqc $(OPTS_FASTQC) --outdir $(@D) $^

# Rule for trim_galore
$(PROJECT)/pull_hmc/trim_fastqs/%_trimmed.fq.gz : $(PROJECT)/pull_hmc/raw_fastqs/%.fastq.gz
	trim_galore $(OPTS_TRIMGALORE) --output_dir $(@D) $^

# Rule for FastQC on trimmed
$(PROJECT)/pull_hmc/trim_fastqcs/%_trimmed.fq_fastqc.zip : $(PROJECT)/pull_hmc/trim_fastqs/%_trimmed.fq.gz
	fastqc $(OPTS_FASTQC) --outdir $(@D) $^

# Rule for bowtie2 alignment
$(PROJECT)/pull_hmc/bowtie2_bams/%_trimmed.fq.gz_aligned.bam : $(PROJECT)/pull_hmc/trim_fastqs/%_trimmed.fq.gz
	bowtie2 $(OPTS_BOWTIE2) $^ | samtools view -bS - > $@
	samtools sort $@ $(patsubst %.bam,%,$@)
	samtools index $@

# Rule for coverage bedGraph
$(PROJECT)/pull_hmc/pulldown_coverages/%_coverage.bdg : $(PROJECT)/pull_hmc/bowtie2_bams/%_trimmed.fq.gz_aligned.bam
	bedtools genomecov -bg -g $(CHROM_PATH) -ibam $^ | sort -T . -k1,1 -k2,2n > $@

# Rule for UCSC bigWig track
$(PROJECT)/$(PROJECT)_hub/$(GENOME)/%_coverage.bw : $(PROJECT)/pull_hmc/pulldown_coverages/%_coverage.bdg
	bedGraphToBigWig $^ $(CHROM_PATH) $@

################################################################################

# Initialize a project
# $(PROJECT) comes from config.mk
.PHONY : init
init :
	mkdir $(PROJECT)
	mkdir $(PROJECT)/{bis_mc,bis_hmc,bis_mc_hmc,pull_mc,pull_hmc}
	mkdir $(PROJECT)/{bis_mc,bis_hmc,bis_mc_hmc}/{raw_fastqs,raw_fastqcs,trim_fastqs,trim_fastqcs,bismark,methylsig_calls}
	mkdir $(PROJECT)/{pull_mc,pull_hmc}/{raw_fastqs,raw_fastqcs,trim_fastqs,trim_fastqcs,bowtie2_bams,pulldown_coverages,macs_peaks,pepr_peaks}
	mkdir $(PROJECT)/{classification_simple,classification_sample,classification_comparison}
	mkdir $(PROJECT)/summary
	mkdir $(PROJECT)/summary/{figures,tables,reports}
	mkdir $(PROJECT)/$(PROJECT)_hub
	mkdir $(PROJECT)/$(PROJECT)_hub/$(GENOME)
