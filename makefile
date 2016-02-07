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
$(PROJECT)/bis_mc_hmc/trim_fastqcs/%_trimmed_fastqc.zip : $(PROJECT)/bis_mc_hmc/trim_fastqs/%_trimmed.fq.gz
	fastqc $(OPTS_FASTQC) --outdir $(@D) $^

# Rule for bismark alignment
$(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.bam : $(PROJECT)/bis_mc_hmc/trim_fastqs/%_trimmed.fq.gz
	bismark $(OPTS_BISMARK) --output_dir $(@D) --temp_dir $(@D) $^

# Rule for bismark methylation extractor
$(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt : $(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.bam
	bismark_methylation_extractor $(OPTS_EXTRACTOR) --output $(@D) $^

# Rule for methylSig input
$(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt : $(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt
	awk -v OFS="\t" '$$4 + $$5 > 0 { print $$1 "." $$2, $$1, $$2, $$3, $$4 + $$5, ($$4 / ($$4 + $$5))*100, ($$5 / ($$4 + $$5))*100 }' $^ | sort -T . -k2,2 -k3,3n > $@

# Rule for annotatr input
$(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_annotatr.txt : $(PROJECT)/bis_mc_hmc/bismark/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt
	awk -v OFS="\t" '$$4 + $$5 > 0 { print $$1, $$2, $$2, $$1 "." $$2, $$4 + $$5, $$3, ($$4 / ($$4 + $$5))*100 }' $^ | sort -T . -k1,1 -k2,2n > $@

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
