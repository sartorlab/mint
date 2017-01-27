################################################################################
# MAKEFILE: pulldown_align rules

if(bool_pull_samp) {

make_var_pull_align_prefix = sprintf('
################################################################################
# Workflow for pulldown_align

PULLDOWN_ALIGN_PREFIXES := %s', paste(pulldown_samples$fullHumanID, collapse=' '))

# NOTE: This cannot be indented because they would mess up the makefile
make_rule_pull_align = '
########################################

.PHONY : pulldown_align
pulldown_align :	$(patsubst %,$(DIR_PULL_RAW_FASTQCS)/%_fastqc.zip,$(PULLDOWN_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_PULL_TRIM_FASTQS)/%_trimmed.fq.gz,$(PULLDOWN_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_PULL_TRIM_FASTQCS)/%_trimmed_fastqc.zip,$(PULLDOWN_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_PULL_BOWTIE2)/%_trimmed.fq.gz_aligned.bam,$(PULLDOWN_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_PULL_COVERAGES)/%_coverage_merged.bdg,$(PULLDOWN_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_TRACK)/%_coverage.bw,$(PULLDOWN_ALIGN_PREFIXES)) \\
					$(DIR_MULTIQC)/pulldown/multiqc_report.html

########################################
.PHONY : pulldown_raw_fastqc
pulldown_raw_fastqc : $(patsubst %,$(DIR_PULL_RAW_FASTQCS)/%_fastqc.zip,$(PULLDOWN_ALIGN_PREFIXES))

# Rule for FastQC on raw
$(DIR_PULL_RAW_FASTQCS)/%_fastqc.zip :
	$(PATH_TO_FASTQC) --format fastq --noextract --outdir $(@D) $(DIR_PULL_RAW_FASTQS)/$*.fastq.gz

########################################
.PHONY : pulldown_trim
pulldown_trim : $(patsubst %,$(DIR_PULL_TRIM_FASTQS)/%_trimmed.fq.gz,$(PULLDOWN_ALIGN_PREFIXES))

# Rule for trim_galore
$(DIR_PULL_TRIM_FASTQS)/%_trimmed.fq.gz :
	$(PATH_TO_TRIMGALORE) $(OPTS_TRIMGALORE_PULLDOWN) --output_dir $(@D) $(DIR_PULL_RAW_FASTQS)/$*.fastq.gz

########################################
.PHONY : pulldown_trim_fastqc
pulldown_trim_fastqc : $(patsubst %,$(DIR_PULL_TRIM_FASTQCS)/%_trimmed_fastqc.zip,$(PULLDOWN_ALIGN_PREFIXES))

# Rule for FastQC on trimmed
$(DIR_PULL_TRIM_FASTQCS)/%_trimmed_fastqc.zip : $(DIR_PULL_TRIM_FASTQS)/%_trimmed.fq.gz
	$(PATH_TO_FASTQC) --format fastq --noextract --outdir $(@D) $<

########################################
.PHONY : pulldown_bowtie2
pulldown_bowtie2 : $(patsubst %,$(DIR_PULL_BOWTIE2)/%_trimmed.fq.gz_aligned.bam,$(PULLDOWN_ALIGN_PREFIXES))

# Rule for bowtie2 alignment
$(DIR_PULL_BOWTIE2)/%_trimmed.fq.gz_aligned.bam : $(DIR_PULL_TRIM_FASTQS)/%_trimmed.fq.gz $(DIR_PULL_TRIM_FASTQCS)/%_trimmed_fastqc.zip
	$(PATH_TO_BOWTIE2) $(OPTS_BOWTIE2) $< 2> $(@D)/$(@F).txt | $(PATH_TO_SAMTOOLS) view -bS - > $@
	$(PATH_TO_SAMTOOLS) sort -o $@ $@
	$(PATH_TO_SAMTOOLS) index $@

########################################
.PHONY : pulldown_coverage
pulldown_coverage : $(patsubst %,$(DIR_PULL_COVERAGES)/%_coverage_merged.bdg,$(PULLDOWN_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_TRACK)/%_coverage.bw,$(PULLDOWN_ALIGN_PREFIXES))

# Rule for coverage bedGraph
.INTERMEDIATE : $(DIR_PULL_COVERAGES)/%_coverage.bdg
$(DIR_PULL_COVERAGES)/%_coverage.bdg : $(DIR_PULL_BOWTIE2)/%_trimmed.fq.gz_aligned.bam
	$(PATH_TO_BEDTOOLS) genomecov -bg -g $(CHROM_PATH) -ibam $< | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@

# Rule for UCSC bigWig track
$(DIR_TRACK)/%_coverage.bw : $(DIR_PULL_COVERAGES)/%_coverage.bdg
	$(PATH_TO_BDG2BW) $< $(CHROM_PATH) $@

# Rule for merged coverage BED
# For use in signal BEDs downstream
$(DIR_PULL_COVERAGES)/%_coverage_merged.bdg : $(DIR_PULL_COVERAGES)/%_coverage.bdg
	$(PATH_TO_BEDTOOLS) merge -d 20 -i $< | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@

########################################
# Rule to do multiqc on the pulldown_align results
.PHONY : pulldown_multiqc
pulldown_multiqc : $(DIR_MULTIQC)/pulldown/multiqc_report.html

$(DIR_MULTIQC)/pulldown/multiqc_report.html :	 $(patsubst %,$(DIR_PULL_RAW_FASTQCS)/%_fastqc.zip,$(PULLDOWN_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_PULL_TRIM_FASTQCS)/%_trimmed_fastqc.zip,$(PULLDOWN_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_PULL_BOWTIE2)/%_trimmed.fq.gz_aligned.bam,$(PULLDOWN_ALIGN_PREFIXES))
	$(PATH_TO_MULTIQC) --force ./pulldown --outdir $(@D)

################################################################################
'

cat(make_var_pull_align_prefix, file = file_make, sep = '\n', append = TRUE)
cat(make_rule_pull_align, file = file_make, sep = '\n', append = TRUE)

########################################################################
# OPTS for config.mk
config_pull_align = '################################################################################
# pulldown_align configuration options

# trim_galore pulldown
# For trim_galore parameters see http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_User_Guide_v0.4.1.pdf
OPTS_TRIMGALORE_PULLDOWN = --quality 20 --illumina --stringency 6 -e 0.2 --gzip --length 25

# bowtie2
# For bowtie2 parameters see http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
OPTS_BOWTIE2 = -q --no-unal -x $(BOWTIE2_GENOME_PATH) -U
'
cat(config_pull_align, file = file_config, sep='\n', append=T)


#######################################
# PBS script
# pulldown_align_q = c(
#	 '#!/bin/bash',
#	 '#### Begin PBS preamble',
#	 '#PBS -N pull_align',
#	 '#PBS -l nodes=1:ppn=8,walltime=24:00:00,pmem=16gb',
#	 '#PBS -A sartor_lab',
#	 '#PBS -q first',
#	 '#PBS -M rcavalca@umich.edu',
#	 '#PBS -m abe',
#	 '#PBS -j oe',
#	 '#PBS -V',
#	 '#### End PBS preamble',
#	 '# Put your job commands after this line',
#	 sprintf('cd ~/latte/mint/projects/%s/',project),
#	 'make -j 8 pulldown_align')
# cat(pulldown_align_q, file=sprintf('projects/%s/pbs_jobs/pulldown_align.q', project), sep='\n')

for(i in 1:nrow(pulldown_samples)) {
	# trackDb.txt entry for chip/input pulldown coverages
	trackEntry = c(
		sprintf('track %s_cov', pulldown_samples[i,'fullHumanID']),
		sprintf('parent %s_sample', pulldown_samples[i,'humanID']),
		sprintf('bigDataUrl %s_coverage.bw', pulldown_samples[i,'fullHumanID']),
		sprintf('shortLabel %s_cov', pulldown_samples[i,'fullHumanID']),
		sprintf('longLabel %s_coverage', pulldown_samples[i,'fullHumanID']),
		'visibility full',
		'autoScale on',
		'type bigWig',
		'priority 1.6',
		' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

}

}
