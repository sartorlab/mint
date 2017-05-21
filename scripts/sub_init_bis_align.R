################################################################################
# MAKEFILE: bisulfite_align rules

if(bool_bis_samp) {

make_var_bis_align_prefix = sprintf('
################################################################################
# Workflow for bisulfite_align

BISULFITE_ALIGN_PREFIXES := %s', paste(bisulfite_samples$fullHumanID, collapse=' '))

# NOTE: This cannot be indented because they would mess up the makefile
make_rule_bis_align = '########################################

.PHONY : bisulfite_align
bisulfite_align :	 $(patsubst %,$(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_BIS_TRIM_FASTQCS)/%_trimmed_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bam,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(DIR_MULTIQC)/bisulfite_align/multiqc_report.html

########################################
.PHONY : bisulfite_raw_fastqc
bisulfite_raw_fastqc : $(patsubst %,$(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES))

$(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip :
	$(PATH_TO_FASTQC) --format fastq --noextract --outdir $(@D) $(DIR_BIS_RAW_FASTQS)/$*.fastq.gz

########################################
.PHONY : bisulfite_trim
bisulfite_trim : $(patsubst %,$(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz,$(BISULFITE_ALIGN_PREFIXES))

# Rule for trim_galore
$(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz :
	$(PATH_TO_TRIMGALORE) $(OPTS_TRIMGALORE_BISULFITE) --output_dir $(@D) $(DIR_BIS_RAW_FASTQS)/$*.fastq.gz

########################################
.PHONY : bisulfite_trim_fastqc
bisulfite_trim_fastqc : $(patsubst %,$(DIR_BIS_TRIM_FASTQCS)/%_trimmed_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES))

# Rule for FastQC on trimmed
$(DIR_BIS_TRIM_FASTQCS)/%_trimmed_fastqc.zip : $(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz
	$(PATH_TO_FASTQC) --format fastq --noextract --outdir $(@D) $<

########################################
.PHONY : bisulfite_bismark
bisulfite_bismark : $(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bam,$(BISULFITE_ALIGN_PREFIXES))

# Rule for bismark alignment
$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bam : $(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz
	$(PATH_TO_BISMARK) $(OPTS_BISMARK) --output_dir $(@D) --temp_dir $(@D) $<
	$(PATH_TO_SAMTOOLS) sort -o $@ $@
	$(PATH_TO_SAMTOOLS) index $@

########################################
# Rule to do multiqc on the bisulfite_align results
.PHONY : bisulfite_align_multiqc
bisulfite_align_multiqc : $(DIR_MULTIQC)/bisulfite_align/multiqc_report.html

$(DIR_MULTIQC)/bisulfite_align/multiqc_report.html :	 $(patsubst %,$(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_TRIM_FASTQCS)/%_trimmed_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bam,$(BISULFITE_ALIGN_PREFIXES))
	$(PATH_TO_MULTIQC) --force ./bisulfite --outdir $(@D)

################################################################################
'

cat(make_var_bis_align_prefix, file = file_make, sep = '\n', append = TRUE)
cat(make_rule_bis_align, file = file_make, sep = '\n', append = TRUE)

########################################################################
# OPTS for config.mk
config_bis_align = '
################################################################################
# bisulfite_align configuration options

# trim_galore bisulfite
# NOTE: IS YOUR DATA RRBS? If not, remove the --rrbs flag
# For trim_galore parameters see http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_User_Guide_v0.4.1.pdf
OPTS_TRIMGALORE_BISULFITE = --quality 20 --illumina --stringency 6 -e 0.2 --gzip --length 25 --rrbs

# bismark
# For bismark parameters see http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide_v0.15.0.pdf
OPTS_BISMARK = --bowtie2 $(GENOME_PATH)
'
cat(config_bis_align, file = file_config, sep='\n', append=T)

}
