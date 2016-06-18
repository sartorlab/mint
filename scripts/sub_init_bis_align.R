################################################################################
# MAKEFILE: bisulfite_align rules

if(bool_bis_samp) {

make_var_bis_align_prefix = sprintf('
################################################################################
# Workflow for bisulfite_align

BISULFITE_ALIGN_PREFIXES := %s', paste(bisulfite_samples$fullHumanID, collapse=' '))

make_var_bis_align_clean_tmp = 'BISFULITE_ALIGN_CLEAN_TMP := $(patsubst %,$(DIR_CLASS_SIMPLE)/%_bisulfite_simple_class_for_annotatr.txt,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_for_annotatr.txt,$(BISULFITE_ALIGN_PREFIXES))
'

# NOTE: This cannot be indented because they would mess up the makefile
make_rule_bis_align = '########################################

.PHONY : bisulfite_align
bisulfite_align : 	$(patsubst %,$(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_BIS_TRIM_FASTQCS)/%_trimmed_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bam,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph.gz,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report.txt.gz,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_RDATA)/%_trimmed_bismark_annotatr_analysis.RData,$(BISULFITE_ALIGN_PREFIXES))\\
					$(patsubst %,$(DIR_TRACK)/%_trimmed_bismark_bt2.bw,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_CLASS_SIMPLE)/%_simple_classification.bed,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_RDATA)/%_simple_class_annotatr_analysis.RData,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(patsubst %,$(DIR_TRACK)/%_simple_classification.bb,$(BISULFITE_ALIGN_PREFIXES)) \\
					$(DIR_MULTIQC)/bisulfite/multiqc_report.html

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
	$(PATH_TO_SAMTOOLS) sort $@ $(patsubst %.bam,%,$@)
	$(PATH_TO_SAMTOOLS) index $@

########################################
.PHONY : bisulfite_extractor
bisulfite_extractor : 	$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph.gz,$(BISULFITE_ALIGN_PREFIXES)) \\
						$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report.txt.gz,$(BISULFITE_ALIGN_PREFIXES)) \\
						$(patsubst %,$(DIR_RDATA)/%_trimmed_bismark_annotatr_analysis.RData,$(BISULFITE_ALIGN_PREFIXES))\\
						$(patsubst %,$(DIR_TRACK)/%_trimmed_bismark_bt2.bw,$(BISULFITE_ALIGN_PREFIXES)) \\

# Rule for bismark methylation extractor
$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph.gz $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report.txt.gz : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bam
	cd $(DIR_BIS_BISMARK); \\
	$(PATH_TO_EXTRACTOR) $(OPTS_EXTRACTOR) $(<F)

# Rule for annotatr input
.INTERMEDIATE : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_for_annotatr.txt
$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_for_annotatr.txt : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report.txt.gz
	$(PATH_TO_AWK) -f ../../scripts/extractor_to_annotatr.awk <(gunzip -c $<) | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@

# Rule for annotatr CpG universe intput
.INTERMEDIATE : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_CpGs.txt
$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_CpGs.txt : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report.txt.gz
	$(PATH_TO_AWK) -v OFS="\\t" \'{print $$1, $$2, $$2}\' <(gunzip -c $<) | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@

# Rule for annotatr of extractor results
$(DIR_RDATA)/%_trimmed_bismark_annotatr_analysis.RData : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_for_annotatr.txt $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_CpGs.txt
	$(PATH_TO_R) ../../scripts/annotatr_bisulfite.R --file $< --genome $(GENOME) --group1 NULL --group0 NULL

# Rule for temporary unzipping of extractor bedGraph
.INTERMEDIATE : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph
$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph.gz
	gunzip -c $< | $(PATH_TO_AWK) \'NR > 1 {print $$0}\' | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@

# Rules for UCSC bigWig track of extractor
$(DIR_TRACK)/%_trimmed_bismark_bt2.bw : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph
	$(PATH_TO_BDG2BW) $< $(CHROM_PATH) $@

########################################
.PHONY : bisulfite_simple_classification
bisulfite_simple_classification : 	$(patsubst %,$(DIR_CLASS_SIMPLE)/%_simple_classification.bed,$(BISULFITE_ALIGN_PREFIXES)) \\
									$(patsubst %,$(DIR_RDATA)/%_simple_class_annotatr_analysis.RData,$(BISULFITE_ALIGN_PREFIXES)) \\
									$(patsubst %,$(DIR_TRACK)/%_simple_classification.bb,$(BISULFITE_ALIGN_PREFIXES)) \\

# Simple classification for percent methylation
$(DIR_CLASS_SIMPLE)/%_bisulfite_simple_classification.bed : $(DIR_BIS_BISMARK)/%_bisulfite_trimmed_bismark_bt2.bedGraph.gz
	$(PATH_TO_R) ../../scripts/classify_simple.R --project $(PROJECT) --inFile $< --outFile $@

# Rule for annotatr input of bisulfite simple classification
.INTERMEDIATE : $(DIR_CLASS_SIMPLE)/%_bisulfite_simple_class_for_annotatr.txt
$(DIR_CLASS_SIMPLE)/%_bisulfite_simple_class_for_annotatr.txt : $(DIR_CLASS_SIMPLE)/%_bisulfite_simple_classification.bed
	$(PATH_TO_AWK) -v OFS="\\t" \'{ print $$1, $$2, $$3, $$4 }\' $< > $@

# Rule for annotatr of simple classification
$(DIR_RDATA)/%_bisulfite_simple_class_annotatr_analysis.RData : $(DIR_CLASS_SIMPLE)/%_bisulfite_simple_class_for_annotatr.txt
	$(PATH_TO_R) ../../scripts/annotatr_classification.R --file $< --genome $(GENOME) --group1 NULL --group2 NULL

# Rule for UCSC bigBed track of simple classifiation
$(DIR_TRACK)/%_bisulfite_simple_classification.bb : $(DIR_CLASS_SIMPLE)/%_bisulfite_simple_classification.bed
	$(PATH_TO_BDG2BB) $< $(CHROM_PATH) $@

########################################
# Rule to do multiqc on the bisulfite_align results
.PHONY : bisulfite_multiqc
bisulfite_multiqc : $(DIR_MULTIQC)/bisulfite/multiqc_report.html

$(DIR_MULTIQC)/bisulfite/multiqc_report.html : 	$(patsubst %,$(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_TRIM_FASTQCS)/%_trimmed_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph.gz,$(BISULFITE_ALIGN_PREFIXES))
	multiqc ./bisulfite --outdir $(@D)

########################################
# Rule to delete all temporary files from make bis_align
.PHONY : clean_bisulfite_align_tmp
clean_bisulfite_align_tmp :
	rm -f $(BISFULITE_ALIGN_CLEAN_TMP)

################################################################################
'

cat(make_var_bis_align_prefix, file = file_make, sep = '\n', append = TRUE)
cat(make_var_bis_align_clean_tmp, file = file_make, sep = '\n', append = TRUE)
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

# Command line option for minimum coverage required for bismark_methylation_extractor
# and scripts/classify_prepare_bisulfite_sample.awk in the sample classification module
# NOTE: This does not affect methylSig runs
OPT_MIN_COV = 5

# bismark_methylation_extractor
# For methylation extractor parameters see http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide_v0.15.0.pdf
OPTS_EXTRACTOR = --single-end --gzip --bedGraph --cutoff $(OPT_MIN_COV) --cytosine_report --genome_folder $(GENOME_PATH) --multicore 1
'
cat(config_bis_align, file = file_config, sep='\n', append=T)

#######################################
# PBS script
# bisulfite_align_q = c(
# 	'#!/bin/bash',
# 	'#### Begin PBS preamble',
# 	'#PBS -N bis_align',
# 	'#PBS -l nodes=1:ppn=15,walltime=24:00:00,pmem=8gb',
# 	'#PBS -A sartor_lab',
# 	'#PBS -q first',
# 	'#PBS -M rcavalca@umich.edu',
# 	'#PBS -m abe',
# 	'#PBS -j oe',
# 	'#PBS -V',
# 	'#### End PBS preamble',
# 	'# Put your job commands after this line',
# 	sprintf('cd ~/latte/mint/projects/%s/',project),
# 	'make -j 3 bisulfite_align')
# cat(bisulfite_align_q, file=sprintf('projects/%s/pbs_jobs/bisulfite_align.q', project), sep='\n')

for(i in 1:nrow(bisulfite_samples)) {
	# trackDb.txt entry for Bismark methylation calls
	trackEntry = c(
	  sprintf('track %s_pct_meth', bisulfite_samples[i,'fullHumanID']),
	  sprintf('parent %s_sample', bisulfite_samples[i,'humanID']),
	  sprintf('bigDataUrl %s_trimmed_bismark_bt2.bw', bisulfite_samples[i,'fullHumanID']),
	  sprintf('shortLabel %s_pct_meth', bisulfite_samples[i,'fullHumanID']),
	  sprintf('longLabel %s_percent_methylation', bisulfite_samples[i,'fullHumanID']),
	  'visibility full',
	  'viewLimits 0:100',
	  'type bigWig',
	  'priority 1.4',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

	# trackDb.txt entry for bisulfite simple classification results
	trackEntry = c(
	  sprintf('track %s_simple_class', bisulfite_samples[i,'fullHumanID']),
	  sprintf('parent %s_sample', bisulfite_samples[i,'humanID']),
	  sprintf('bigDataUrl %s_simple_classification.bb', bisulfite_samples[i,'fullHumanID']),
	  sprintf('shortLabel %s_simp_class', bisulfite_samples[i,'fullHumanID']),
	  sprintf('longLabel %s_simple_classification', bisulfite_samples[i,'fullHumanID']),
	  'visibility pack',
	  'itemRgb on',
	  'type bigBed 9 .',
	  'priority 1.2',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
}

}
