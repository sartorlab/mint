################################################################################
# MAKEFILE: bisulfite_sample rules

if(bool_bis_samp) {

make_var_bis_sample_prefix = sprintf('
################################################################################
# Workflow for bisulfite_sample

BISULFITE_SAMPLE_PREFIXES := %s', paste(bisulfite_samples$fullHumanID, collapse=' '))

make_var_bis_sample_clean_tmp = 'BISFULITE_SAMPLE_CLEAN_TMP := $(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph,$(BISULFITE_SAMPLE_PREFIXES))
'

# NOTE: This cannot be indented because they would mess up the makefile
make_rule_bis_sample = '########################################

.PHONY : bisulfite_sample
bisulfite_sample :	$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph.gz,$(BISULFITE_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report.txt.gz,$(BISULFITE_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_RDATA)/%_trimmed_bismark_annotatr_analysis.RData,$(BISULFITE_SAMPLE_PREFIXES))\\
					$(patsubst %,$(DIR_TRACK)/%_trimmed_bismark_bt2.bw,$(BISULFITE_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_CLASS_SIMPLE)/%_bismark_simple_classification.bed,$(BISULFITE_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_RDATA)/%_bismark_simple_classification_annotatr_analysis.RData,$(BISULFITE_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_TRACK)/%_bismark_simple_classification.bb,$(BISULFITE_SAMPLE_PREFIXES)) \\
					$(DIR_MULTIQC)/bisulfite_sample/multiqc_report.html

########################################
.PHONY : bisulfite_extractor
bisulfite_extractor :	 $(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph.gz,$(BISULFITE_SAMPLE_PREFIXES)) \\
						$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bismark.cov.gz,$(BISULFITE_SAMPLE_PREFIXES)) \\
						$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report.txt.gz,$(BISULFITE_SAMPLE_PREFIXES)) \\
						$(patsubst %,$(DIR_RDATA)/%_trimmed_bismark_annotatr_analysis.RData,$(BISULFITE_SAMPLE_PREFIXES))\\
						$(patsubst %,$(DIR_TRACK)/%_trimmed_bismark_bt2.bw,$(BISULFITE_SAMPLE_PREFIXES)) \\

# Rule for bismark methylation extractor
$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bismark.cov.gz $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph.gz $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report.txt.gz : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bam
	cd $(DIR_BIS_BISMARK); \\
	$(PATH_TO_EXTRACTOR) $(OPTS_EXTRACTOR) $(<F)

# Rule for temporary extractor results for annotatr
# NOTE: bismark.cov.gz file is 1-based start and end, need to subtract
# 1 from start in order for annotatr to properly interpret it
.INTERMEDIATE : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bismark.cov
$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bismark.cov : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bismark.cov.gz
	gunzip -c $< | $(PATH_TO_AWK) -v OFS="\\t" \'{print $$1, $$2 - 1, $$3, ".", $$4, ".", $$5 + $$6}\' > $@

# Rule for annotatr of extractor results
$(DIR_RDATA)/%_trimmed_bismark_annotatr_analysis.RData : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bismark.cov
	$(PATH_TO_R) ../../scripts/annotatr_annotations.R --file $< --genome $(GENOME) --annot_type bismark --group1 NULL --group0 NULL

# Rule for temporary unzipping of extractor bedGraph
.INTERMEDIATE : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph
$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph.gz
	gunzip -c $< | $(PATH_TO_AWK) \'NR > 1 {print $$0}\' | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@

# Rules for UCSC bigWig track of extractor
$(DIR_TRACK)/%_trimmed_bismark_bt2.bw : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph
	$(PATH_TO_BDG2BW) $< $(CHROM_PATH) $@

########################################
.PHONY : bisulfite_simple_classification
bisulfite_simple_classification :	 $(patsubst %,$(DIR_CLASS_SIMPLE)/%_bisulfite_bismark_simple_classification.bed,$(BISULFITE_SAMPLE_PREFIXES)) \\
									$(patsubst %,$(DIR_RDATA)/%_bismark_simple_classification_annotatr_analysis.RData,$(BISULFITE_SAMPLE_PREFIXES)) \\
									$(patsubst %,$(DIR_TRACK)/%_bismark_simple_classification.bb,$(BISULFITE_SAMPLE_PREFIXES)) \\

# Simple classification for percent methylation
$(DIR_CLASS_SIMPLE)/%_bisulfite_bismark_simple_classification.bed : $(DIR_BIS_BISMARK)/%_bisulfite_trimmed_bismark_bt2.bedGraph.gz
	$(PATH_TO_R) ../../scripts/classify_simple.R --project $(PROJECT) --inFile $< --outFile $@ --group1 NULL --group0 NULL

# Rule for annotatr of simple classification
$(DIR_RDATA)/%_bisulfite_bismark_simple_classification_annotatr_analysis.RData : $(DIR_CLASS_SIMPLE)/%_bisulfite_bismark_simple_classification.bed
	$(PATH_TO_R) ../../scripts/annotatr_annotations.R --file $< --genome $(GENOME) --annot_type simple_bisulfite_bismark --group1 NULL --group0 NULL

# Rule for UCSC bigBed track of simple classifiation
$(DIR_TRACK)/%_bisulfite_bismark_simple_classification.bb : $(DIR_CLASS_SIMPLE)/%_bisulfite_bismark_simple_classification.bed
	$(PATH_TO_BED2BB) $< $(CHROM_PATH) $@

########################################
# Rule to do multiqc on the bisulfite_sample results
.PHONY : bisulfite_sample_multiqc
bisulfite_sample_multiqc : $(DIR_MULTIQC)/bisulfite_sample/multiqc_report.html

$(DIR_MULTIQC)/bisulfite_sample/multiqc_report.html :	 $(patsubst %,$(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip,$(BISULFITE_SAMPLE_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_TRIM_FASTQCS)/%_trimmed_fastqc.zip,$(BISULFITE_SAMPLE_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.bedGraph.gz,$(BISULFITE_SAMPLE_PREFIXES))
	$(PATH_TO_MULTIQC) --force ./bisulfite --outdir $(@D)

########################################
# Rule to delete all temporary files from make bis_sample
.PHONY : clean_bisulfite_sample_tmp
clean_bisulfite_sample_tmp :
	rm -f $(BISFULITE_SAMPLE_CLEAN_TMP)

################################################################################
'

cat(make_var_bis_sample_prefix, file = file_make, sep = '\n', append = TRUE)
cat(make_var_bis_sample_clean_tmp, file = file_make, sep = '\n', append = TRUE)
cat(make_rule_bis_sample, file = file_make, sep = '\n', append = TRUE)

########################################################################
# OPTS for config.mk
config_bis_sample = '
################################################################################
# bisulfite_sample configuration options

# Command line option for minimum coverage required for bismark_methylation_extractor
# and scripts/classify_prepare_bisulfite_sample.awk in the sample classification module
OPT_MIN_COV = 5

# bismark_methylation_extractor
# For methylation extractor parameters see http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide_v0.15.0.pdf
OPTS_EXTRACTOR = --single-end --gzip --bedGraph --cutoff $(OPT_MIN_COV) --cytosine_report --genome_folder $(GENOME_PATH) --multicore 1
'
cat(config_bis_sample, file = file_config, sep='\n', append=T)

#######################################
# PBS script
# bisulfite_sample_q = c(
#	 '#!/bin/bash',
#	 '#### Begin PBS preamble',
#	 '#PBS -N bis_sample',
#	 '#PBS -l nodes=1:ppn=15,walltime=24:00:00,pmem=8gb',
#	 '#PBS -A sartor_lab',
#	 '#PBS -q first',
#	 '#PBS -M rcavalca@umich.edu',
#	 '#PBS -m abe',
#	 '#PBS -j oe',
#	 '#PBS -V',
#	 '#### End PBS preamble',
#	 '# Put your job commands after this line',
#	 sprintf('cd ~/latte/mint/projects/%s/',project),
#	 'make -j 3 bisulfite_sample')
# cat(bisulfite_sample_q, file=sprintf('projects/%s/pbs_jobs/bisulfite_sample.q', project), sep='\n')

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
		sprintf('track %s_bismark_simple_class', bisulfite_samples[i,'fullHumanID']),
		sprintf('parent %s_sample', bisulfite_samples[i,'humanID']),
		sprintf('bigDataUrl %s_bismark_simple_classification.bb', bisulfite_samples[i,'fullHumanID']),
		sprintf('shortLabel %s_bismark_simp_class', bisulfite_samples[i,'fullHumanID']),
		sprintf('longLabel %s_bismark_simple_classification', bisulfite_samples[i,'fullHumanID']),
		'visibility pack',
		'itemRgb on',
		'type bigBed 9 .',
		'priority 1.2',
		' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
}

}
