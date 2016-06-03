################################################################################
# MAKEFILE: pulldown_sample rules

if(bool_pull_samp) {

make_var_pull_samp_prefix = sprintf('
################################################################################
# Workflow for pulldown_sample

PULLDOWN_SAMPLE_PREFIXES := %s', paste(pulldown_samples_noinput$fullHumanID, collapse=' '))

make_var_pull_samp_clean_tmp = 'PULLDOWN_SAMPLE_CLEAN_TMP := $(patsubst %,$(DIR_CLASS_SIMPLE)/%_pulldown_simple_class_for_annotatr.txt,$(PULLDOWN_SAMPLE_PREFIXES)) \\
						$(patsubst %,$(DIR_PULL_MACS)/%_macs2_peaks_tmp.narrowPeak,$(PULLDOWN_SAMPLE_PREFIXES)) \\
						$(patsubst %,$(DIR_PULL_MACS)/%_pulldown_macs2_peaks_for_annotatr.txt,$(PULLDOWN_SAMPLE_PREFIXES)) \\
'

# NOTE: This cannot be indented because they would mess up the makefile
make_rule_pull_samp = '########################################

.PHONY : pulldown_sample
pulldown_sample : pulldown_align pulldown_macs2 pulldown_simple_classification

########################################
.PHONY : pulldown_macs2
pulldown_macs2 : 	$(patsubst %,$(DIR_PULL_MACS)/%_macs2_peaks.narrowPeak,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_PULL_MACS)/%_macs2_model.r,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_PULL_MACS)/%_macs2_model.pdf,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_RDATA)/%_macs2_peaks_annotatr_analysis.RData,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_TRACK)/%_macs2_peaks.bb,$(PULLDOWN_SAMPLE_PREFIXES))

# Rule for macs2 peaks
$(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak $(DIR_PULL_MACS)/%_pulldown_macs2_model.r : 	$(DIR_PULL_BOWTIE2)/%_pulldown_trimmed.fq.gz_aligned.bam \\
														$(DIR_PULL_BOWTIE2)/%_input_pulldown_trimmed.fq.gz_aligned.bam
	$(PATH_TO_MACS) callpeak -t $(word 1, $^) -c $(word 2, $^) -f BAM --name $(patsubst %_peaks.narrowPeak,%,$(@F)) --outdir $(@D) $(OPTS_MACS)

# Rule for macs2 model
$(DIR_PULL_MACS)/%_pulldown_macs2_model.pdf : $(DIR_PULL_MACS)/%_pulldown_macs2_model.r
	cd $(DIR_PULL_MACS); \\
	Rscript $(<F)

# Rule for annotatr input of macs2 peaks
# NOTE: Using fold change ($7) and p-value ($8)
.INTERMEDIATE : $(DIR_PULL_MACS)/%_pulldown_macs2_peaks_for_annotatr.txt
$(DIR_PULL_MACS)/%_pulldown_macs2_peaks_for_annotatr.txt : $(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak
	$(PATH_TO_AWK) -v OFS="\\t" \'{ print $$1, $$2, $$3, $$4, $$7, "*", $$8 }\' $< > $@

# Rule for annotatr of macs2 narrowPeak
$(DIR_RDATA)/%_pulldown_macs2_peaks_annotatr_analysis.RData : $(DIR_PULL_MACS)/%_pulldown_macs2_peaks_for_annotatr.txt
	$(PATH_TO_R) ../../scripts/annotatr_classification.R --file $< --genome $(GENOME)

# Rule to cap macs2 narrowPeaks at 1000 for bigBed
.INTERMEDIATE : $(DIR_PULL_MACS)/%_macs2_peaks_tmp.narrowPeak
$(DIR_PULL_MACS)/%_macs2_peaks_tmp.narrowPeak : $(DIR_PULL_MACS)/%_macs2_peaks.narrowPeak
	$(PATH_TO_AWK) -f ../../scripts/macs_fix_narrowPeak.awk $^ | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@

# Rule for UCSC bigBed track
$(DIR_TRACK)/%_macs2_peaks.bb : $(DIR_PULL_MACS)/%_macs2_peaks_tmp.narrowPeak
	$(PATH_TO_BDG2BB) -type=bed6+4 -as=narrowPeak.as $^ $(CHROM_PATH) $@

########################################
.PHONY : pulldown_simple_classification
pulldown_simple_classification : 	$(patsubst %,$(DIR_CLASS_SIMPLE)/%_simple_classification.bed,$(PULLDOWN_SAMPLE_PREFIXES)) \\
									$(patsubst %,$(DIR_RDATA)/%_simple_class_annotatr_analysis.RData,$(PULLDOWN_SAMPLE_PREFIXES)) \\
									$(patsubst %,$(DIR_TRACK)/%_simple_classification.bb,$(PULLDOWN_SAMPLE_PREFIXES))

# Rule for simple classification of macs2 peaks
$(DIR_CLASS_SIMPLE)/%_pulldown_simple_classification.bed : $(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak
	$(PATH_TO_R) ../../scripts/classify_simple.R --project $(PROJECT) --inFile $< --outFile $@

# Rule for annotatr input of simple classification
.INTERMEDIATE : $(DIR_CLASS_SIMPLE)/%_pulldown_simple_class_for_annotatr.txt
$(DIR_CLASS_SIMPLE)/%_pulldown_simple_class_for_annotatr.txt : $(DIR_CLASS_SIMPLE)/%_pulldown_simple_classification.bed
	$(PATH_TO_AWK) -v OFS="\\t" \'{ print $$1, $$2, $$3, $$4 }\' $< > $@

# Rule for annotatr of simple classification
$(DIR_RDATA)/%_pulldown_simple_class_annotatr_analysis.RData : $(DIR_CLASS_SIMPLE)/%_pulldown_simple_class_for_annotatr.txt
	$(PATH_TO_R) ../../scripts/annotatr_classification.R --file $< --genome $(GENOME)

# Rule for UCSC bigBed track of simple classifiation
$(DIR_TRACK)/%_pulldown_simple_classification.bb : $(DIR_CLASS_SIMPLE)/%_pulldown_simple_classification.bed
	$(PATH_TO_BDG2BB) $< $(CHROM_PATH) $@

########################################
# Rule to delete all temporary files from make pulldown_sample
.PHONY : clean_pulldown_sample_tmp
clean_pulldown_sample_tmp :
	rm -f $(PULLDOWN_SAMPLE_CLEAN_TMP)

################################################################################'

cat(make_var_pull_samp_prefix, file = file_make, sep = '\n', append = TRUE)
cat(make_var_pull_samp_clean_tmp, file = file_make, sep = '\n', append = TRUE)
cat(make_rule_pull_samp, file = file_make, sep = '\n', append = TRUE)

########################################################################
# OPTS for config.mk
if(genome == 'hg19' || genome == 'hg38') {
	macs_genome = 'hs'
} else if (genome == 'mm9' || genome == 'mm10') {
	macs_genome = 'mm'
} else {
	macs_genome = 'SPECIFY'
}

config_pull_sample = sprintf('################################################################################
# pulldown_sample configuration options

# macs2
# NOTE: Please ensure the genome size matches the organism in the study
# hg19, hg38, mm9, and mm10 are automatically populated.
# For documentation about parameters see https://github.com/taoliu/MACS
OPTS_MACS = --gsize %s --qvalue 0.01 --mfold 5 50
',
	macs_genome)
cat(config_pull_sample, file = file_config, sep='\n', append=T)

#######################################
# PBS script
# pulldown_sample_q = c(
# 	'#!/bin/bash',
# 	'#### Begin PBS preamble',
# 	'#PBS -N pull_sample',
# 	'#PBS -l nodes=1:ppn=6,walltime=24:00:00,pmem=8gb',
# 	'#PBS -A sartor_lab',
# 	'#PBS -q first',
# 	'#PBS -M rcavalca@umich.edu',
# 	'#PBS -m abe',
# 	'#PBS -j oe',
# 	'#PBS -V',
# 	'#### End PBS preamble',
# 	'# Put your job commands after this line',
# 	sprintf('cd ~/latte/mint/projects/%s/',project),
# 	'make -j 6 pulldown_sample')
# cat(pulldown_sample_q, file=sprintf('projects/%s/pbs_jobs/pulldown_sample.q', project), sep='\n')

for(i in 1:nrow(pulldown_samples_noinput)) {
	# trackDb.txt entry for MACS2 output
	trackEntry = c(
	  sprintf('track %s_peaks', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('parent %s_sample', pulldown_samples_noinput[i,'humanID']),
	  sprintf('bigDataUrl %s_macs2_peaks.bb', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('shortLabel %s_peaks', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('longLabel %s_macs2_peaks', pulldown_samples_noinput[i,'fullHumanID']),
	  'visibility dense',
	  'type bigBed',
	  'priority 1.5',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

	# trackDb.txt entry for pulldown simple classification output
	trackEntry = c(
	  sprintf('track %s_simple_class', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('parent %s_sample', pulldown_samples_noinput[i,'humanID']),
	  sprintf('bigDataUrl %s_simple_classification.bb', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('shortLabel %s_simp_class', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('longLabel %s_simple_classification', pulldown_samples_noinput[i,'fullHumanID']),
	  'visibility pack',
	  'itemRgb on',
	  'type bigBed 9 .',
	  'priority 1.2',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
}

}
