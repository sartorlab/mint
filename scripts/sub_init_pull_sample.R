################################################################################
# MAKEFILE: pulldown_sample rules

if(bool_pull_samp) {

make_var_pull_samp_prefix = sprintf('
################################################################################
# Workflow for pulldown_sample

PULLDOWN_SAMPLE_PREFIXES := %s', paste(pulldown_samples_noinput$fullHumanID, collapse=' '))

# NOTE: This cannot be indented because they would mess up the makefile
make_rule_pull_samp = '########################################

.PHONY : pulldown_sample
pulldown_sample :	 $(patsubst %,$(DIR_PULL_MACS)/%_macs2_peaks.narrowPeak,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_PULL_MACS)/%_macs2_model.r,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_PULL_MACS)/%_macs2_model.pdf,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_RDATA)/%_macs2_annotatr_analysis.RData,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_CLASS_SIMPLE)/%_macs2_simple_classification.bed,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_RDATA)/%_macs2_simple_classification_annotatr_analysis.RData,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_TRACK)/%_macs2_simple_classification.bb,$(PULLDOWN_SAMPLE_PREFIXES))

########################################
.PHONY : pulldown_macs2
pulldown_macs2 :	 $(patsubst %,$(DIR_PULL_MACS)/%_macs2_peaks.narrowPeak,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_PULL_MACS)/%_macs2_model.r,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_PULL_MACS)/%_macs2_model.pdf,$(PULLDOWN_SAMPLE_PREFIXES)) \\
					$(patsubst %,$(DIR_RDATA)/%_macs2_annotatr_analysis.RData,$(PULLDOWN_SAMPLE_PREFIXES)) \\

# Rule for macs2 peaks
$(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak $(DIR_PULL_MACS)/%_pulldown_macs2_model.r :	 $(DIR_PULL_BOWTIE2)/%_pulldown_trimmed.fq.gz_aligned.bam \\
														$(DIR_PULL_BOWTIE2)/%_input_pulldown_trimmed.fq.gz_aligned.bam
	$(PATH_TO_MACS) callpeak -t $(word 1, $^) -c $(word 2, $^) -f BAM --name $(patsubst %_peaks.narrowPeak,%,$(@F)) --outdir $(@D) $(OPTS_MACS)

# Rule for macs2 model
$(DIR_PULL_MACS)/%_pulldown_macs2_model.pdf : $(DIR_PULL_MACS)/%_pulldown_macs2_model.r
	cd $(DIR_PULL_MACS); \\
	Rscript $(<F)

# Rule for annotatr of macs2 narrowPeak
$(DIR_RDATA)/%_pulldown_macs2_annotatr_analysis.RData : $(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak
	$(PATH_TO_R) ../../scripts/annotatr_annotations.R --file $< --genome $(GENOME) --annot_type macs2 --group1 NULL --group0 NULL

########################################
.PHONY : pulldown_simple_classification
pulldown_simple_classification :	 $(patsubst %,$(DIR_CLASS_SIMPLE)/%_macs2_simple_classification.bed,$(PULLDOWN_SAMPLE_PREFIXES)) \\
									$(patsubst %,$(DIR_RDATA)/%_macs2_simple_classification_annotatr_analysis.RData,$(PULLDOWN_SAMPLE_PREFIXES)) \\
									$(patsubst %,$(DIR_TRACK)/%_macs2_simple_classification.bb,$(PULLDOWN_SAMPLE_PREFIXES))

# Rule for simple classification of macs2 peaks
$(DIR_CLASS_SIMPLE)/%_pulldown_macs2_simple_classification.bed : $(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak
	$(PATH_TO_R) ../../scripts/classify_simple.R --project $(PROJECT) --inFile $< --outFile $@ --group1 NULL --group0 NULL

# Rule for annotatr of simple classification
$(DIR_RDATA)/%_pulldown_macs2_simple_classification_annotatr_analysis.RData : $(DIR_CLASS_SIMPLE)/%_pulldown_macs2_simple_classification.bed
	$(PATH_TO_R) ../../scripts/annotatr_annotations.R --file $< --genome $(GENOME) --annot_type simple_pulldown_macs2 --group1 NULL --group0 NULL

# Rule for UCSC bigBed track of simple classifiation
$(DIR_TRACK)/%_pulldown_macs2_simple_classification.bb : $(DIR_CLASS_SIMPLE)/%_pulldown_macs2_simple_classification.bed
	$(PATH_TO_BED2BB) $< $(CHROM_PATH) $@

################################################################################'

cat(make_var_pull_samp_prefix, file = file_make, sep = '\n', append = TRUE)
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
#	 '#!/bin/bash',
#	 '#### Begin PBS preamble',
#	 '#PBS -N pull_sample',
#	 '#PBS -l nodes=1:ppn=6,walltime=24:00:00,pmem=8gb',
#	 '#PBS -A sartor_lab',
#	 '#PBS -q first',
#	 '#PBS -M rcavalca@umich.edu',
#	 '#PBS -m abe',
#	 '#PBS -j oe',
#	 '#PBS -V',
#	 '#### End PBS preamble',
#	 '# Put your job commands after this line',
#	 sprintf('cd ~/latte/mint/projects/%s/',project),
#	 'make -j 6 pulldown_sample')
# cat(pulldown_sample_q, file=sprintf('projects/%s/pbs_jobs/pulldown_sample.q', project), sep='\n')

for(i in 1:nrow(pulldown_samples_noinput)) {

	# trackDb.txt entry for pulldown simple classification output
	trackEntry = c(
		sprintf('track %s_simple_class', pulldown_samples_noinput[i,'fullHumanID']),
		sprintf('parent %s_sample', pulldown_samples_noinput[i,'humanID']),
		sprintf('bigDataUrl %s_macs2_simple_classification.bb', pulldown_samples_noinput[i,'fullHumanID']),
		sprintf('shortLabel %s_macs2_simp_class', pulldown_samples_noinput[i,'fullHumanID']),
		sprintf('longLabel %s_macs2_simple_classification', pulldown_samples_noinput[i,'fullHumanID']),
		'visibility pack',
		'itemRgb on',
		'type bigBed 9 .',
		'priority 1.2',
		' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
}

}
