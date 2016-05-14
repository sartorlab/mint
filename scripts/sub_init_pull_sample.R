################################################################################
# MAKEFILE: pulldown_sample rules

if(bool_pull_samp) {

make_var_pull_samp_prefix = sprintf('
################################################################################
# Workflow for pulldown_sample
PULLDOWN_SAMPLE_PREFIXES := %s', paste(pulldown_samples_noinput$fullHumanID, collapse=' '))


make_var_pull_samp = 'PULLDOWN_SAMPLE_PREREQS :=	$(patsubst %,$(DIR_TRACK)/%_simple_classification.bb,$(PULLDOWN_SAMPLE_PREFIXES)) \\
						$(patsubst %,$(DIR_SUM_FIGURES)/%_simple_class_counts.png,$(PULLDOWN_SAMPLE_PREFIXES)) \\
						$(patsubst %,$(DIR_CLASS_SIMPLE)/%_simple_classification.bed,$(PULLDOWN_SAMPLE_PREFIXES)) \\
						$(patsubst %,$(DIR_TRACK)/%_macs2_peaks.bb,$(PULLDOWN_SAMPLE_PREFIXES)) \\
						$(patsubst %,$(DIR_PULL_MACS)/%_macs2_peaks.narrowPeak,$(PULLDOWN_SAMPLE_PREFIXES))'

make_var_pull_samp_clean_tmp = 'PULLDOWN_SAMPLE_CLEAN_TMP := $(patsubst %,$(DIR_CLASS_SIMPLE)/%_pulldown_simple_class_for_annotatr.txt,$(PULLDOWN_SAMPLE_PREFIXES)) \\
						$(patsubst %,$(DIR_PULL_MACS)/%_macs2_peaks_tmp.narrowPeak,$(PULLDOWN_SAMPLE_PREFIXES))
'

# NOTE: This cannot be indented because they would mess up the makefile
make_rule_pull_samp = '
.PHONY : pulldown_sample
pulldown_sample : pulldown_align $(PULLDOWN_SAMPLE_PREREQS)

# Rule for UCSC bigBed track of simple classifiation
$(DIR_TRACK)/%_pulldown_simple_classification.bb : $(DIR_CLASS_SIMPLE)/%_pulldown_simple_classification.bed
	$(PATH_TO_BDG2BB) $< $(CHROM_PATH) $@

# Rule for annotatr of simple classification
$(DIR_SUM_FIGURES)/%_pulldown_simple_class_counts.png : $(DIR_CLASS_SIMPLE)/%_pulldown_simple_class_for_annotatr.txt
	$(PATH_TO_R) ../../scripts/annotatr_classification.R --file $< --genome $(GENOME)

.INTERMEDIATE : $(DIR_CLASS_SIMPLE)/%_pulldown_simple_class_for_annotatr.txt
$(DIR_CLASS_SIMPLE)/%_pulldown_simple_class_for_annotatr.txt : $(DIR_CLASS_SIMPLE)/%_pulldown_simple_classification.bed
	$(PATH_TO_AWK) -v OFS="\\t" \'{ print $$1, $$2, $$3, $$4 }\' $< > $@

# Simple classification
$(DIR_CLASS_SIMPLE)/%_pulldown_simple_classification.bed : $(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak
	$(PATH_TO_R) ../../scripts/classify_simple.R --project $(PROJECT) --inFile $< --outFile $@

# Rule for UCSC bigBed track
$(DIR_TRACK)/%_macs2_peaks.bb : $(DIR_PULL_MACS)/%_macs2_peaks_tmp.narrowPeak
	$(PATH_TO_BDG2BB) -type=bed6+4 -as=narrowPeak.as $^ $(CHROM_PATH) $@

# Rule for macs2 peak fix
.INTERMEDIATE : $(DIR_PULL_MACS)/%_macs2_peaks_tmp.narrowPeak
$(DIR_PULL_MACS)/%_macs2_peaks_tmp.narrowPeak : $(DIR_PULL_MACS)/%_macs2_peaks.narrowPeak
	$(PATH_TO_AWK) -f ../../scripts/macs_fix_narrowPeak.awk $^ | sort -T . -k1,1 -k2,2n > $@

# Rule for macs2 peaks
$(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak : 	$(DIR_PULL_BOWTIE2)/%_pulldown_trimmed.fq.gz_aligned.bam \\
														$(DIR_PULL_BOWTIE2)/%_input_pulldown_trimmed.fq.gz_aligned.bam
	$(PATH_TO_MACS) callpeak -t $(word 1, $^) -c $(word 2, $^) -f BAM -g hs --name $(patsubst %_peaks.narrowPeak,%,$(@F)) --outdir $(@D)

# Rule to delete all temporary files from make bis_align
.PHONY : clean_pulldown_sample_tmp
clean_pulldown_sample_tmp :
	rm -f $(PULLDOWN_SAMPLE_CLEAN_TMP)
'

cat(make_var_pull_samp_prefix, file = file_make, sep = '\n', append = TRUE)
cat(make_var_pull_samp, file = file_make, sep = '\n', append = TRUE)
cat(make_var_pull_samp_clean_tmp, file = file_make, sep = '\n', append = TRUE)
cat(make_rule_pull_samp, file = file_make, sep = '\n', append = TRUE)

#######################################
# PBS script
# pulldown_sample_q = c(
# 	'#!/bin/bash',
# 	'#### Begin PBS preamble',
# 	'#PBS -N pull_sample',
# 	'#PBS -l procs=4,mem=32gb,walltime=6:00:00',
# 	'#PBS -A sartor_lab',
# 	'#PBS -q first',
# 	'#PBS -M rcavalca@umich.edu',
# 	'#PBS -m abe',
# 	'#PBS -j oe',
# 	'#PBS -V',
# 	'#### End PBS preamble',
# 	'# Put your job commands after this line',
# 	sprintf('cd ~/latte/mint/projects/%s/',project),
# 	'make -j 4 pulldown_sample')
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
