################################################################################
# MAKEFILE: sample_classification rules

if(bool_bis_samp || bool_pull_samp) {

make_rule_sample_class_bis_module = '
.INTERMEDIATE : $(DIR_BIS_BISMARK)/%_bisulfite_trimmed_bismark_bt2.CpG_report.txt
$(DIR_BIS_BISMARK)/%_bisulfite_trimmed_bismark_bt2.CpG_report.txt : $(DIR_BIS_BISMARK)/%_bisulfite_trimmed_bismark_bt2.CpG_report.txt.gz
	gunzip -c $< > $@

# Intermediates for the bisulfite piece
$(DIR_BIS_BISMARK)/%_bisulfite_highmeth.txt $(DIR_BIS_BISMARK)/%_bisulfite_lowmeth.txt $(DIR_BIS_BISMARK)/%_bisulfite_nometh_signal.txt $(DIR_BIS_BISMARK)/%_bisulfite_nometh_nosignal.txt : $(DIR_BIS_BISMARK)/%_bisulfite_trimmed_bismark_bt2.CpG_report.txt
	$(PATH_TO_AWK) -f ../../scripts/classify_prepare_bisulfite_sample.awk $<
'

make_rule_sample_class_pull_module = '
# Intermediates for the pulldown piece
.INTERMEDIATE : $(DIR_PULL_MACS)/%_pulldown_peak.txt
$(DIR_PULL_MACS)/%_pulldown_peak.txt : $(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak
	$(PATH_TO_AWK) -v OFS="\\t" \'{print $$1, $$2, $$3}\' $< \\
	| sort -T . -k1,1 -k2,2n \\
	> $@

.INTERMEDIATE : $(DIR_PULL_MACS)/%_pulldown_nopeak_signal.txt
$(DIR_PULL_MACS)/%_pulldown_nopeak_signal.txt : $(DIR_PULL_MACS)/%_pulldown_nopeak.txt $(DIR_PULL_MACS)/%_pulldown_signal.txt
	$(PATH_TO_BEDTOOLS) intersect -a $(word 1, $^) -b $(word 2, $^) \\
	| sort -T . -k1,1 -k2,2n \\
	> $@

.INTERMEDIATE : $(DIR_PULL_MACS)/%_pulldown_nopeak_nosignal.txt
$(DIR_PULL_MACS)/%_pulldown_nopeak_nosignal.txt : $(DIR_PULL_MACS)/%_pulldown_nopeak.txt $(DIR_PULL_MACS)/%_pulldown_nosignal.txt
	$(PATH_TO_BEDTOOLS) intersect -a $(word 1, $^) -b $(word 2, $^) \\
	| sort -T . -k1,1 -k2,2n \\
	> $@

.INTERMEDIATE : $(DIR_PULL_MACS)/%_pulldown_nopeak.txt
$(DIR_PULL_MACS)/%_pulldown_nopeak.txt : $(DIR_PULL_MACS)/%_pulldown_peak.txt
	$(PATH_TO_BEDTOOLS) complement -g <(sort -T . -k1,1 $(CHROM_PATH)) -i $< \\
	| sort -T . -k1,1 -k2,2n \\
	> $@

.INTERMEDIATE : $(DIR_PULL_MACS)/%_pulldown_signal.txt
$(DIR_PULL_MACS)/%_pulldown_signal.txt : $(DIR_PULL_COVERAGES)/%_input_pulldown_merged_coverage.bdg
	cp $< $@

.INTERMEDIATE : $(DIR_PULL_MACS)/%_pulldown_nosignal.txt
$(DIR_PULL_MACS)/%_pulldown_nosignal.txt : $(DIR_PULL_MACS)/%_pulldown_signal.txt
	$(PATH_TO_BEDTOOLS) complement -g <(sort -T . -k1,1 $(CHROM_PATH)) -i $< \\
	| sort -T . -k1,1 -k2,2n \\
	> $@
'

make_var_sample_class_prefix = sprintf('
################################################################################
# Workflow for sample_classification
SAMPLE_CLASS_PREFIXES := %s', paste(unique(samples$humanID), collapse=' '))
cat(make_var_sample_class_prefix, file = file_make, sep = '\n', append = TRUE)

# The sample class type depends on the type of samples present
if(bool_bis_samp && bool_pull_samp) {
	sample_class_tmps = '$(SAMPLE_CLASS_CLEAN_TMP) := $(patsubst %,$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_highmeth.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_lowmeth.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_nometh_signal.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_nometh_nosignal.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_MACS)/%_hmc_pulldown_peak.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_signal.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_nosignal.txt,$(SAMPLE_CLASS_PREFIXES))'
	sample_class_target = '$(DIR_CLASS_SAMPLE)/%_sample_classification.bed : 	$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_highmeth.txt \\
								$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_lowmeth.txt \\
								$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_nometh_signal.txt \\
								$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_nometh_nosignal.txt \\
								$(DIR_PULL_MACS)/%_hmc_pulldown_peak.txt \\
								$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_signal.txt \\
								$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_nosignal.txt'
	rule1 = make_rule_sample_class_bis_module
	rule2 = make_rule_sample_class_pull_module
	class_script = '../../scripts/classify_hybrid_sample.sh'
} else if (bool_bis_samp && !bool_pull_samp) {
	############################################################
	# NOTE: THIS IS NOT EXPLICITLY SUPPORTED RIGHT NOW
	############################################################
	sample_class_tmps = '$(SAMPLE_CLASS_CLEAN_TMP) := $(patsubst %,$(DIR_BIS_BISMARK)/%_mc_bisulfite_highmeth.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_BISMARK)/%_mc_bisulfite_lowmeth.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_BISMARK)/%_mc_bisulfite_nometh_signal.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_BISMARK)/%_mc_bisulfite_nometh_nosignal.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_BISMARK)/%_hmc_bisulfite_highmeth.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_BISMARK)/%_hmc_bisulfite_lowmeth.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_BISMARK)/%_hmc_bisulfite_nometh_signal.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_BISMARK)/%_hmc_bisulfite_nometh_nosignal.txt,$(SAMPLE_CLASS_PREFIXES))'
	sample_class_target = '$(DIR_CLASS_SAMPLE)/%_sample_classification.bed : 	$(DIR_BIS_BISMARK)/%_mc_bisulfite_highmeth.txt \\
								$(DIR_BIS_BISMARK)/%_mc_bisulfite_lowmeth.txt \\
								$(DIR_BIS_BISMARK)/%_mc_bisulfite_nometh_signal.txt \\
								$(DIR_BIS_BISMARK)/%_mc_bisulfite_nometh_nosignal.txt \\
								$(DIR_BIS_BISMARK)/%_hmc_bisulfite_highmeth.txt \\
								$(DIR_BIS_BISMARK)/%_hmc_bisulfite_lowmeth.txt \\
								$(DIR_BIS_BISMARK)/%_hmc_bisulfite_nometh_signal.txt \\
								$(DIR_BIS_BISMARK)/%_hmc_bisulfite_nometh_nosignal.txt'
	rule1 = make_rule_sample_class_bis_module
	rule2 = ''
	class_script = '../../scripts/classify_bisulfite_sample.sh'
} else {
	sample_class_tmps = '$(SAMPLE_CLASS_CLEAN_TMP) := $(patsubst %,$(DIR_PULL_MACS)/%_mc_pulldown_peak.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_MACS)/%_mc_pulldown_nopeak_signal.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_MACS)/%_mc_pulldown_nopeak_nosignal.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_MACS)/%_hmc_pulldown_peak.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_signal.txt,$(SAMPLE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_nosignal.txt,$(SAMPLE_CLASS_PREFIXES))'
	sample_class_target = '$(DIR_CLASS_SAMPLE)/%_sample_classification.bed : 	$(DIR_PULL_MACS)/%_mc_pulldown_peak.txt \\
								$(DIR_PULL_MACS)/%_mc_pulldown_nopeak_signal.txt \\
								$(DIR_PULL_MACS)/%_mc_pulldown_nopeak_nosignal.txt \\
								$(DIR_PULL_MACS)/%_hmc_pulldown_peak.txt \\
								$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_signal.txt \\
								$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_nosignal.txt'
	rule1 = make_rule_sample_class_pull_module
	rule2 = ''
	class_script = '../../scripts/classify_pulldown_sample.sh'
}

make_rule_class_sample = sprintf('
# Master rule
.PHONY : sample_classification
sample_classification : 	$(patsubst %%,$(DIR_TRACK)/%%_sample_classification.bb,$(SAMPLE_CLASS_PREFIXES)) \\
		$(patsubst %%,$(DIR_SUM_FIGURES)/%%_sample_class_counts.png,$(SAMPLE_CLASS_PREFIXES)) \\
		$(patsubst %%,$(DIR_CLASS_SAMPLE)/%%_sample_classification.bed,$(SAMPLE_CLASS_PREFIXES))

# Rule for sample classification bigBed
$(DIR_TRACK)/%%_sample_classification.bb : $(DIR_CLASS_SAMPLE)/%%_sample_classification.bed
	$(PATH_TO_BDG2BB) $^ $(CHROM_PATH) $@

# Rule for annotatr of sample classification
$(DIR_SUM_FIGURES)/%%_sample_class_counts.png : $(DIR_CLASS_SAMPLE)/%%_sample_class_for_annotatr.txt
	$(PATH_TO_R) ../../scripts/annotatr_classification.R --file $< --genome $(GENOME)

.INTERMEDIATE : $(DIR_CLASS_SAMPLE)/%%_sample_class_for_annotatr.txt
$(DIR_CLASS_SAMPLE)/%%_sample_class_for_annotatr.txt : $(DIR_CLASS_SAMPLE)/%%_sample_classification.bed
	cut -f 1-4 $< > $@

# Classification BED
%s
	bash %s $(PATH_TO_BEDTOOLS) $(PATH_TO_AWK) $(CHROM_PATH) $@ $^

%s
%s

# Clean temporary files that make does not clean up
%s

.PHONY : sample_classification_clean_tmp
sample_classification_clean_tmp :
	rm -f $(SAMPLE_CLASS_CLEAN_TMP)',
	sample_class_target, class_script, rule1, rule2, sample_class_tmps)
cat(make_rule_class_sample, file = file_make, sep = '\n', append = TRUE)

#######################################
# PBS script
# pulldown_sample_q = c(
# 	'#!/bin/bash',
# 	'#### Begin PBS preamble',
# 	'#PBS -N class_sample',
# 	'#PBS -l procs=4,mem=48gb,walltime=6:00:00',
# 	'#PBS -A sartor_lab',
# 	'#PBS -q first',
# 	'#PBS -M rcavalca@umich.edu',
# 	'#PBS -m abe',
# 	'#PBS -j oe',
# 	'#PBS -V',
# 	'#### End PBS preamble',
# 	'# Put your job commands after this line',
# 	sprintf('cd ~/latte/mint/projects/%s/',project),
# 	'make -j 4 sample_classification')
# cat(pulldown_sample_q, file=sprintf('projects/%s/pbs_jobs/classify_sample.q', project), sep='\n')

for(sample in unique(samples$humanID)) {
	# trackDb.txt entry for sample classification
	trackEntry = c(
	  sprintf('track %s_sample_classification', sample),
	  sprintf('parent %s_sample', sample),
	  sprintf('bigDataUrl %s_sample_classification.bb', sample),
	  sprintf('shortLabel %s_sample_class', sample),
	  sprintf('longLabel %s_sample_classification', sample),
	  'visibility pack',
	  'itemRgb on',
	  'type bigBed 9 .',
	  'priority 1.1',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
}

}
