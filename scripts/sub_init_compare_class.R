################################################################################
# MAKEFILE: compare_classification rules

if(bool_bis_comp || bool_pull_comp) {

make_rule_compare_class_bis_module = '
# Intermediates for the bisulfite piece
.INTERMEDIATE : $(DIR_BIS_MSIG)/%_bisulfite_DMup.txt
$(DIR_BIS_MSIG)/%_bisulfite_DMup.txt : $(DIR_BIS_MSIG)/%_bisulfite_methylSig.txt
	awk -v OFS="\\t" \'NR > 1 && $$5 < 0.05 && $$7 > 0 { print $$1, $$2, $$3 }\' $< | sort -T . -k1,1 -k2,2n > $@

.INTERMEDIATE : $(DIR_BIS_MSIG)/%_bisulfite_DMdown.txt
$(DIR_BIS_MSIG)/%_bisulfite_DMdown.txt : $(DIR_BIS_MSIG)/%_bisulfite_methylSig.txt
	awk -v OFS="\\t" \'NR > 1 && $$5 < 0.05 && $$7 < 0 { print $$1, $$2, $$3 }\' $< | sort -T . -k1,1 -k2,2n > $@

.INTERMEDIATE : $(DIR_BIS_MSIG)/%_bisulfite_noDM_signal.txt
$(DIR_BIS_MSIG)/%_bisulfite_noDM_signal.txt : $(DIR_BIS_MSIG)/%_bisulfite_methylSig.txt
	awk -v OFS="\\t" \'NR > 1 && $$5 > 0.05 { print $$1, $$2, $$3 }\' $< | sort -T . -k1,1 -k2,2n > $@

.INTERMEDIATE : $(DIR_BIS_MSIG)/%_bisulfite_noDM_nosignal.txt
$(DIR_BIS_MSIG)/%_bisulfite_noDM_nosignal.txt : $(DIR_BIS_MSIG)/%_bisulfite_methylSig.txt
	bedtools complement -i <(awk -v OFS="\\t" \'NR > 1 { print $$1, $$2, $$3 }\' $<) -g <(sort -T . -k1,1 $(CHROM_PATH)) | sort -T . -k1,1 -k2,2n > $@
'

make_rule_compare_class_pull_module = '
# Intermediates for the pulldown piece
.INTERMEDIATE : $(DIR_PULL_PEPR)/%_pulldown_tmp_up.txt
$(DIR_PULL_PEPR)/%_pulldown_tmp_up.txt : $(DIR_PULL_PEPR)/%_pulldown__PePr_up_peaks.bed
	awk -v OFS="\\t" \'{print $$1, $$2, $$3}\' $< \\
	| sort -T . -k1,1 -k2,2n \\
	> $@

.INTERMEDIATE : $(DIR_PULL_PEPR)/%_pulldown_tmp_down.txt
$(DIR_PULL_PEPR)/%_pulldown_tmp_down.txt : $(DIR_PULL_PEPR)/%_pulldown__PePr_down_peaks.bed
	awk -v OFS="\\t" \'{print $$1, $$2, $$3}\' $< \\
	| sort -T . -k1,1 -k2,2n \\
	> $@

.INTERMEDIATE : $(DIR_PULL_PEPR)/%_pulldown_DMup.txt
$(DIR_PULL_PEPR)/%_pulldown_DMup.txt : $(DIR_PULL_PEPR)/%_pulldown_tmp_up.txt $(DIR_PULL_PEPR)/%_pulldown_tmp_down.txt
	bedops --difference $^ > $@

.INTERMEDIATE : $(DIR_PULL_PEPR)/%_pulldown_DMdown.txt
$(DIR_PULL_PEPR)/%_pulldown_DMdown.txt : $(DIR_PULL_PEPR)/%_pulldown_tmp_down.txt $(DIR_PULL_PEPR)/%_pulldown_tmp_up.txt
	bedops --difference $^ > $@

.INTERMEDIATE : $(DIR_PULL_PEPR)/%_pulldown_tmp_disjoint_DM.txt
$(DIR_PULL_PEPR)/%_pulldown_tmp_disjoint_DM.txt : $(DIR_PULL_PEPR)/%_pulldown_DMup.txt $(DIR_PULL_PEPR)/%_pulldown_DMdown.txt
	cat $^ | sort -T . -k1,1 -k2,2n > $@

.INTERMEDIATE : $(DIR_PULL_PEPR)/%_pulldown_tmp_disjoint_noDM.txt
$(DIR_PULL_PEPR)/%_pulldown_tmp_disjoint_noDM.txt: $(DIR_PULL_PEPR)/%_pulldown_tmp_disjoint_DM.txt
	bedtools complement -i $< -g <(sort -T . -k1,1 $(CHROM_PATH)) > $@

.INTERMEDIATE : $(DIR_PULL_PEPR)/%_pulldown_tmp_signal.txt
$(DIR_PULL_PEPR)/%_pulldown_tmp_signal.txt : $(DIR_PULL_PEPR)/%_pulldown_merged_signal.bed
	cp $< $@

.INTERMEDIATE : $(DIR_PULL_PEPR)/%_pulldown_tmp_nosignal.txt
$(DIR_PULL_PEPR)/%_pulldown_tmp_nosignal.txt : $(DIR_PULL_PEPR)/%_pulldown_tmp_signal.txt
	bedtools complement -i $< -g <(sort -T . -k1,1 $(CHROM_PATH)) > $@

.INTERMEDIATE : $(DIR_PULL_PEPR)/%_pulldown_DMdown.txt
$(DIR_PULL_PEPR)/%_pulldown_noDM_signal.txt : $(DIR_PULL_PEPR)/%_pulldown_tmp_disjoint_noDM.txt $(DIR_PULL_PEPR)/%_pulldown_tmp_signal.txt
	bedtools intersect -a $(word 1, $^) -b $(word 2, $^) | sort -T . -k1,1 -k2,2n > $@

.INTERMEDIATE : $(DIR_PULL_PEPR)/%_pulldown_noDM_nosignal.txt
$(DIR_PULL_PEPR)/%_pulldown_noDM_nosignal.txt : $(DIR_PULL_PEPR)/%_pulldown_tmp_disjoint_noDM.txt $(DIR_PULL_PEPR)/%_pulldown_tmp_nosignal.txt
	bedtools intersect -a $(word 1, $^) -b $(word 2, $^) | sort -T . -k1,1 -k2,2n > $@
'

# Collect

make_var_compare_class_prefix = sprintf('
################################################################################
# Workflow for compare_classification
COMPARE_CLASS_PREFIXES := %s',
	paste(unique(comparisons$humanID), collapse=' '))
cat(make_var_compare_class_prefix, file = file_make, sep = '\n', append = TRUE)

# The compare class type depends on the type of compares present
if(bool_bis_comp && bool_pull_comp) {
	compare_class_target = '$(DIR_CLASS_COMPARE)/%_compare_classification.bed : 	$(DIR_BIS_MSIG)/%_mc_hmc_bisulfite_DMup.txt \\
														$(DIR_BIS_MSIG)/%_mc_hmc_bisulfite_DMdown.txt \\
														$(DIR_BIS_MSIG)/%_mc_hmc_bisulfite_noDM_signal.txt \\
														$(DIR_BIS_MSIG)/%_mc_hmc_bisulfite_noDM_nosignal.txt \\
														$(DIR_PULL_PEPR)/%_hmc_pulldown_DMup.txt \\
														$(DIR_PULL_PEPR)/%_hmc_pulldown_DMdown.txt \\
														$(DIR_PULL_PEPR)/%_hmc_pulldown_noDM_signal.txt \\
														$(DIR_PULL_PEPR)/%_hmc_pulldown_noDM_nosignal.txt'
	rule1 = make_rule_compare_class_bis_module
	rule2 = make_rule_compare_class_pull_module
	class_script = '../../scripts/classify_compare.sh'
} else if (bool_bis_comp && !bool_pull_comp) {
	############################################################
	# NOTE: THIS IS NOT EXPLICITLY SUPPORTED RIGHT NOW
	############################################################
	compare_class_target = '$(DIR_CLASS_COMPARE)/%_compare_classification.bed : 	$(DIR_BIS_MSIG)/%_mc_bisulfite_DMup.txt \\
														$(DIR_BIS_MSIG)/%_mc_bisulfite_DMdown.txt \\
														$(DIR_BIS_MSIG)/%_mc_bisulfite_noDM_signal.txt \\
														$(DIR_BIS_MSIG)/%_mc_bisulfite_noDM_nosignal.txt \\
														$(DIR_BIS_MSIG)/%_hmc_bisulfite_DMup.txt \\
														$(DIR_BIS_MSIG)/%_hmc_bisulfite_DMdown.txt \\
														$(DIR_BIS_MSIG)/%_hmc_bisulfite_noDM_signal.txt \\
														$(DIR_BIS_MSIG)/%_hmc_bisulfite_noDM_nosignal.txt'
	rule1 = make_rule_compare_class_bis_module
	rule2 = ''
	class_script = '../../scripts/classify_compare.sh'
} else if (!bool_bis_comp && bool_pull_comp) {
	compare_class_target = '$(DIR_CLASS_COMPARE)/%_compare_classification.bed : 	$(DIR_PULL_PEPR)/%_mc_pulldown_DMup.txt \\
														$(DIR_PULL_PEPR)/%_mc_pulldown_DMdown.txt \\
														$(DIR_PULL_PEPR)/%_mc_pulldown_noDM_signal.txt \\
														$(DIR_PULL_PEPR)/%_mc_pulldown_noDM_nosignal.txt \\
														$(DIR_PULL_PEPR)/%_hmc_pulldown_DMup.txt \\
														$(DIR_PULL_PEPR)/%_hmc_pulldown_DMdown.txt \\
														$(DIR_PULL_PEPR)/%_hmc_pulldown_noDM_signal.txt \\
														$(DIR_PULL_PEPR)/%_hmc_pulldown_noDM_nosignal.txt'
	rule1 = make_rule_compare_class_pull_module
	rule2 = ''
	class_script = '../../scripts/classify_compare.sh'
}

make_rule_class_compare = sprintf('
# Master rule
.PHONY : compare_classification
compare_classification : 	$(patsubst %%,$(DIR_TRACK)/%%_compare_classification.bb,$(COMPARE_CLASS_PREFIXES)) \\
		$(patsubst %%,$(DIR_SUM_FIGURES)/%%_compare_class_counts.png,$(COMPARE_CLASS_PREFIXES)) \\
		$(patsubst %%,$(DIR_CLASS_COMPARE)/%%_compare_classification.bed,$(COMPARE_CLASS_PREFIXES))

# Rule for compare classification bigBed
$(DIR_TRACK)/%%_compare_classification.bb : $(DIR_CLASS_COMPARE)/%%_compare_classification.bed
	bedToBigBed $^ $(CHROM_PATH) $@

# Rule for annotatr of compare classification
$(DIR_SUM_FIGURES)/%%_compare_class_counts.png : $(DIR_CLASS_COMPARE)/%%_compare_class_for_annotatr.txt
	Rscript ../../scripts/annotatr_classification.R --file $< --genome $(GENOME)

.INTERMEDIATE : $(DIR_CLASS_COMPARE)/%%_compare_class_for_annotatr.txt
$(DIR_CLASS_COMPARE)/%%_compare_class_for_annotatr.txt : $(DIR_CLASS_COMPARE)/%%_compare_classification.bed
	cut -f 1-4 $< > $@

# NOTE: There is a known bug in make that incorrectly determines implicit intermediate
# files when they occur in a list of multiple targets and prerequisites.
# https://savannah.gnu.org/bugs/index.php?32042
# The easiest workaround is to make them precious and remove them

# Classification BED
.PRECIOUS : $(DIR_CLASS_COMPARE)/%%_compare_classification.bed
%s
	bash %s $(CHROM_PATH) $@ $^

%s
%s',
	compare_class_target, class_script, rule1, rule2)
cat(make_rule_class_compare, file = file_make, sep = '\n', append = TRUE)

#######################################
# PBS script
bisulfite_compare_q = c(
	'#!/bin/bash',
	'#### Begin PBS preamble',
	'#PBS -N class_compare',
	'#PBS -l procs=4,mem=32gb,walltime=6:00:00',
	'#PBS -A sartor_lab',
	'#PBS -q first',
	'#PBS -M rcavalca@umich.edu',
	'#PBS -m abe',
	'#PBS -j oe',
	'#PBS -V',
	'#### End PBS preamble',
	'# Put your job commands after this line',
	sprintf('cd ~/latte/mint/projects/%s/',project),
	'make -j 4 compare_classification')
cat(bisulfite_compare_q, file=sprintf('projects/%s/pbs_jobs/classify_compare.q', project), sep='\n')

for(comparison in unique(comparisons$humanID)) {
	# trackDb.txt entry for comparison classification
	trackEntry = c(
	  sprintf('track %s_compare_classification', comparison),
	  sprintf('parent %s_group_comparison', comparison),
	  sprintf('bigDataUrl %s_compare_classification.bb', comparison),
	  sprintf('shortLabel %s_comp_class', comparison),
	  sprintf('longLabel %s_compare_classification', comparison),
	  'visibility pack',
	  'itemRgb on',
	  'type bigBed 9 .',
	  'priority 1.1',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
}

}
