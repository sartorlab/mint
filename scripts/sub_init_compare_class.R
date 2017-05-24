################################################################################
# MAKEFILE: compare_classification rules

if(bool_bis_comp || bool_pull_comp) {

make_rule_compare_class_bis_module = '
# Intermediates for the bisulfite piece
# Each needs 0-based start and 1-based end to match other files
.INTERMEDIATE : $(DIR_BIS_DSS)/%_bisulfite_DMup.txt
$(DIR_BIS_DSS)/%_bisulfite_DMup.txt : $(DIR_BIS_DSS)/%_bisulfite_dss_significant.txt
	$(PATH_TO_AWK) -v OFS="\\t" \'NR > 1 && $$7 >= 0 { print $$1, $$2 - 1, $$3 }\' $< | sort -T $(DIR_TMP) -k1,1 -k2,2n | bedtools merge > $@

.INTERMEDIATE : $(DIR_BIS_DSS)/%_bisulfite_DMdown.txt
$(DIR_BIS_DSS)/%_bisulfite_DMdown.txt : $(DIR_BIS_DSS)/%_bisulfite_dss_significant.txt
	$(PATH_TO_AWK) -v OFS="\\t" \'NR > 1 && $$7 < 0 { print $$1, $$2 - 1, $$3 }\' $< | sort -T $(DIR_TMP) -k1,1 -k2,2n | bedtools merge > $@

.INTERMEDIATE : $(DIR_BIS_DSS)/%_bisulfite_noDM_signal.txt
$(DIR_BIS_DSS)/%_bisulfite_noDM_signal.txt : $(DIR_BIS_DSS)/%_bisulfite_dss_all.txt $(DIR_BIS_DSS)/%_bisulfite_dss_significant.txt
	$(PATH_TO_BEDTOOLS) intersect -a <($(PATH_TO_AWK) -v OFS="\\t" \'NR > 1 { print $$1, $$2 - 1, $$3 }\' $(word 1, $^)) -b <($(PATH_TO_AWK) -v OFS="\\t" \'NR > 1 { print $$1, $$2 - 1, $$3 }\' $(word 2, $^)) -v | sort -T $(DIR_TMP) -k1,1 -k2,2n | bedtools merge > $@

.INTERMEDIATE : $(DIR_BIS_DSS)/%_bisulfite_noDM_nosignal.txt
$(DIR_BIS_DSS)/%_bisulfite_noDM_nosignal.txt : $(DIR_BIS_DSS)/%_bisulfite_dss_all.txt
	$(PATH_TO_BEDTOOLS) complement -i <($(PATH_TO_AWK) -v OFS="\\t" \'NR > 1 { print $$1, $$2 - 1, $$3 }\' <(sort -T $(DIR_TMP) -k1,1 -k2,2n $<)) -g <(sort -T $(DIR_TMP) -k1,1 $(CHROM_PATH)) | sort -T $(DIR_TMP) -k1,1 -k2,2n | bedtools merge > $@
'

make_rule_compare_class_pull_module = '
# Intermediates for the pulldown piece
.INTERMEDIATE : $(DIR_PULL_CSAW)/%_pulldown_tmp_up.txt
$(DIR_PULL_CSAW)/%_pulldown_tmp_up.txt : $(DIR_PULL_CSAW)/%_pulldown_csaw_significant.txt
	$(PATH_TO_AWK) -v OFS="\\t" \'NR > 1 && $$8 >= 0 {print $$1, $$2, $$3}\' $< \\
	| sort -T $(DIR_TMP) -k1,1 -k2,2n \\
	> $@

.INTERMEDIATE : $(DIR_PULL_CSAW)/%_pulldown_tmp_down.txt
$(DIR_PULL_CSAW)/%_pulldown_tmp_down.txt : $(DIR_PULL_CSAW)/%_pulldown_csaw_significant.txt
	$(PATH_TO_AWK) -v OFS="\\t" \'NR > 1 && $$8 < 0 {print $$1, $$2, $$3}\' $< \\
	| sort -T $(DIR_TMP) -k1,1 -k2,2n \\
	> $@


# Each needs 0-based start and 1-based end to match other files
.INTERMEDIATE : $(DIR_PULL_CSAW)/%_pulldown_DMup.txt
$(DIR_PULL_CSAW)/%_pulldown_DMup.txt : $(DIR_PULL_CSAW)/%_pulldown_tmp_up.txt $(DIR_PULL_CSAW)/%_pulldown_tmp_down.txt
	$(PATH_TO_BEDOPS) --difference $^ | bedtools merge > $@

.INTERMEDIATE : $(DIR_PULL_CSAW)/%_pulldown_DMdown.txt
$(DIR_PULL_CSAW)/%_pulldown_DMdown.txt : $(DIR_PULL_CSAW)/%_pulldown_tmp_down.txt $(DIR_PULL_CSAW)/%_pulldown_tmp_up.txt
	$(PATH_TO_BEDOPS) --difference $^ | bedtools merge > $@

.INTERMEDIATE : $(DIR_PULL_CSAW)/%_pulldown_tmp_disjoint_DM.txt
$(DIR_PULL_CSAW)/%_pulldown_tmp_disjoint_DM.txt : $(DIR_PULL_CSAW)/%_pulldown_DMup.txt $(DIR_PULL_CSAW)/%_pulldown_DMdown.txt
	cat $^ | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@

.INTERMEDIATE : $(DIR_PULL_CSAW)/%_pulldown_tmp_disjoint_noDM.txt
$(DIR_PULL_CSAW)/%_pulldown_tmp_disjoint_noDM.txt: $(DIR_PULL_CSAW)/%_pulldown_tmp_disjoint_DM.txt
	$(PATH_TO_BEDTOOLS) complement -i $< -g <(sort -T $(DIR_TMP) -k1,1 $(CHROM_PATH)) \\
	| $(PATH_TO_AWK) -v OFS="\\t" \'$$2 != $$3 {print $$0}\' \\
	> $@

.INTERMEDIATE : $(DIR_PULL_CSAW)/%_pulldown_tmp_signal.txt
$(DIR_PULL_CSAW)/%_pulldown_tmp_signal.txt : $(DIR_PULL_CSAW)/%_pulldown_merged_signal.bed
	cp $< $@

.INTERMEDIATE : $(DIR_PULL_CSAW)/%_pulldown_tmp_nosignal.txt
$(DIR_PULL_CSAW)/%_pulldown_tmp_nosignal.txt : $(DIR_PULL_CSAW)/%_pulldown_tmp_signal.txt
	$(PATH_TO_BEDTOOLS) complement -i $< -g <(sort -T $(DIR_TMP) -k1,1 $(CHROM_PATH)) \\
	| $(PATH_TO_AWK) -v OFS="\\t" \'$$2 != $$3 {print $$0}\' \\
	> $@

.INTERMEDIATE : $(DIR_PULL_CSAW)/%_pulldown_noDM_signal.txt
$(DIR_PULL_CSAW)/%_pulldown_noDM_signal.txt : $(DIR_PULL_CSAW)/%_pulldown_tmp_disjoint_noDM.txt $(DIR_PULL_CSAW)/%_pulldown_tmp_signal.txt
	$(PATH_TO_BEDTOOLS) intersect -a $(word 1, $^) -b $(word 2, $^) | sort -T $(DIR_TMP) -k1,1 -k2,2n | bedtools merge > $@

.INTERMEDIATE : $(DIR_PULL_CSAW)/%_pulldown_noDM_nosignal.txt
$(DIR_PULL_CSAW)/%_pulldown_noDM_nosignal.txt : $(DIR_PULL_CSAW)/%_pulldown_tmp_disjoint_noDM.txt $(DIR_PULL_CSAW)/%_pulldown_tmp_nosignal.txt
	$(PATH_TO_BEDTOOLS) intersect -a $(word 1, $^) -b $(word 2, $^) | sort -T $(DIR_TMP) -k1,1 -k2,2n | bedtools merge > $@
'

# Collect

make_var_compare_class_prefix = sprintf('
################################################################################
# Workflow for compare_classification
COMPARE_CLASS_PREFIXES := %s',
	paste(unique(comparisons$comparison), collapse=' '))
cat(make_var_compare_class_prefix, file = file_make, sep = '\n', append = TRUE)

# The compare class type depends on the type of compares present
if(bool_bis_comp && bool_pull_comp) {
	compare_class_tmps = 'COMPARE_CLASS_CLEAN_TMP := $(patsubst %,$(DIR_BIS_DSS)/%_mc_hmc_bisulfite_DMup.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_DSS)/%_mc_hmc_bisulfite_DMdown.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_DSS)/%_mc_hmc_bisulfite_noDM_signal.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_DSS)/%_mc_hmc_bisulfite_noDM_nosignal.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_CSAW)/%_hmc_pulldown_DMup.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_CSAW)/%_hmc_pulldown_DMdown.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_CSAW)/%_hmc_pulldown_noDM_signal.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_CSAW)/%_hmc_pulldown_noDM_nosignal.txt,$(COMPARE_CLASS_PREFIXES)) \\'
	compare_class_target = '$(DIR_CLASS_COMPARE)/%_compare_classification.bed :	 $(DIR_BIS_DSS)/%_mc_hmc_bisulfite_DMup.txt \\
								$(DIR_BIS_DSS)/%_mc_hmc_bisulfite_DMdown.txt \\
								$(DIR_BIS_DSS)/%_mc_hmc_bisulfite_noDM_signal.txt \\
								$(DIR_BIS_DSS)/%_mc_hmc_bisulfite_noDM_nosignal.txt \\
								$(DIR_PULL_CSAW)/%_hmc_pulldown_DMup.txt \\
								$(DIR_PULL_CSAW)/%_hmc_pulldown_DMdown.txt \\
								$(DIR_PULL_CSAW)/%_hmc_pulldown_noDM_signal.txt \\
								$(DIR_PULL_CSAW)/%_hmc_pulldown_noDM_nosignal.txt'
	rule1 = make_rule_compare_class_bis_module
	rule2 = make_rule_compare_class_pull_module
	class_script = '../../scripts/classify_compare.sh'
} else if (bool_bis_comp && !bool_pull_comp) {
	############################################################
	# NOTE: THIS IS NOT EXPLICITLY SUPPORTED RIGHT NOW
	############################################################
	compare_class_tmps = 'COMPARE_CLASS_CLEAN_TMP := $(patsubst %,$(DIR_BIS_DSS)/%_mc_bisulfite_DMup.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_DSS)/%_mc_bisulfite_DMdown.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_DSS)/%_mc_bisulfite_noDM_signal.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_DSS)/%_mc_bisulfite_noDM_nosignal.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_DSS)/%_hmc_bisulfite_DMup.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_DSS)/%_hmc_bisulfite_DMdown.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_DSS)/%_hmc_bisulfite_noDM_signal.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_BIS_DSS)/%_hmc_bisulfite_noDM_nosignal.txt,$(COMPARE_CLASS_PREFIXES)) \\'
	compare_class_target = '$(DIR_CLASS_COMPARE)/%_compare_classification.bed :	 $(DIR_BIS_DSS)/%_mc_bisulfite_DMup.txt \\
								$(DIR_BIS_DSS)/%_mc_bisulfite_DMdown.txt \\
								$(DIR_BIS_DSS)/%_mc_bisulfite_noDM_signal.txt \\
								$(DIR_BIS_DSS)/%_mc_bisulfite_noDM_nosignal.txt \\
								$(DIR_BIS_DSS)/%_hmc_bisulfite_DMup.txt \\
								$(DIR_BIS_DSS)/%_hmc_bisulfite_DMdown.txt \\
								$(DIR_BIS_DSS)/%_hmc_bisulfite_noDM_signal.txt \\
								$(DIR_BIS_DSS)/%_hmc_bisulfite_noDM_nosignal.txt'
	rule1 = make_rule_compare_class_bis_module
	rule2 = ''
	class_script = '../../scripts/classify_compare.sh'
} else if (!bool_bis_comp && bool_pull_comp) {
	compare_class_tmps = 'COMPARE_CLASS_CLEAN_TMP := $(patsubst %,$(DIR_PULL_CSAW)/%_mc_pulldown_DMup.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_CSAW)/%_mc_pulldown_DMdown.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_CSAW)/%_mc_pulldown_noDM_signal.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_CSAW)/%_mc_pulldown_noDM_nosignal.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_CSAW)/%_hmc_pulldown_DMup.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_CSAW)/%_hmc_pulldown_DMdown.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_CSAW)/%_hmc_pulldown_noDM_signal.txt,$(COMPARE_CLASS_PREFIXES)) \\
								$(patsubst %,$(DIR_PULL_CSAW)/%_hmc_pulldown_noDM_nosignal.txt,$(COMPARE_CLASS_PREFIXES)) \\'
	compare_class_target = '$(DIR_CLASS_COMPARE)/%_compare_classification.bed :	 $(DIR_PULL_CSAW)/%_mc_pulldown_DMup.txt \\
								$(DIR_PULL_CSAW)/%_mc_pulldown_DMdown.txt \\
								$(DIR_PULL_CSAW)/%_mc_pulldown_noDM_signal.txt \\
								$(DIR_PULL_CSAW)/%_mc_pulldown_noDM_nosignal.txt \\
								$(DIR_PULL_CSAW)/%_hmc_pulldown_DMup.txt \\
								$(DIR_PULL_CSAW)/%_hmc_pulldown_DMdown.txt \\
								$(DIR_PULL_CSAW)/%_hmc_pulldown_noDM_signal.txt \\
								$(DIR_PULL_CSAW)/%_hmc_pulldown_noDM_nosignal.txt'
	rule1 = make_rule_compare_class_pull_module
	rule2 = ''
	class_script = '../../scripts/classify_compare.sh'
}

make_rule_class_compare = sprintf('
# Master rule
.PHONY : compare_classification
compare_classification : $(patsubst %%,$(DIR_TRACK)/%%_compare_classification.bb,$(COMPARE_CLASS_PREFIXES)) \\
		$(patsubst %%,$(DIR_RDATA)/%%_compare_classification_annotatr_analysis.RData,$(COMPARE_CLASS_PREFIXES)) \\
		$(patsubst %%,$(DIR_CLASS_COMPARE)/%%_compare_classification.bed,$(COMPARE_CLASS_PREFIXES))

# Rule for compare classification bigBed
$(DIR_TRACK)/%%_compare_classification.bb : $(DIR_CLASS_COMPARE)/%%_compare_classification.bed
	$(PATH_TO_BED2BB) $^ $(CHROM_PATH) $@

# Rule for annotatr of compare classification
$(DIR_RDATA)/%%_compare_classification_annotatr_analysis.RData : $(DIR_CLASS_COMPARE)/%%_compare_classification.bed
	$(PATH_TO_R) ../../scripts/annotatr_annotations.R --file $< --genome $(GENOME) --annot_type compare_class --group1 NULL --group0 NULL

# Classification BED
.PRECIOUS : $(DIR_CLASS_COMPARE)/%%_compare_classification.bed
%s
	bash %s $(PATH_TO_BEDTOOLS) $(PATH_TO_AWK) $(CHROM_PATH) $@ $^

%s
%s

# Clean temporary files that make does not clean up
%s

.PHONY : clean_compare_classification_tmp
clean_compare_classification_tmp :
	rm -f $(COMPARE_CLASS_CLEAN_TMP)',
	compare_class_target, class_script, rule1, rule2, compare_class_tmps)
cat(make_rule_class_compare, file = file_make, sep = '\n', append = TRUE)

#######################################
# PBS script
# bisulfite_compare_q = c(
#	 '#!/bin/bash',
#	 '#### Begin PBS preamble',
#	 '#PBS -N class_compare',
#	 '#PBS -l nodes=1:ppn=4,walltime=24:00:00,pmem=16gb',
#	 '#PBS -A sartor_lab',
#	 '#PBS -q first',
#	 '#PBS -M rcavalca@umich.edu',
#	 '#PBS -m abe',
#	 '#PBS -j oe',
#	 '#PBS -V',
#	 '#### End PBS preamble',
#	 '# Put your job commands after this line',
#	 sprintf('cd ~/latte/mint/projects/%s/',project),
#	 'make -j 4 compare_classification')
# cat(bisulfite_compare_q, file=sprintf('projects/%s/pbs_jobs/classify_compare.q', project), sep='\n')

for(comparison in unique(comparisons$comparison)) {
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
