################################################################################
# MAKEFILE: pulldown_compare rules

if(bool_pull_comp) {
	# Keep track of the compares for the master make rules
	pulldown_compares = c()
	pulldown_clean_tmps = c()
	for(i in 1:nrow(pulldown_comparisons)) {
		# Establish row variables
		projectID = pulldown_comparisons[i,'projectID']
		sampleID = pulldown_comparisons[i,'sampleID']
		humanID = pulldown_comparisons[i,'humanID']
		pull = pulldown_comparisons[i,'pulldown']
		bis = pulldown_comparisons[i,'bisulfite']
		mc_stat = pulldown_comparisons[i,'mc']
		hmc_stat = pulldown_comparisons[i,'hmc']
		input = pulldown_comparisons[i,'input']
		group = pulldown_comparisons[i,'group']
		fullHumanID = pulldown_comparisons[i,'fullHumanID']

		if( pull == 1 ) {
		  platform = "pulldown"
		} else if ( bis == 1 ) {
		  platform = "bisulfite"
		}

		if( mc_stat == 1 && hmc_stat == 1 ) {
		  mark = "mc_hmc"
		} else if ( mc_stat == 1 && hmc_stat == 0 ) {
		  mark = "mc"
		} else if ( mc_stat == 0 && hmc_stat == 1 ) {
		  mark = "hmc"
		}

		# Sorting this way ensures the higher group number is groupB
		# NOTE: This makes the PePr DM test match that of methylSig,
		# where control is the lower number (often 0) and the treatment
		# is the higher number (often 1). PePr up peaks mean up in chip1/input1
		# with respect to chip2/input2.
		groups = sort(as.integer(unlist(strsplit(group, ','))), decreasing = TRUE)

		groupA = subset(pulldown_samples,
			grepl(groups[1], pulldown_samples$group) &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc_stat &
			pulldown_samples$hmc == hmc_stat &
			pulldown_samples$input == 0)
		groupB = subset(pulldown_samples,
			grepl(groups[2], pulldown_samples$group) &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc_stat &
			pulldown_samples$hmc == hmc_stat &
			pulldown_samples$input == 0)

		inputGroupA = subset(pulldown_samples,
			grepl(groups[1], pulldown_samples$group) &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc_stat &
			pulldown_samples$hmc == hmc_stat &
			pulldown_samples$input == 1)
		inputGroupB = subset(pulldown_samples,
			grepl(groups[2], pulldown_samples$group) &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc_stat &
			pulldown_samples$hmc == hmc_stat &
			pulldown_samples$input == 1)

		########################################################################
		# Setup variables to put into the makefile

		# For the PePr call from projects/project/pepr_peaks/
		var_input1 = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', inputGroupA$fullHumanID), sep='', collapse=',')
		var_input2 = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', inputGroupB$fullHumanID), sep='', collapse=',')
		var_chip1 = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', groupA$fullHumanID), sep='', collapse=',')
		var_chip2 = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', groupB$fullHumanID), sep='', collapse=',')

		# For the prerequisites in the make rule
		var_merged_input1_pre = paste(sprintf('$(DIR_PULL_COVERAGES)/%s_merged_coverage.bdg', inputGroupA$fullHumanID), sep='', collapse=' ')
		var_merged_input2_pre = paste(sprintf('$(DIR_PULL_COVERAGES)/%s_merged_coverage.bdg', inputGroupB$fullHumanID), sep='', collapse=' ')
		var_input1_pre = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', inputGroupA$fullHumanID), sep='', collapse=' ')
		var_input2_pre = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', inputGroupB$fullHumanID), sep='', collapse=' ')
		var_chip1_pre = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', groupA$fullHumanID), sep='', collapse=' ')
		var_chip2_pre = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', groupB$fullHumanID), sep='', collapse=' ')
		var_name = fullHumanID

		# Targets
		input_signal = sprintf('$(DIR_PULL_PEPR)/%s_merged_signal.bed', var_name)
		chip1_bed = sprintf('$(DIR_PULL_PEPR)/%s__PePr_chip1_peaks.bed', var_name)
		chip2_bed = sprintf('$(DIR_PULL_PEPR)/%s__PePr_chip2_peaks.bed', var_name)
		combined_bed = sprintf('$(DIR_PULL_PEPR)/%s_PePr_combined.bed', var_name)
		annotatr_bed = sprintf('$(DIR_PULL_PEPR)/%s_PePr_for_annotatr.txt', var_name)
		annotatr_png = sprintf('$(DIR_SUM_FIGURES)/%s_PePr_counts.png', var_name)
		bigbed = sprintf('$(DIR_TRACK)/%s_PePr_peaks.bb', var_name)

		########################################################################
		# Variables for the makefile

		# Write the pulldown_compare variables for this comparison
		make_var_pull_compare = c(
			'################################################################################',
			'# Workflow for pulldown_compare',
			sprintf('PULLDOWN_COMPARE_%s_PREREQS := %s %s %s %s', i, chip1_bed, bigbed, input_signal, annotatr_png),
			sprintf('PULLDOWN_COMPARE_%s_INPUT1 := %s', i, var_input1),
			sprintf('PULLDOWN_COMPARE_%s_INPUT2 := %s', i, var_input2),
			sprintf('PULLDOWN_COMPARE_%s_CHIP1 := %s', i, var_chip1),
			sprintf('PULLDOWN_COMPARE_%s_CHIP2 := %s', i, var_chip2),
			sprintf('PULLDOWN_COMPARE_%s_NAME := %s', i, var_name),
			sprintf('PULLDOWN_COMPARE_%s_CLEAN_TMP := %s', i, annotatr_bed))
		cat(make_var_pull_compare, file = file_make, sep='\n', append=T)

		# Write the pulldown_compare rule for this comparison
		make_rule_pull_compare = c(
			sprintf('.PHONY : pulldown_compare_%s', i),
			sprintf('pulldown_compare_%s : $(PULLDOWN_COMPARE_%s_PREREQS)', i, i),
			'',
			sprintf('%s : %s %s %s %s', chip1_bed, var_input1_pre, var_input2_pre, var_chip1_pre, var_chip2_pre),
			sprintf('	$(PATH_TO_PEPR) --input1=$(PULLDOWN_COMPARE_%s_INPUT1) --input2=$(PULLDOWN_COMPARE_%s_INPUT2) --chip1=$(PULLDOWN_COMPARE_%s_CHIP1) --chip2=$(PULLDOWN_COMPARE_%s_CHIP2) --name=$(PULLDOWN_COMPARE_%s_NAME) --output-directory=$(DIR_PULL_PEPR) $(OPTS_PEPR_%s)', i, i, i, i, i, var_name),
			sprintf('%s : %s', chip2_bed, chip1_bed),
			'',
			sprintf('%s : %s %s', combined_bed, chip1_bed, chip2_bed),
			'	bash ../../scripts/pepr_combine.sh $(word 1,$^) $(word 2,$^) $@',
			'',
			sprintf('.INTERMEDIATE : %s', annotatr_bed),
			sprintf('%s : %s', annotatr_bed, combined_bed),
			'	cut -f 1-4 $< > $@',
			'',
			sprintf('%s : %s', annotatr_png, annotatr_bed),
			'	$(PATH_TO_R) ../../scripts/annotatr_classification.R --file $< --genome $(GENOME)',
			'',
			sprintf('%s : %s %s', input_signal, var_merged_input1_pre, var_merged_input2_pre),
			'	cat $^ | sort -T . -k1,1 -k2,2n | bedtools merge -d 20 | sort -T . -k1,1 -k2,2n > $@',
			'',
			sprintf('%s : %s', bigbed, combined_bed),
			'	$(PATH_TO_BDG2BB) $^ $(CHROM_PATH) $@',
			'',
			sprintf('.PHONY : clean_pulldown_compare_tmp_%s', i),
			sprintf('clean_pulldown_compare_tmp_%s :
				rm -f $(PULLDOWN_COMPARE_%s_CLEAN_TMP)', i, i),
			'')
		cat(make_rule_pull_compare, file = file_make, sep='\n', append=T)

		# Track the number of pulldown compares
		pulldown_compares = c(pulldown_compares, sprintf('pulldown_compare_%s', i))
		pulldown_clean_tmps = c(pulldown_clean_tmps, sprintf('clean_pulldown_compare_tmp_%s', i))

		########################################################################
		# OPTS for config.mk
		config_pull_compare = sprintf(
			'OPTS_PEPR_%s = --file-format=bam --peaktype=sharp --diff --threshold=1e-05 --num-processors=1
			',
			var_name)
		cat(config_pull_compare, file = file_config, sep='\n', append=T)

		# trackDb.txt entry for PePr output
		trackEntry = c(
		  sprintf('track %s_DM', fullHumanID),
		  sprintf('parent %s_group_comparison', humanID),
		  sprintf('bigDataUrl %s_PePr_peaks.bb', fullHumanID),
		  sprintf('shortLabel %s_DM', fullHumanID),
		  sprintf('longLabel %s_DM_PePr_peaks', fullHumanID),
		  'visibility pack',
		  'itemRgb on',
		  'type bigBed 9 .',
		  'priority 1.3',
		  ' ')
		cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
	}

	############################################################################
	# Write the master pulldown_compare rule that will call all created
	make_rule_master_pull_compare = c(
		'.PHONY : pulldown_compare',
		sprintf('pulldown_compare : pulldown_align %s', paste(pulldown_compares, collapse=' ')),
		'',
		'.PHONY : clean_pulldown_compare_tmp',
		sprintf('clean_pulldown_compare_tmp : %s', paste(pulldown_clean_tmps, collapse=' ')),
		'')
	cat(make_rule_master_pull_compare, file = file_make, sep='\n', append=T)

	#######################################
	# PBS script
	pulldown_compare_q = c(
		'#!/bin/bash',
		'#### Begin PBS preamble',
		'#PBS -N pull_compare',
		'#PBS -l nodes=1:ppn=8,walltime=24:00:00,pmem=8gb',
		'#PBS -A sartor_lab',
		'#PBS -q first',
		'#PBS -M rcavalca@umich.edu',
		'#PBS -m abe',
		'#PBS -j oe',
		'#PBS -V',
		'#### End PBS preamble',
		'# Put your job commands after this line',
		sprintf('cd ~/latte/mint/projects/%s/',project),
		'make pulldown_compare')
	cat(pulldown_compare_q, file=sprintf('projects/%s/pbs_jobs/pulldown_compare.q', project), sep='\n')

}
