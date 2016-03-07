################################################################################
# MAKEFILE: bisulfite_compare rules

if(bool_bis_comp) {

	bisulfite_compares = c()
	for(i in 1:nrow(bisulfite_comparisons)) {
		# Establish row variables
		sampleID = bisulfite_comparisons[i,'sampleID']
		humanID = bisulfite_comparisons[i,'humanID']
		pull = bisulfite_comparisons[i,'pulldown']
		bis = bisulfite_comparisons[i,'bisulfite']
		mc_stat = bisulfite_comparisons[i,'mc']
		hmc_stat = bisulfite_comparisons[i,'hmc']
		group = bisulfite_comparisons[i,'group']
		fullHumanID = bisulfite_comparisons[i,'fullHumanID']

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

		# Sorting this ensures that the lower group number is A
		groups = sort(as.integer(unlist(strsplit(group, ','))))

		groupA = subset(bisulfite_samples,
			grepl(groups[1], bisulfite_samples$group) &
			bisulfite_samples$pulldown == pull &
			bisulfite_samples$bisulfite == bis &
			bisulfite_samples$mc == mc_stat &
			bisulfite_samples$hmc == hmc_stat &
			bisulfite_samples$input == 0)
		groupB = subset(bisulfite_samples,
			grepl(groups[2], bisulfite_samples$group) &
			bisulfite_samples$pulldown == pull &
			bisulfite_samples$bisulfite == bis &
			bisulfite_samples$mc == mc_stat &
			bisulfite_samples$hmc == hmc_stat &
			bisulfite_samples$input == 0)

		########################################################################
		# Setup variables to put into the makefile

		var_cytfiles = paste(c(
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt', groupA$fullHumanID), sep=''),
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt', groupB$fullHumanID), sep='')), collapse=',')
		var_sampleids = paste(c(
			groupA$fullHumanID,
			groupB$fullHumanID), collapse=',')
		var_treatment = paste(c(
			rep.int(groups[1],nrow(groupA)),
			rep.int(groups[2],nrow(groupB))), collapse=',')
		var_comparison = fullHumanID

		var_cytfiles_pre = paste(c(
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt', groupA$fullHumanID), sep=''),
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt', groupB$fullHumanID), sep='')), collapse=' ')

		# Targets
		msig_results = sprintf('$(DIR_BIS_MSIG)/%s_methylSig.txt', var_comparison)
		msig_tmp_results = sprintf('$(DIR_BIS_MSIG)/%s_methylSig_tmp.txt', var_comparison)
		msig_bigwig = sprintf('$(DIR_TRACK)/%s_methylSig.bw', var_comparison)

		########################################################################
		# Variables for the makefile

		# Write the bisulfite_compare variables for this comparison
		make_vars_bis_compare = c(
			'################################################################################',
			'# Workflow for bisulfite_compare',
			sprintf('BISULFITE_COMPARE_%s_PREREQS := %s %s', i, msig_results, msig_bigwig),
			sprintf('BISULFITE_COMPARE_%s_CYTFILES := %s', i, var_cytfiles),
			sprintf('BISULFITE_COMPARE_%s_SAMPLEIDS := %s', i, var_sampleids),
			sprintf('BISULFITE_COMPARE_%s_TREATMENT := %s', i, var_treatment),
			sprintf('BISULFITE_COMPARE_%s_COMPARISON := %s', i, var_comparison))
		cat(make_vars_bis_compare, file = file_make, sep='\n', append=T)

		########################################################################
		# Rules for the makefile

		# Write the bisulfite_compare rule for this comparison
		make_rule_bis_compare = c(
			sprintf('.PHONY : bisulfite_compare_%s', i),
			sprintf('bisulfite_compare_%s : $(BISULFITE_COMPARE_%s_PREREQS)', i, i),
			'',
			sprintf('%s : %s', msig_results, var_cytfiles_pre),
			sprintf('	Rscript ../../scripts/process_bisulfite_comparison_run_methylSig.R --project $(PROJECT) --cytfiles $(BISULFITE_COMPARE_%s_CYTFILES) --sampleids $(BISULFITE_COMPARE_%s_SAMPLEIDS) --treatment $(BISULFITE_COMPARE_%s_TREATMENT) --assembly $(GENOME) --pipeline mint --comparison $(BISULFITE_COMPARE_%s_COMPARISON) $(OPTS_METHYLSIG)', i, i, i, i),
			'',
			sprintf('.INTERMEDIATE : %s', msig_tmp_results),
			sprintf('%s : %s', msig_tmp_results, msig_results), # THIS IS CUSTOMIZABLE to p-value or FDR and the threshold
			"	awk -v OFS='\\t' '$$5 < 0.05 {print $$1, $$2, $$3, $$7 }' $^ | sort -T . -k1,1 -k2,2n > $@",
			'',
			sprintf('%s : %s', msig_bigwig, msig_tmp_results),
			'	bedGraphToBigWig $^ $(CHROM_PATH) $@',
			'')
		cat(make_rule_bis_compare, file = file_make, sep='\n', append=T)

		# Track the number of bisulfite compares
		bisulfite_compares = c(bisulfite_compares, sprintf('bisulfite_compare_%s', i))

		# trackDb.txt entry for methylSig result
		trackEntry = c(
		  sprintf('track %s_DM', fullHumanID),
		  sprintf('parent %s_group_comparison', humanID),
		  sprintf('bigDataUrl %s_methylSig.bw', fullHumanID),
		  sprintf('shortLabel %s_DM', fullHumanID),
		  sprintf('longLabel %s_DM_methylSig', fullHumanID),
		  'visibility full',
		  'autoScale on',
		  'alwaysZero on',
		  'at y=0.0 on',
		  'type bigWig',
		  'priority 1.2',
		  ' ')
		cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
	}

	############################################################################
	# Write the master bisulfite_compare rule that will call all created
	make_rule_master_bis_compare = c(
		'.PHONY : bisulfite_compare',
		sprintf('bisulfite_compare : bisulfite_align %s', paste(bisulfite_compares, collapse=' ')),
		'')
	cat(make_rule_master_bis_compare, file = file_make, sep='\n', append=T)

	#######################################
	# PBS script
	bisulfite_compare_q = c(
		'#!/bin/bash',
		'#### Begin PBS preamble',
		'#PBS -N bis_compare',
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
		'make -j 4 bisulfite_compare')
	cat(bisulfite_compare_q, file=sprintf('projects/%s/pbs_jobs/bisulfite_compare.q', project), sep='\n')
}
