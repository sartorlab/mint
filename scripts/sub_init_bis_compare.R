################################################################################
# MAKEFILE: bisulfite_compare rules

if(bool_bis_comp) {

	########################################################################
	# OPTS for config.mk
	config_bis_compare = '################################################################################
# bisulfite_compare configuration options

# DMC for CpG resolution, and DMR for region resolution (window size parameter
# used in the methylSig options below).
OPT_DM_TYPE = DMR

# Thresholds to use for DMCs or DMRs (above) methylSig output
# FDR significance level
OPT_MSIG_DM_FDR_THRESHOLD = 0.05
# Desired absolute value of methylation difference
OPT_MSIG_DM_DIFF_THRESHOLD = 10
'
	cat(config_bis_compare, file = file_config, sep='\n', append=T)

	bisulfite_compares = c()
	bisulfite_compare_rules = c()
	bisulfite_clean_tmps = c()
	bisulfite_compare_configs = c()

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
		# From scripts/methylSig_run.R, the bigger number is considered
		# the treatment and the smaller is considered the control
		# And meth.diff is always treatment - control
		# group = c('Treatment' = max(treatment),'Control' = min(treatment))
		# So hyper and hypo is with respect to the higher group number
		groups = sort(as.integer(unlist(strsplit(group, ','))), decreasing=TRUE)

		# Create a list
		sample_groups = lapply(bisulfite_samples$group, function(g){
			as.integer(unlist(strsplit(g, ',')))
		})

		# Get the correct indices
		groupAindices = sapply(sample_groups, function(sg){
			groups[1] %in% sg
		})
		groupBindices = sapply(sample_groups, function(sg){
			groups[2] %in% sg
		})

		groupA = subset(bisulfite_samples,
			groupAindices &
			bisulfite_samples$pulldown == pull &
			bisulfite_samples$bisulfite == bis &
			bisulfite_samples$mc == mc_stat &
			bisulfite_samples$hmc == hmc_stat &
			bisulfite_samples$input == 0)
		groupB = subset(bisulfite_samples,
			groupBindices &
			bisulfite_samples$pulldown == pull &
			bisulfite_samples$bisulfite == bis &
			bisulfite_samples$mc == mc_stat &
			bisulfite_samples$hmc == hmc_stat &
			bisulfite_samples$input == 0)
		groupAname = subset(group_names, group == groups[1])$name
		groupBname = subset(group_names, group == groups[2])$name

		if(groupAname == '' || groupBname == '') {
			stop('Groups used for comparisons must be named. Check your _groups.txt file.')
		}

		var_comparison = fullHumanID

		# Targets after aggregation of per chromosome methylSig analyses
		msig_results = sprintf('$(DIR_BIS_MSIG)/%s_$(OPT_DM_TYPE)_methylSig.txt', var_comparison)
		msig_tmp_results = sprintf('$(DIR_BIS_MSIG)/%s_$(OPT_DM_TYPE)_methylSig_tmp.txt', var_comparison)
		annotatr_bed = sprintf('$(DIR_BIS_MSIG)/%s_$(OPT_DM_TYPE)_methylSig_for_annotatr.txt', var_comparison)
		annotatr_rdata = sprintf('$(DIR_RDATA)/%s_$(OPT_DM_TYPE)_methylSig_annotatr_analysis.RData', var_comparison)
		msig_bigwig = sprintf('$(DIR_TRACK)/%s_methylSig.bw', var_comparison)

		########################################################################
		# Setup variables to put into the makefile

		# NOTE: The collapse here is a comma
		var_cytfiles = paste(c(
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_methylSig.txt', groupA$fullHumanID), sep=''),
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_methylSig.txt', groupB$fullHumanID), sep='')), collapse=',')
		var_sampleids = paste(c(
			groupA$fullHumanID,
			groupB$fullHumanID), collapse=',')
		var_treatment = paste(c(
			rep.int(groups[1],nrow(groupA)),
			rep.int(groups[2],nrow(groupB))), collapse=',')

		# NOTE: The collapse here is a space
		var_cytfiles_pre = paste(c(
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_methylSig.txt', groupA$fullHumanID), sep=''),
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_methylSig.txt', groupB$fullHumanID), sep='')), collapse=' ')

		########################################################################
		# Variables for the makefile

		# Write the bisulfite_compare variables for this comparison
		make_vars_bis_compare = c(
			'########################################',
			sprintf('# Workflow for bisulfite_compare_%s', i),
			'',
			sprintf('BISULFITE_COMPARE_%s_PREREQS := %s %s %s', i, msig_results, annotatr_rdata, msig_bigwig),
			sprintf('BISULFITE_COMPARE_%s_CYTFILES := %s', i, var_cytfiles),
			sprintf('BISULFITE_COMPARE_%s_SAMPLEIDS := %s', i, var_sampleids),
			sprintf('BISULFITE_COMPARE_%s_TREATMENT := %s', i, var_treatment),
			sprintf('BISULFITE_COMPARE_%s_COMPARISON := %s_$(OPT_DM_TYPE)_methylSig', i, var_comparison),
			sprintf('BISULFITE_COMPARE_%s_CLEAN_TMP := %s', i, var_cytfiles_pre))

		########################################################################
		# Rules for the makefile

		# Write the bisulfite_compare rule for this comparison
		make_rule_bis_compare = c(
			'',
			'########################################',
			sprintf('.PHONY : bisulfite_compare_%s', i),
			sprintf('bisulfite_compare_%s : $(BISULFITE_COMPARE_%s_PREREQS)', i, i),
			'',
			'# Rule for methylSig input of bismark CpG report',
			'.INTERMEDIATE : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_for_methylSig.txt',
			'$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_for_methylSig.txt : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report.txt.gz',
			'	$(PATH_TO_AWK) -f ../../scripts/extractor_to_methylSig.awk <(gunzip -c $<) | sort -T $(DIR_TMP) -k2,2 -k3,3n > $@',
			'',
			'# Rule for methylSig',
			sprintf('%s : %s', msig_results, var_cytfiles_pre),
			sprintf('	$(PATH_TO_R) ../../scripts/methylSig_run.R --project $(PROJECT) --cytfiles $(BISULFITE_COMPARE_%s_CYTFILES) --sampleids $(BISULFITE_COMPARE_%s_SAMPLEIDS) --treatment $(BISULFITE_COMPARE_%s_TREATMENT) --assembly $(GENOME) --pipeline mint --outprefix $(BISULFITE_COMPARE_%s_COMPARISON) $(OPTS_METHYLSIG_%s)', i, i, i, i, var_comparison))

		make_rule_bis_compare_post = c(
			'',
			'# Rule for methylSig filtering for annotatr and bigWig',
			'# NOTE: Files used for methylSig are 1-based start and end',
			'# subtract 1 from start in order for annotatr to properly interpret it',
			sprintf('.INTERMEDIATE : %s', msig_tmp_results),
			sprintf('%s : %s', msig_tmp_results, msig_results), # THIS IS CUSTOMIZABLE
			"	$(PATH_TO_AWK) -v OFS='\\t' -v FDR=$(OPT_MSIG_DM_FDR_THRESHOLD) -v DIFF=$(OPT_MSIG_DM_DIFF_THRESHOLD) '$$6 < FDR && sqrt($$7^2) > DIFF { print $$1, $$2 - 1, $$3, $$7 }' $^ | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@",
			'',
			'# Rule for annotatr input of methylSig filtered results',
			'# NOTE: Files used for methylSig are 1-based start and end',
			'# subtract 1 from start (in awk script) in order for annotatr to properly interpret it',
			sprintf('.INTERMEDIATE : %s', annotatr_bed),
			sprintf('%s : %s', annotatr_bed, msig_results),
			sprintf('	$(PATH_TO_AWK) -v FDR=$(OPT_MSIG_DM_FDR_THRESHOLD) -v DIFF=$(OPT_MSIG_DM_DIFF_THRESHOLD) -v GROUP1=$(GROUP1_NAME_%s) -v GROUP0=$(GROUP0_NAME_%s) -f ../../scripts/methylSig_to_annotatr.awk $< > $@', i, i),
			'',
			'# Rule for annotatr of methylSig filtered results',
			sprintf('%s : %s', annotatr_rdata, annotatr_bed),
			sprintf('	$(PATH_TO_R) ../../scripts/annotatr_annotations.R --file $< --genome $(GENOME) --annot_type methylSig --group1 $(GROUP1_NAME_%s) --group0 $(GROUP0_NAME_%s)', i, i),
			'',
			'# Rule for UCSC bigWig of filtered methylSig results',
			sprintf('%s : %s', msig_bigwig, msig_tmp_results),
			'	$(PATH_TO_BDG2BW) $^ $(CHROM_PATH) $@',
			'',
			'########################################',
			sprintf('# Rule to delete all temporary files from make bisulfite_compare_%s',i),
			sprintf('BISULFITE_COMPARE_%s_CLEAN_TMP := %s', i, annotatr_bed),
			'',
			sprintf('.PHONY : clean_bisulfite_compare_tmp_%s', i),
			sprintf('clean_bisulfite_compare_tmp_%s :
				rm -f $(BISULFITE_COMPARE_%s_CLEAN_TMP)', i, i),
			'')

		# Track all the rules for the bisulfite compares
		bisulfite_compare_rules = c(
			bisulfite_compare_rules,
			make_vars_bis_compare,
			make_rule_bis_compare,
			bisulfite_compare_rules,
			make_rule_bis_compare_post)

		########################################################################
		# OPTS for config.mk
		config_bis_compare = sprintf('########################################
# bisulfite_compare_%s configuration options

# Informative names for group1 and group0
# GROUP1_NAME should be the group with higher group number in the project annotation
# file and GROUP0_NAME should be the group with the lower group number
# If unsure, check the "workflow for bisulfite_compare_%s" section of the makefile
GROUP1_NAME_%s := %s
GROUP0_NAME_%s := %s

# For methylSig parameters, in R, after installing methylSig do:
# ?methylSig::methylSigReadData and ?methylSig::methylSigCalc
OPTS_METHYLSIG_%s = --context CpG --resolution base --destranded TRUE --maxcount 500 --mincount 5 --filterSNPs TRUE --dmtype $(OPT_DM_TYPE) --winsize.tile 50 --dispersion both --local.disp FALSE --winsize.disp 200 --local.meth FALSE --winsize.meth 200 --minpergroup 2,2 --T.approx TRUE --ncores 1 --quiet FALSE
',
			i, i, i, groupAname, i, groupBname, var_comparison)

		bisulfite_compare_configs = c(
			bisulfite_compare_configs,
			config_bis_compare)

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

		# Track the number of bisulfite compares and bisulfite clean tmps
		bisulfite_compares = c(bisulfite_compares, sprintf('bisulfite_compare_%s', i))
		bisulfite_clean_tmps = c(bisulfite_clean_tmps, sprintf('clean_bisulfite_compare_tmp_%s', i))
	}

	############################################################################
	# Write the master bisulfite_compare rule that will call all created
	make_rule_master_bis_compare = c(
		'################################################################################',
		'# Workflow for bisulfite_compare',
		'',
		'########################################',
		'',
		'.PHONY : bisulfite_compare',
		sprintf('bisulfite_compare : %s', paste(bisulfite_compares, collapse=' ')),
		'',
		bisulfite_compare_rules,
		'',
		'########################################',
		'# Rule to delete all temporary files from make bisulfite_compare',
		'.PHONY : clean_bisulfite_compare_tmp',
		sprintf('clean_bisulfite_compare_tmp : %s', bisulfite_clean_tmps),
		'',
		'################################################################################',
		'',
		'')
	cat(make_rule_master_bis_compare, file = file_make, sep='\n', append=T)

	config_master_bis_compare = c(
		bisulfite_compare_configs
	)
	cat(config_master_bis_compare, file = file_config, sep='\n', append=T)

	#######################################
	# PBS script
	# bisulfite_compare_q = c(
	#	 '#!/bin/bash',
	#	 '#### Begin PBS preamble',
	#	 '#PBS -N bis_compare',
	#	 '#PBS -l nodes=1:ppn=4,walltime=24:00:00,pmem=8gb',
	#	 '#PBS -A sartor_lab',
	#	 '#PBS -q first',
	#	 '#PBS -M rcavalca@umich.edu',
	#	 '#PBS -m abe',
	#	 '#PBS -j oe',
	#	 '#PBS -V',
	#	 '#### End PBS preamble',
	#	 '# Put your job commands after this line',
	#	 sprintf('cd ~/latte/mint/projects/%s/',project),
	#	 'make -j 4 bisulfite_compare')
	# cat(bisulfite_compare_q, file=sprintf('projects/%s/pbs_jobs/bisulfite_compare.q', project), sep='\n')
}
