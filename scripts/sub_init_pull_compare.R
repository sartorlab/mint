################################################################################
# MAKEFILE: pulldown_compare rules

if(bool_pull_comp) {

	########################################################################
	# OPTS for config.mk
	config_pull_compare = '################################################################################
# pulldown_compare configuration options

# Thresholds to use for csaw output
# FDR significance level
OPT_CSAW_DM_FDR_THRESHOLD = 0.05
# Alternative p-value significance level if no results FDR < 0.25
OPT_CSAW_DM_PVAL_THRESHOLD = 0.005
	'
	cat(config_pull_compare, file = file_config, sep='\n', append=T)

	# Keep track of the compares for the master make rules
	pulldown_compares = c()
	pulldown_compare_rules = c()
	pulldown_clean_tmps = c()
	pulldown_compare_configs = c()
	for(i in 1:nrow(pulldown_comparisons)) {
		# Establish row variables
		projectID = pulldown_comparisons[i,'projectID']
		comparison = pulldown_comparisons[i,'comparison']
		pull = pulldown_comparisons[i,'pulldown']
		bis = pulldown_comparisons[i,'bisulfite']
		mc_stat = pulldown_comparisons[i,'mc']
		hmc_stat = pulldown_comparisons[i,'hmc']
		input = pulldown_comparisons[i,'input']
		model = pulldown_comparisons[i,'model']
		contrast = pulldown_comparisons[i,'contrast']
		covariates = pulldown_comparisons[i,'covariates']
		covIsNumeric = pulldown_comparisons[i,'covIsNumeric']
		groups = pulldown_comparisons[i,'groups']
		interpretation = pulldown_comparisons[i,'interpretation']
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

		# Sorting this way ensures the higher group number is groupA
		group_order = order(as.integer(unlist(strsplit(groups, ','))), decreasing = TRUE)
		groups = as.integer(unlist(strsplit(groups, ',')))[group_order]
		interps = unlist(strsplit(interpretation, ','))[rev(group_order)]

		# Create a list
		sample_groups = lapply(pulldown_samples$group, function(g){
			as.integer(unlist(strsplit(as.character(g), ',')))
		})

		# Get the correct indices
		groupAindices = sapply(sample_groups, function(sg){
			groups[1] %in% sg
		})
		groupBindices = sapply(sample_groups, function(sg){
			groups[2] %in% sg
		})

		groupA = subset(pulldown_samples,
			groupAindices &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc_stat &
			pulldown_samples$hmc == hmc_stat &
			pulldown_samples$input == 0)
		groupB = subset(pulldown_samples,
			groupBindices &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc_stat &
			pulldown_samples$hmc == hmc_stat &
			pulldown_samples$input == 0)
		groupAname = interps[1]
		groupBname = interps[2]

		# if(groupAname == '' || groupBname == '') {
		# 	stop('Groups used for comparisons must be named. Check your _groups.txt file.')
		# }

		inputGroupA = subset(pulldown_samples,
			groupAindices &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc_stat &
			pulldown_samples$hmc == hmc_stat &
			pulldown_samples$input == 1)
		inputGroupB = subset(pulldown_samples,
			groupBindices &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc_stat &
			pulldown_samples$hmc == hmc_stat &
			pulldown_samples$input == 1)

		########################################################################
		# Deal with covariates
		if(!is.na(covariates)) {
			covariates = unlist(strsplit(covariates, ','))
			var_covariates = c()
			for(j in seq_along(covariates)) {
				covariate_values = c(groupA[, covariates[j]], groupB[, covariates[j]])
				var_covariate = paste(covariates[j], paste(covariate_values, collapse=','), sep=':')
				var_covariates = c(var_covariates, var_covariate)
			}
			var_covariates = paste(var_covariates, collapse=';')
		} else {
			var_covariates = NA
		}

		########################################################################
		# Setup variables to put into the makefile

		# For the csaw call from projects/project/
		var_projectID = projectID
		var_input = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', c(inputGroupA$fullHumanID, inputGroupB$fullHumanID)), sep='', collapse=',')
		var_chip = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', c(groupA$fullHumanID, groupB$fullHumanID)), sep='', collapse=',')
		var_chipnames = paste(sprintf('%s', c(groupA$humanID, groupB$humanID)), sep='', collapse=',')
		var_useinput = input
		var_model = model
		var_groups = paste(sprintf('%s', c(groupA$group, groupB$group)), sep='', collapse=',')
		var_interpretation = interpretation
		var_contrast = contrast
		var_covIsNumeric = covIsNumeric

		# For the prerequisites in the make rule
		var_merged_input_pre = paste(sprintf('$(DIR_PULL_COVERAGES)/%s_coverage_merged.bdg', c(inputGroupA$fullHumanID, inputGroupB$fullHumanID)), sep='', collapse=' ')
		var_input_pre = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', c(inputGroupA$fullHumanID, inputGroupB$fullHumanID)), sep='', collapse=' ')
		var_chip_pre = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', c(groupA$fullHumanID,groupB$fullHumanID)), sep='', collapse=' ')
		var_name = fullHumanID

		# Targets
		input_signal = sprintf('$(DIR_PULL_CSAW)/%s_merged_signal.bed', var_name)
		csaw_significant = sprintf('$(DIR_PULL_CSAW)/%s_csaw_significant.txt', var_name)
		annotatr_txt = sprintf('$(DIR_PULL_CSAW)/%s_csaw_for_annotatr.txt', var_name)
		bigBed_bed = sprintf('$(DIR_PULL_CSAW)/%s_csaw_for_bigBed.bed', var_name)
		annotatr_rdata = sprintf('$(DIR_RDATA)/%s_csaw_annotatr_analysis.RData', var_name)
		bigbed = sprintf('$(DIR_TRACK)/%s_csaw_peaks.bb', var_name)
		simple_bed = sprintf('$(DIR_CLASS_SIMPLE)/%s_csaw_simple_classification.bed', var_name)
		simple_bb = sprintf('$(DIR_TRACK)/%s_csaw_simple_classification.bb', var_name)
		simple_rdata = sprintf('$(DIR_RDATA)/%s_csaw_simple_classification_annotatr_analysis.RData', var_name)

		########################################################################
		# Variables for the makefile

		# Write the pulldown_compare variables for this comparison
		make_var_pull_compare = c(
			'########################################',
			sprintf('# Workflow for pulldown_compare_%s', i),
			'',
			sprintf('PULLDOWN_COMPARE_%s_PREREQS := %s %s %s %s', i, csaw_significant, bigbed, input_signal, annotatr_rdata),
			sprintf('PULLDOWN_COMPARE_SIMPLE_%s_PREREQS := %s %s %s', i, simple_bed, simple_bb, simple_rdata),
			sprintf('PULLDOWN_COMPARE_%s_INPUT := %s', i, var_input),
			sprintf('PULLDOWN_COMPARE_%s_CHIP := %s', i, var_chip),
			sprintf('PULLDOWN_COMPARE_%s_NAME := %s', i, var_name),
			sprintf('PULLDOWN_COMPARE_%s_CLEAN_TMP := %s', i, bigBed_bed))

		# Write the pulldown_compare rule for this comparison
		make_rule_pull_compare = c(
			'',
			'########################################',
			sprintf('.PHONY : pulldown_compare_%s', i),
			sprintf('pulldown_compare_%s : $(PULLDOWN_COMPARE_%s_PREREQS)', i, i),
			'',
			'# Rule for csaw peaks',
			sprintf('%s : %s %s', csaw_significant, var_input_pre, var_chip_pre),
			sprintf('	$(PATH_TO_R) ../../scripts/csaw_run.R --project %s --chipfiles $(PULLDOWN_COMPARE_%s_CHIP) --inputfiles $(PULLDOWN_COMPARE_%s_INPUT) --chipnames %s --useinput %s --model %s --groups %s --contrast %s --covariates %s --covIsNumeric %s --interpretation %s --quiet FALSE --outprefix $(PULLDOWN_COMPARE_%s_NAME) $(OPTS_CSAW_%s)', var_projectID, i, i, var_chipnames, var_useinput, var_model, var_groups, var_contrast, var_covariates, var_covIsNumeric, var_interpretation, i, var_name),
			sprintf('%s : %s', annotatr_txt, csaw_significant),
			sprintf('%s : %s', bigBed_bed, csaw_significant),
			'',
			'# Rule for annotatr of csaw peaks',
			sprintf('%s : %s', annotatr_rdata, annotatr_txt),
			sprintf('	$(PATH_TO_R) ../../scripts/annotatr_annotations.R --file $< --genome $(GENOME) --annot_type csaw --group1 $(CHIP1_NAME_%s) --group0 $(CHIP0_NAME_%s)', i, i),
			'',
			'# Rule to merge input signals from the two groups',
			sprintf('%s : %s', input_signal, var_merged_input_pre),
			'	cat $^ | sort -T $(DIR_TMP) -k1,1 -k2,2n | bedtools merge -d 20 | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@',
			'',
			'# Rule for UCSC bigBed track of csaw peaks',
			sprintf('%s : %s', bigbed, bigBed_bed),
			'	$(PATH_TO_BED2BB) $^ $(CHROM_PATH) $@',
			'',
			'########################################',
			sprintf('.PHONY : pulldown_compare_simple_classification_%s', i),
			sprintf('pulldown_compare_simple_classification_%s : $(PULLDOWN_COMPARE_SIMPLE_%s_PREREQS)', i, i),
			'',
			'# Rule for simple classification of combined csaw peaks',
			sprintf('%s : %s', simple_bed, annotatr_txt),
			sprintf('	$(PATH_TO_R) ../../scripts/classify_simple.R --project $(PROJECT) --inFile $< --outFile $@ --group1 $(CHIP1_NAME_%s) --group0 $(CHIP0_NAME_%s)', i, i),
			'',
			'# Rule for annotatr of csaw simple classifications',
			sprintf('%s : %s', simple_rdata, simple_bed),
			sprintf('	$(PATH_TO_R) ../../scripts/annotatr_annotations.R --file $< --genome $(GENOME) --annot_type simple_pulldown_csaw --group1 $(CHIP1_NAME_%s) --group0 $(CHIP0_NAME_%s)', i, i),
			'',
			'# Rule for UCSC bigBed track',
			sprintf('%s : %s', simple_bb, simple_bed),
			'	$(PATH_TO_BED2BB) $< $(CHROM_PATH) $@',
			'',
			'########################################',
			sprintf('# Rule to delete all temporary files from make pulldown_compare_%s',i),
			sprintf('.PHONY : clean_pulldown_compare_tmp_%s', i),
			sprintf('clean_pulldown_compare_tmp_%s :', i),
			sprintf('	rm -f $(PULLDOWN_COMPARE_%s_CLEAN_TMP)', i))

		# Track all the rules for the pulldown compares
		pulldown_compare_rules = c(
			pulldown_compare_rules,
			make_var_pull_compare,
			make_rule_pull_compare)

		########################################################################
		# OPTS for config.mk
		config_pull_compare = sprintf('########################################
# pulldown_compare_%s configuration options

# Informative names for group1 and group0
# GROUP1_NAME should be the group with higher group number in the project annotation
# file and GROUP0_NAME should be the group with the lower group number
# If unsure, check the "workflow for pulldown_compare_%s" section of the makefile
CHIP1_NAME_%s := %s
CHIP0_NAME_%s := %s

# For csaw parameters see http://www.bioconductor.org/packages/csaw
OPTS_CSAW_%s = --fraglength 110 --winwidth 100 --winspacing 50 --prior.count 5 --chipfold 2 --mergewithin 500 --maxmerged 2000 --FDRthreshold $(OPT_CSAW_DM_FDR_THRESHOLD) --pvalthreshold $(OPT_CSAW_DM_PVAL_THRESHOLD)
',
			i, i, i, groupAname, i, groupBname, var_name)

		# Track all the configs for the pulldown compares
		pulldown_compare_configs = c(
			pulldown_compare_configs,
			config_pull_compare)

		# trackDb.txt entry for csaw output
		trackEntry = c(
			sprintf('track %s_DM', fullHumanID),
			sprintf('parent %s_group_comparison', comparison),
			sprintf('bigDataUrl %s_csaw_peaks.bb', fullHumanID),
			sprintf('shortLabel %s_DM', fullHumanID),
			sprintf('longLabel %s_DM_csaw_peaks', fullHumanID),
			'visibility pack',
			'itemRgb on',
			'type bigBed 9 .',
			'priority 1.3',
			' ')
		cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

		# Track the pulldown compares
		pulldown_compares = c(pulldown_compares, sprintf('pulldown_compare_%s', i), sprintf('pulldown_compare_simple_classification_%s', i))
		pulldown_clean_tmps = c(pulldown_clean_tmps, sprintf('clean_pulldown_compare_tmp_%s', i))
	}

	############################################################################
	# Write the master pulldown_compare rule that will call all created
	make_rule_master_pull_compare = c(
		'################################################################################',
		'# Workflow for pulldown_compare',
		'',
		'########################################',
		'',
		'.PHONY : pulldown_compare',
		sprintf('pulldown_compare : %s', paste(pulldown_compares, collapse=' ')),
		'',
		pulldown_compare_rules,
		'',
		'########################################',
		'# Rule to delete all temporary files from make pulldown_compare',
		'.PHONY : clean_pulldown_compare_tmp',
		sprintf('clean_pulldown_compare_tmp : %s', paste(pulldown_clean_tmps, collapse=' ')),
		'',
		'################################################################################',
		'')
	cat(make_rule_master_pull_compare, file = file_make, sep='\n', append=T)

	config_master_pull_compare = c(
		pulldown_compare_configs
	)
	cat(config_master_pull_compare, file = file_config, sep='\n', append=T)

	#######################################
	# PBS script
	# pulldown_compare_q = c(
	#	 '#!/bin/bash',
	#	 '#### Begin PBS preamble',
	#	 '#PBS -N pull_compare',
	#	 '#PBS -l nodes=1:ppn=8,pmem=16gb,walltime=24:00:00',
	#	 '#PBS -A sartor_lab',
	#	 '#PBS -q first',
	#	 '#PBS -M rcavalca@umich.edu',
	#	 '#PBS -m abe',
	#	 '#PBS -j oe',
	#	 '#PBS -V',
	#	 '#### End PBS preamble',
	#	 '# Put your job commands after this line',
	#	 sprintf('cd ~/latte/mint/projects/%s/',project),
	#	 'make pulldown_compare')
	# cat(pulldown_compare_q, file=sprintf('projects/%s/pbs_jobs/pulldown_compare.q', project), sep='\n')

}
