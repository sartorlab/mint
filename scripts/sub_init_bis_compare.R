################################################################################
# MAKEFILE: bisulfite_compare rules

if(bool_bis_comp) {

	########################################################################
	# OPTS for config.mk
	config_bis_compare = '################################################################################
# bisulfite_compare configuration options

# Thresholds to use for DMCs or DMRs (above) dss output
# FDR significance level
OPT_DSS_DM_FDR_THRESHOLD = 0.05
# Alternative p-value significance level if no results FDR < 0.25
OPT_DSS_DM_PVAL_THRESHOLD = 0.005
# Desired absolute value of methylation difference
OPT_DSS_DM_DIFF_THRESHOLD = 10
'
	cat(config_bis_compare, file = file_config, sep='\n', append=T)

	bisulfite_compares = c()
	bisulfite_compare_rules = c()
	bisulfite_clean_tmps = c()
	bisulfite_compare_configs = c()

	for(i in 1:nrow(bisulfite_comparisons)) {
		# Establish row variables
		projectID = bisulfite_comparisons[i,'projectID']
		comparison = bisulfite_comparisons[i,'comparison']
		pull = bisulfite_comparisons[i,'pulldown']
		bis = bisulfite_comparisons[i,'bisulfite']
		mc_stat = bisulfite_comparisons[i,'mc']
		hmc_stat = bisulfite_comparisons[i,'hmc']
		input = bisulfite_comparisons[i,'input']
		model = bisulfite_comparisons[i,'model']
		contrast = bisulfite_comparisons[i,'contrast']
		covariates = bisulfite_comparisons[i,'covariates']
		covIsNumeric = bisulfite_comparisons[i,'covIsNumeric']
		groups = bisulfite_comparisons[i,'groups']
		interpretation = bisulfite_comparisons[i,'interpretation']
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

		# Sorting this ensures that the higher group number is A
		# From scripts/dss_run.R, the bigger number is considered
		# the treatment and the smaller is considered the control
		# And meth.diff is always treatment - control
		# group = c('Treatment' = max(treatment),'Control' = min(treatment))
		# So hyper and hypo is with respect to the higher group number
		group_order = order(as.integer(unlist(strsplit(groups, ','))), decreasing = TRUE)
		groups = as.integer(unlist(strsplit(groups, ',')))[group_order]
		interps = unlist(strsplit(interpretation, ','))[rev(group_order)]

		# Create a list
		sample_groups = lapply(bisulfite_samples$group, function(g){
			as.integer(unlist(strsplit(as.character(g), ',')))
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
		groupAname = interps[1]
		groupBname = interps[2]

		if(groupAname == '' || groupBname == '') {
			stop('Groups used for comparisons must be named. Check your _groups.txt file.')
		}

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

		var_comparison = fullHumanID

		# Targets after aggregation of per chromosome dss analyses
		dss_results = sprintf('$(DIR_BIS_DSS)/%s_dss_significant.txt', var_comparison)
		dss_bedgraph = sprintf('$(DIR_BIS_DSS)/%s_dss_for_bigWig.bedGraph', var_comparison)
		dss_bedgraph_sorted = sprintf('$(DIR_BIS_DSS)/%s_dss_for_bigWig_sorted.bedGraph', var_comparison)
		annotatr_bed = sprintf('$(DIR_BIS_DSS)/%s_dss_for_annotatr.txt', var_comparison)
		annotatr_rdata = sprintf('$(DIR_RDATA)/%s_dss_annotatr_analysis.RData', var_comparison)
		dss_bigwig = sprintf('$(DIR_TRACK)/%s_dss.bw', var_comparison)

		########################################################################
		# Setup variables to put into the makefile

		# NOTE: The collapse here is a comma
		var_cytfiles = paste(c(
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_dss.txt', groupA$fullHumanID), sep=''),
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_dss.txt', groupB$fullHumanID), sep='')), collapse=',')
		var_sampleids = paste(c(
			groupA$fullHumanID,
			groupB$fullHumanID), collapse=',')
		var_groups = paste(c(
			rep.int(groups[1],nrow(groupA)),
			rep.int(groups[2],nrow(groupB))), collapse=',')

		# NOTE: The collapse here is a space
		var_cytfiles_pre = paste(c(
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_dss.txt', groupA$fullHumanID), sep=''),
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_dss.txt', groupB$fullHumanID), sep='')), collapse=' ')

		########################################################################
		# Variables for the makefile

		# Write the bisulfite_compare variables for this comparison
		make_vars_bis_compare = c(
			'########################################',
			sprintf('# Workflow for bisulfite_compare_%s', i),
			'',
			sprintf('BISULFITE_COMPARE_%s_PREREQS := %s %s %s', i, dss_results, annotatr_rdata, dss_bigwig),
			sprintf('BISULFITE_COMPARE_%s_CYTFILES := %s', i, var_cytfiles),
			sprintf('BISULFITE_COMPARE_%s_SAMPLEIDS := %s', i, var_sampleids),
			sprintf('BISULFITE_COMPARE_%s_GROUPS := %s', i, var_groups),
			sprintf('BISULFITE_COMPARE_%s_COMPARISON := %s', i, var_comparison),
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
			'# Rule for dss input of bismark CpG report',
			'.INTERMEDIATE : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_for_dss.txt',
			'$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_for_dss.txt : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report.txt.gz',
			'	$(PATH_TO_AWK) -v OFS="\\t" -v MIN_COV=$(OPT_MIN_COV) \'$$4 + $$5 >= MIN_COV {print $$1, $$2, $$3, $$4 + $$5, $$4}\' <(gunzip -c $<) | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@',
			'',
			'# Rule for dss',
			sprintf('%s : %s', dss_results, var_cytfiles_pre),
			sprintf('	$(PATH_TO_R) ../../scripts/dss_run.R --project $(PROJECT) --genome $(GENOME) --files $(BISULFITE_COMPARE_%s_CYTFILES) --samplenames $(BISULFITE_COMPARE_%s_SAMPLEIDS) --model %s --groups $(BISULFITE_COMPARE_%s_GROUPS) --contrast %s --covariates %s --covIsNumeric %s --interpretation %s --outprefix $(BISULFITE_COMPARE_%s_COMPARISON) $(OPTS_DSS_%s)', i, i, model, i, contrast, var_covariates, covIsNumeric, interpretation, i, var_comparison),
			sprintf('%s : %s', dss_bedgraph, dss_results),
			sprintf('%s : %s', annotatr_bed, dss_results))

		make_rule_bis_compare_post = c(
			'',
			'# Rule for annotatr of dss filtered results',
			sprintf('%s : %s', annotatr_rdata, annotatr_bed),
			sprintf('	$(PATH_TO_R) ../../scripts/annotatr_annotations.R --file $< --genome $(GENOME) --annot_type dss --group1 $(GROUP1_NAME_%s) --group0 $(GROUP0_NAME_%s)', i, i),
			'',
			'# Rule to sort bedGraph to match CHROM_PATH sorting',
			sprintf('%s : %s', dss_bedgraph_sorted, dss_bedgraph),
			'	sort -T $(DIR_TMP) -k1,1 -k2,2n $^ > $@',
			'',
			'# Rule for UCSC bigWig of filtered dss_bedgraph_sorted results',
			sprintf('%s : %s', dss_bigwig, dss_bedgraph_sorted),
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

# For dss parameters see http://www.bioconductor.org/packages/DSS
OPTS_DSS_%s = --destrand TRUE --tilewidth 50 --methdiffthreshold $(OPT_DSS_DM_DIFF_THRESHOLD) --FDRthreshold $(OPT_DSS_DM_FDR_THRESHOLD)  --pvalthreshold $(OPT_DSS_DM_PVAL_THRESHOLD) --quiet FALSE
',
			i, i, i, groupAname, i, groupBname, var_comparison)

		bisulfite_compare_configs = c(
			bisulfite_compare_configs,
			config_bis_compare)

		# trackDb.txt entry for dss result
		trackEntry = c(
			sprintf('track %s_DM', fullHumanID),
			sprintf('parent %s_group_comparison', comparison),
			sprintf('bigDataUrl %s_dss.bw', fullHumanID),
			sprintf('shortLabel %s_DM', fullHumanID),
			sprintf('longLabel %s_DM_dss', fullHumanID),
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
