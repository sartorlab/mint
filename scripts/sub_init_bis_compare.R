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
OPT_MSIG_DM_DIFF_THRESHOLD = 10'
	cat(config_bis_compare, file = file_config, sep='\n', append=T)

	bisulfite_compares = c()
	bisulfite_clean_tmps = c()

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

		########################################################################
		# Setup variables to put into the makefile

		var_cytfiles = paste(c(
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_methylSig.txt', groupA$fullHumanID), sep=''),
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_methylSig.txt', groupB$fullHumanID), sep='')), collapse=',')
		var_sampleids = paste(c(
			groupA$fullHumanID,
			groupB$fullHumanID), collapse=',')
		var_treatment = paste(c(
			rep.int(groups[1],nrow(groupA)),
			rep.int(groups[2],nrow(groupB))), collapse=',')
		var_comparison = fullHumanID

		var_cytfiles_pre = paste(c(
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_methylSig.txt', groupA$fullHumanID), sep=''),
			paste(sprintf('$(DIR_BIS_BISMARK)/%s_trimmed_bismark_bt2.CpG_report_for_methylSig.txt', groupB$fullHumanID), sep='')), collapse=' ')

		# Targets
		msig_results = sprintf('$(DIR_BIS_MSIG)/%s_$(OPT_DM_TYPE)_methylSig.txt', var_comparison)
		msig_tmp_results = sprintf('$(DIR_BIS_MSIG)/%s_$(OPT_DM_TYPE)_methylSig_tmp.txt', var_comparison)
		annotatr_bed = sprintf('$(DIR_BIS_MSIG)/%s_$(OPT_DM_TYPE)_methylSig_for_annotatr.txt', var_comparison)
		annotatr_rdata = sprintf('$(DIR_RDATA)/%s_$(OPT_DM_TYPE)_methylSig_annotatr_analysis.RData', var_comparison)
		msig_bigwig = sprintf('$(DIR_TRACK)/%s_methylSig.bw', var_comparison)

		########################################################################
		# Variables for the makefile

		# Write the bisulfite_compare variables for this comparison
		make_vars_bis_compare = c(
			'################################################################################',
			'# Workflow for bisulfite_compare',
			sprintf('BISULFITE_COMPARE_%s_PREREQS := %s %s %s', i, msig_results, msig_bigwig, annotatr_rdata),
			sprintf('BISULFITE_COMPARE_%s_CYTFILES := %s', i, var_cytfiles),
			sprintf('BISULFITE_COMPARE_%s_SAMPLEIDS := %s', i, var_sampleids),
			sprintf('BISULFITE_COMPARE_%s_TREATMENT := %s', i, var_treatment),
			sprintf('BISULFITE_COMPARE_%s_COMPARISON := %s_$(OPT_DM_TYPE)_methylSig', i, var_comparison),
			sprintf('BISULFITE_COMPARE_%s_CLEAN_TMP := %s %s %s', i, msig_tmp_results, annotatr_bed, var_cytfiles_pre))
		cat(make_vars_bis_compare, file = file_make, sep='\n', append=T)

		########################################################################
		# Rules for the makefile

		# Write the bisulfite_compare rule for this comparison
		make_rule_bis_compare = c(
			sprintf('.PHONY : bisulfite_compare_%s', i),
			sprintf('bisulfite_compare_%s : bisulfite_align $(BISULFITE_COMPARE_%s_PREREQS)', i, i),
			'',
'# Rule for methylSig input
.INTERMEDIATE : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_for_methylSig.txt
$(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report_for_methylSig.txt : $(DIR_BIS_BISMARK)/%_trimmed_bismark_bt2.CpG_report.txt.gz
	$(PATH_TO_AWK) -f ../../scripts/extractor_to_methylSig.awk <(gunzip -c $<) | sort -T $(DIR_TMP) -k2,2 -k3,3n > $@',
			sprintf('%s : %s', msig_results, var_cytfiles_pre),
			sprintf('	$(PATH_TO_R) ../../scripts/methylSig_run.R --project $(PROJECT) --cytfiles $(BISULFITE_COMPARE_%s_CYTFILES) --sampleids $(BISULFITE_COMPARE_%s_SAMPLEIDS) --treatment $(BISULFITE_COMPARE_%s_TREATMENT) --assembly $(GENOME) --pipeline mint --outprefix $(BISULFITE_COMPARE_%s_COMPARISON) $(OPTS_METHYLSIG_%s)', i, i, i, i, var_comparison),
			'',
			sprintf('.INTERMEDIATE : %s', msig_tmp_results),
			sprintf('%s : %s', msig_tmp_results, msig_results), # THIS IS CUSTOMIZABLE
			"	$(PATH_TO_AWK) -v OFS='\\t' -v FDR=$(OPT_MSIG_DM_FDR_THRESHOLD) -v DIFF=$(OPT_MSIG_DM_DIFF_THRESHOLD) 'NR > 1 && $$6 < FDR && sqrt($$7^2) > DIFF { print $$1, $$2, $$3, $$7 }' $^ | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@",
			'',
			sprintf('.INTERMEDIATE : %s', annotatr_bed),
			sprintf('%s : %s', annotatr_bed, msig_results),
			'	$(PATH_TO_AWK) -v FDR=$(OPT_MSIG_DM_FDR_THRESHOLD) -v DIFF=$(OPT_MSIG_DM_DIFF_THRESHOLD) -f ../../scripts/methylSig_to_annotatr.awk $< > $@',
			'',
			sprintf('%s : %s', annotatr_rdata, annotatr_bed),
			'	$(PATH_TO_R) ../../scripts/annotatr_bisulfite.R --file $< --genome $(GENOME)',
			'',
			sprintf('%s : %s', msig_bigwig, msig_tmp_results),
			'	$(PATH_TO_BDG2BW) $^ $(CHROM_PATH) $@',
			'',
			sprintf('.PHONY : clean_bisulfite_compare_tmp_%s', i),
			sprintf('clean_bisulfite_compare_tmp_%s :
				rm -f $(BISULFITE_COMPARE_%s_CLEAN_TMP)', i, i),
			'')
		cat(make_rule_bis_compare, file = file_make, sep='\n', append=T)

		# Track the number of bisulfite compares and bisulfite clean tmps
		bisulfite_compares = c(bisulfite_compares, sprintf('bisulfite_compare_%s', i))
		bisulfite_clean_tmps = c(bisulfite_clean_tmps, sprintf('clean_bisulfite_compare_tmp_%s', i))

		########################################################################
		# OPTS for config.mk
		config_bis_compare = sprintf('# See ?methylSig::methylSigReadData and ?methylSig::methylSigCalc after installing methylSig in R for parameter information
OPTS_METHYLSIG_%s = --context CpG --resolution base --destranded TRUE --maxcount 500 --mincount 5 --filterSNPs TRUE --dmtype $(OPT_DM_TYPE) --winsize.tile 50 --dispersion both --local.disp FALSE --winsize.disp 200 --local.meth FALSE --winsize.meth 200 --minpergroup 2,2 --T.approx TRUE --ncores 4 --quiet FALSE
',
			var_comparison)
		cat(config_bis_compare, file = file_config, sep='\n', append=T)

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
		'',
		'.PHONY : clean_bisulfite_compare_tmp',
		sprintf('clean_bisulfite_compare_tmp : %s', paste(bisulfite_clean_tmps, collapse=' ')),
		'')
	cat(make_rule_master_bis_compare, file = file_make, sep='\n', append=T)

	#######################################
	# PBS script
	bisulfite_compare_q = c(
		'#!/bin/bash',
		'#### Begin PBS preamble',
		'#PBS -N bis_compare',
		'#PBS -l nodes=1:ppn=4,walltime=24:00:00,pmem=8gb',
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
