################################################################################
# MAKEFILE: pulldown_compare rules

if(bool_pull_comp) {
    # Keep track of the compares for the master make rules
    pulldown_compares = c()
    pulldown_compare_rules = c()
    pulldown_clean_tmps = c()
    pulldown_compare_configs = c()
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

        # Create a list
        sample_groups = lapply(pulldown_samples$group, function(g){
            as.integer(unlist(strsplit(g, ',')))
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
        # Setup variables to put into the makefile

        # For the PePr call from projects/project/pepr_peaks/
        var_input1 = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', inputGroupA$fullHumanID), sep='', collapse=',')
        var_input2 = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', inputGroupB$fullHumanID), sep='', collapse=',')
        var_chip1 = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', groupA$fullHumanID), sep='', collapse=',')
        var_chip2 = paste(sprintf('$(DIR_PULL_BOWTIE2)/%s_trimmed.fq.gz_aligned.bam', groupB$fullHumanID), sep='', collapse=',')

        # For the prerequisites in the make rule
        var_merged_input1_pre = paste(sprintf('$(DIR_PULL_COVERAGES)/%s_coverage_merged.bdg', inputGroupA$fullHumanID), sep='', collapse=' ')
        var_merged_input2_pre = paste(sprintf('$(DIR_PULL_COVERAGES)/%s_coverage_merged.bdg', inputGroupB$fullHumanID), sep='', collapse=' ')
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
        bigBed_bed = sprintf('$(DIR_PULL_PEPR)/%s_PePr_for_bigBed.bed', var_name)
        annotatr_rdata = sprintf('$(DIR_RDATA)/%s_PePr_annotatr_analysis.RData', var_name)
        bigbed = sprintf('$(DIR_TRACK)/%s_PePr_peaks.bb', var_name)
        simple_bed = sprintf('$(DIR_CLASS_SIMPLE)/%s_PePr_simple_classification.bed', var_name)
        simple_bb = sprintf('$(DIR_TRACK)/%s_PePr_simple_classification.bb', var_name)
        simple_rdata = sprintf('$(DIR_RDATA)/%s_PePr_simple_classification_annotatr_analysis.RData', var_name)

        ########################################################################
        # Variables for the makefile

        # Write the pulldown_compare variables for this comparison
        make_var_pull_compare = c(
            '########################################',
            sprintf('# Workflow for pulldown_compare_%s', i),
            '',
            sprintf('PULLDOWN_COMPARE_%s_PREREQS := %s %s %s %s', i, chip1_bed, bigbed, input_signal, annotatr_rdata),
            sprintf('PULLDOWN_COMPARE_SIMPLE_%s_PREREQS := %s %s %s', i, simple_bed, simple_bb, simple_rdata),
            sprintf('PULLDOWN_COMPARE_%s_INPUT1 := %s', i, var_input1),
            sprintf('PULLDOWN_COMPARE_%s_INPUT2 := %s', i, var_input2),
            sprintf('PULLDOWN_COMPARE_%s_CHIP1 := %s', i, var_chip1),
            sprintf('PULLDOWN_COMPARE_%s_CHIP2 := %s', i, var_chip2),
            sprintf('PULLDOWN_COMPARE_%s_NAME := %s', i, var_name),
            sprintf('PULLDOWN_COMPARE_%s_CLEAN_TMP := %s', i, bigBed_bed))

        # Write the pulldown_compare rule for this comparison
        make_rule_pull_compare = c(
            '',
            '########################################',
            sprintf('.PHONY : pulldown_compare_%s', i),
            sprintf('pulldown_compare_%s : $(PULLDOWN_COMPARE_%s_PREREQS)', i, i),
            '',
            '# Rule for PePr peaks',
            sprintf('%s : %s %s %s %s', chip1_bed, var_input1_pre, var_input2_pre, var_chip1_pre, var_chip2_pre),
            sprintf('    $(PATH_TO_PEPR) --input1=$(PULLDOWN_COMPARE_%s_INPUT1) --input2=$(PULLDOWN_COMPARE_%s_INPUT2) --chip1=$(PULLDOWN_COMPARE_%s_CHIP1) --chip2=$(PULLDOWN_COMPARE_%s_CHIP2) --name=$(PULLDOWN_COMPARE_%s_NAME) --output-directory=$(DIR_PULL_PEPR) $(OPTS_PEPR_%s)', i, i, i, i, i, var_name),
            sprintf('%s : %s', chip2_bed, chip1_bed),
            '',
            '# Rule to combine PePr peaks for bigBed',
            '# NOTE: This script ensures chip1 and chip2 peaks do not overlap',
            '# and then combines the peaks and keeps track of their source',
            sprintf('.INTERMEDIATE : %s', bigBed_bed),
            sprintf('%s : %s %s', bigBed_bed, chip1_bed, chip2_bed),
            sprintf('    bash ../../scripts/pepr_combine.sh $(word 1,$^) $(word 2,$^) $@ $(CHIP1_NAME_%s) $(CHIP2_NAME_%s)', i, i),
            '',
            '# Rule to combine PePr peaks (to save and for annotatr)',
            '# NOTE: Using fold change ($7) and p-value ($8)',
            sprintf('%s : %s %s', combined_bed, chip1_bed, chip2_bed),
            sprintf('    cat <(awk -v OFS="\\t" -v CHIP1=$(CHIP1_NAME_%s) \'{print $$1, $$2, $$3, CHIP1, $$7, "*", $$8}\' $(word 1,$^)) <(awk -v OFS="\\t" -v CHIP2=$(CHIP2_NAME_%s) \'{print $$1, $$2, $$3, CHIP2, $$7, "*", $$8}\' $(word 2,$^)) > $@', i, i),
            '',
            '# Rule for annotatr of PePr peaks',
            sprintf('%s : %s', annotatr_rdata, combined_bed),
            sprintf('    $(PATH_TO_R) ../../scripts/annotatr_annotations.R --file $< --genome $(GENOME) --annot_type PePr --group1 $(CHIP1_NAME_%s) --group0 $(CHIP2_NAME_%s)', i, i),
            '',
            '# Rule to merge input signals from the two groups',
            sprintf('%s : %s %s', input_signal, var_merged_input1_pre, var_merged_input2_pre),
            '    cat $^ | sort -T $(DIR_TMP) -k1,1 -k2,2n | bedtools merge -d 20 | sort -T $(DIR_TMP) -k1,1 -k2,2n > $@',
            '',
            '# Rule for UCSC bigBed track of PePr peaks',
            sprintf('%s : %s', bigbed, bigBed_bed),
            '    $(PATH_TO_BDG2BB) $^ $(CHROM_PATH) $@',
            '',
            '########################################',
            sprintf('.PHONY : pulldown_compare_simple_classification_%s', i),
            sprintf('pulldown_compare_simple_classification_%s : $(PULLDOWN_COMPARE_SIMPLE_%s_PREREQS)', i, i),
            '',
            '# Rule for simple classification of combined PePr peaks',
            sprintf('%s : %s', simple_bed, combined_bed),
            sprintf('    $(PATH_TO_R) ../../scripts/classify_simple.R --project $(PROJECT) --inFile $< --outFile $@ --group1 $(CHIP1_NAME_%s) --group0 $(CHIP2_NAME_%s)', i, i),
            '',
            '# Rule for annotatr of PePr simple classifications',
            sprintf('%s : %s', simple_rdata, simple_bed),
            sprintf('    $(PATH_TO_R) ../../scripts/annotatr_annotations.R --file $< --genome $(GENOME) --annot_type simple_pulldown_PePr --group1 $(CHIP1_NAME_%s) --group0 $(CHIP2_NAME_%s)', i, i),
            '',
            '# Rule for UCSC bigBed track',
            sprintf('%s : %s', simple_bb, simple_bed),
            '    $(PATH_TO_BDG2BB) $< $(CHROM_PATH) $@',
            '',
            '########################################',
            sprintf('# Rule to delete all temporary files from make pulldown_compare_%s',i),
            sprintf('.PHONY : clean_pulldown_compare_tmp_%s', i),
            sprintf('clean_pulldown_compare_tmp_%s :', i),
            sprintf('    rm -f $(PULLDOWN_COMPARE_%s_CLEAN_TMP)', i))

        # Track all the rules for the pulldown compares
        pulldown_compare_rules = c(
            pulldown_compare_rules,
            make_var_pull_compare,
            make_rule_pull_compare)

        ########################################################################
        # OPTS for config.mk
        config_pull_compare = sprintf('########################################
# pulldown_compare_%s configuration options

# Informative names for chip1 and chip2 groups
# CHIP1_NAME should be for the group with higher group number in the project annotation
# file and CHIP2_NAME should be for the group with the lower group number
# If unsure, check the "workflow for pulldown_compare_%s" section of the makefile
CHIP1_NAME_%s := chip1
CHIP2_NAME_%s := chip2

# For PePr parameters see https://ones.ccmb.med.umich.edu/wiki/PePr/
OPTS_PEPR_%s = --file-format=bam --peaktype=sharp --diff --threshold=1e-05 --num-processors=1
',
            i, i, i, i, var_name)

        # Track all the configs for the pulldown compares
        pulldown_compare_configs = c(
            pulldown_compare_configs,
            config_pull_compare)

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
        '################################################################################',
        '# pulldown_compare configuration options',
        '',
        pulldown_compare_configs
    )
    cat(config_master_pull_compare, file = file_config, sep='\n', append=T)

    #######################################
    # PBS script
    # pulldown_compare_q = c(
    #     '#!/bin/bash',
    #     '#### Begin PBS preamble',
    #     '#PBS -N pull_compare',
    #     '#PBS -l nodes=1:ppn=8,pmem=16gb,walltime=24:00:00',
    #     '#PBS -A sartor_lab',
    #     '#PBS -q first',
    #     '#PBS -M rcavalca@umich.edu',
    #     '#PBS -m abe',
    #     '#PBS -j oe',
    #     '#PBS -V',
    #     '#### End PBS preamble',
    #     '# Put your job commands after this line',
    #     sprintf('cd ~/latte/mint/projects/%s/',project),
    #     'make pulldown_compare')
    # cat(pulldown_compare_q, file=sprintf('projects/%s/pbs_jobs/pulldown_compare.q', project), sep='\n')

}
