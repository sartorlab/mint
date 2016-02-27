# Working directory is mint/

library(optparse)

# Parse arguments
option_list = list(
  make_option('--project', type='character'),
  make_option('--genome', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

project = opt$project
genome = opt$genome

########################################################
# Read annotations

annots = read.table(file=sprintf('projects/%s/data/%s_annotation.txt', project, project), header=T, sep='\t', stringsAsFactors=F)

# Split by sample and comparisons
sample_annots = subset(annots, !grepl('comparison', annots$sampleID))
comparison_annots = subset(annots, grepl('comparison', annots$sampleID))

########################################################
# Deal with samples

# Create paths to the sampleID.fastq.gz in projects/project/data/raw_fastqs/
sample_annots$samplePath = apply(sample_annots, 1, function(r){
	sprintf("%s/projects/%s/data/raw_fastqs/%s.fastq.gz", getwd(), project, r['sampleID'])
})

# Create paths and file names to humanID_mark_platform.fastq.gz in
# projects/project/bisulfite/raw_fastqs and projects/project/pulldown/raw_fastqs
# as determined by the platform
humanPath = c()
humanFile = c()
for(i in 1:nrow(sample_annots)) {
	humanID = sample_annots[i,'humanID']
	pull = sample_annots[i,'pulldown']
	bis = sample_annots[i,'bisulfite']
	mc = sample_annots[i,'mc']
	hmc = sample_annots[i,'hmc']
	input = sample_annots[i,'input']

	if( pull == 1 ) {
	  platform = "pulldown"
	} else if ( bis == 1 ) {
	  platform = "bisulfite"
	}

	if( mc == 1 && hmc == 1 ) {
	  mark = "mc_hmc"
	} else if ( mc == 1 && hmc == 0 ) {
	  mark = "mc"
	} else if ( mc == 0 && hmc == 1 ) {
	  mark = "hmc"
	}

	if ( input == 1 ) {
	  humanPath = c(humanPath, sprintf("projects/%s/%s/raw_fastqs/%s_%s_input_%s.fastq.gz", project, platform, humanID, mark, platform))
	  humanFile = c(humanFile, sprintf("%s_%s_input_%s.fastq.gz", humanID, mark, platform))
	} else {
	  humanPath = c(humanPath, sprintf("projects/%s/%s/raw_fastqs/%s_%s_%s.fastq.gz", project, platform, humanID, mark, platform))
	  humanFile = c(humanFile, sprintf("%s_%s_%s.fastq.gz", humanID, mark, platform))
	}
}

# Add the humanPath and humanFile columns to sample_annots
sample_annots$humanPath = humanPath
sample_annots$humanFile = humanFile

# Create the variables for the project specific variables.mk that are
# included in the makefile
pullAlignFiles = subset(sample_annots, pulldown == 1)$humanFile
bisAlignFiles = subset(sample_annots, bisulfite == 1)$humanFile
pullSampleFiles = subset(sample_annots, pulldown == 1 & input == 0)$humanFile

# Write the variables to projects/project/variables.mk
message(sprintf('Writing projects/%s/variables.mk ...', project))
vars = c(
	sprintf('BISULFITE_FASTQ_FILES := %s', paste(bisAlignFiles, collapse=' ')),
	sprintf('PULLDOWN_FASTQ_FILES := %s', paste(pullAlignFiles, collapse=' ')),
	sprintf('PULLDOWN_SAMPLE_FILES := %s', paste(pullSampleFiles, collapse=' '))
	)
cat(vars, file = sprintf('projects/%s/variables.mk', project), sep='\n')

# Create the symlinks between sampleID.fastq.gz and humanID_mark_platform.fastq.gz
# All downstream analysis is done with respect to the humanID_mark_platform.fastq.gz
# files so that context is readable from the file name
for(i in 1:nrow(sample_annots)) {
	# symlink data/raw_fastq to pulldown/raw_fastq or bisulfite/raw_fastq
	command = sprintf('ln -s %s %s', sample_annots[i,'samplePath'], sample_annots[i,'humanPath'])
	message(sprintf('Creating symlink: %s', command))
	system(command)
}

########################################################
# Deal with comparisons if there are any

if(nrow(comparison_annots) > 0) {
	# Add parameter columns to comparison_annots for PePr and methylSig
	comparison_annots$pepr_input1 = NA
	comparison_annots$pepr_input2 = NA
	comparison_annots$pepr_chip1 = NA
	comparison_annots$pepr_chip2 = NA
	comparison_annots$pepr_name = NA
	comparison_annots$msig_project = NA
	comparison_annots$msig_cytfiles = NA
	comparison_annots$msig_sampleids = NA
	comparison_annots$msig_treatment = NA
	comparison_annots$msig_assembly = NA
	comparison_annots$msig_comparison = NA

	pulldown_compares = c()
	bisulfite_compares = c()

	for(i in 1:nrow(comparison_annots)) {
		# The order of the groups is important for biological interpretation

		# Establish row variables
		projectID = comparison_annots[i,'projectID']
		sampleID = comparison_annots[i,'sampleID']
		humanID = comparison_annots[i,'humanID']
		pull = comparison_annots[i,'pulldown']
		bis = comparison_annots[i,'bisulfite']
		mc = comparison_annots[i,'mc']
		hmc = comparison_annots[i,'hmc']
		input = comparison_annots[i,'input']
		group = comparison_annots[i,'group']

		if( pull == 1 ) {
		  platform = "pulldown"
		} else if ( bis == 1 ) {
		  platform = "bisulfite"
		}

		if( mc == 1 && hmc == 1 ) {
		  mark = "mc_hmc"
		} else if ( mc == 1 && hmc == 0 ) {
		  mark = "mc"
		} else if ( mc == 0 && hmc == 1 ) {
		  mark = "hmc"
		}

		# Sorting this ensures that the lower group number is A
		groups = sort(as.integer(unlist(strsplit(group, ','))))

		groupA = subset(sample_annots,
			grepl(groups[1], sample_annots$group) &
			sample_annots$pulldown == pull &
			sample_annots$bisulfite == bis &
			sample_annots$mc == mc &
			sample_annots$hmc == hmc &
			sample_annots$input == 0)
		groupB = subset(sample_annots,
			grepl(groups[2], sample_annots$group) &
			sample_annots$pulldown == pull &
			sample_annots$bisulfite == bis &
			sample_annots$mc == mc &
			sample_annots$hmc == hmc &
			sample_annots$input == 0)

		inputGroupA = subset(sample_annots,
			grepl(groups[1], sample_annots$group) &
			sample_annots$pulldown == pull &
			sample_annots$bisulfite == bis &
			sample_annots$mc == mc &
			sample_annots$hmc == hmc &
			sample_annots$input == 1)
		inputGroupB = subset(sample_annots,
			grepl(groups[2], sample_annots$group) &
			sample_annots$pulldown == pull &
			sample_annots$bisulfite == bis &
			sample_annots$mc == mc &
			sample_annots$hmc == hmc &
			sample_annots$input == 1)

		if(platform == 'pulldown') {

			# For the PePr call from projects/project/pepr_peaks/
			input1 = paste('../bowtie2_bams/', gsub('.fastq.gz', '_trimmed.fq.gz_aligned.bam', inputGroupA$humanFile), sep='', collapse=',')
			input2 = paste('../bowtie2_bams/', gsub('.fastq.gz', '_trimmed.fq.gz_aligned.bam', inputGroupB$humanFile), sep='', collapse=',')
			chip1 = paste('../bowtie2_bams/', gsub('.fastq.gz', '_trimmed.fq.gz_aligned.bam', groupA$humanFile), sep='', collapse=',')
			chip2 = paste('../bowtie2_bams/', gsub('.fastq.gz', '_trimmed.fq.gz_aligned.bam', groupB$humanFile), sep='', collapse=',')

			# For the prerequisites in the make rule
			input1_pre = paste('$(DIR_PULL_BOWTIE2)/', gsub('.fastq.gz', '_trimmed.fq.gz_aligned.bam', inputGroupA$humanFile), sep='', collapse=' ')
			input2_pre = paste('$(DIR_PULL_BOWTIE2)/', gsub('.fastq.gz', '_trimmed.fq.gz_aligned.bam', inputGroupB$humanFile), sep='', collapse=' ')
			chip1_pre = paste('$(DIR_PULL_BOWTIE2)/', gsub('.fastq.gz', '_trimmed.fq.gz_aligned.bam', groupA$humanFile), sep='', collapse=' ')
			chip2_pre = paste('$(DIR_PULL_BOWTIE2)/', gsub('.fastq.gz', '_trimmed.fq.gz_aligned.bam', groupB$humanFile), sep='', collapse=' ')
			name = sprintf('%s_%s_%s', humanID, mark, platform)

			up_bed = sprintf('$(DIR_PULL_PEPR)/%s__PePr_up_peaks.bed', name)
			down_bed = sprintf('$(DIR_PULL_PEPR)/%s__PePr_down_peaks.bed', name)
			combined_bed = sprintf('$(DIR_PULL_PEPR)/%s_PePr_combined.bed', name)
			bigbed = sprintf('$(DIR_TRACK)/%s_PePr_peaks.bb', name)

			# Write the PePr variables for this comparison
			pepr_vars = c(
				sprintf('PULLDOWN_COMPARE_%s_PREREQS := %s %s', i, up_bed, bigbed),
				sprintf('PULLDOWN_COMPARE_%s_INPUT1 := %s', i, input1),
				sprintf('PULLDOWN_COMPARE_%s_INPUT2 := %s', i, input2),
				sprintf('PULLDOWN_COMPARE_%s_CHIP1 := %s', i, chip1),
				sprintf('PULLDOWN_COMPARE_%s_CHIP2 := %s', i, chip2),
				sprintf('PULLDOWN_COMPARE_%s_NAME := %s', i, name))
			cat(pepr_vars, file = sprintf('projects/%s/variables.mk', project), sep='\n', append=T)

			# Write the PePr rule for this comparison
			pepr_rule = c(
				sprintf('.PHONY : pulldown_compare_%s', i),
				sprintf('pulldown_compare_%s : $(PULLDOWN_COMPARE_%s_PREREQS)', i, i),
				'',
				sprintf('%s : %s %s %s %s', up_bed, input1_pre, input2_pre, chip1_pre, chip2_pre),
				'	cd projects/$(PROJECT)/pepr_peaks; \\',
				sprintf('	$(PATH_TO_PEPR) --input1=$(PULLDOWN_COMPARE_%s_INPUT1) --input2=$(PULLDOWN_COMPARE_%s_INPUT2) --chip1=$(PULLDOWN_COMPARE_%s_CHIP1) --chip2=$(PULLDOWN_COMPARE_%s_CHIP2) --name=$(PULLDOWN_COMPARE_%s_NAME) $(OPTS_PEPR)', i, i, i, i, i),
				sprintf('%s : %s', down_bed, up_bed),
				'',
				sprintf('%s : %s %s', combined_bed, up_bed, down_bed),
				'	bash ../../scripts/combine_pepr.sh $(word 1,$^) $(word 2,$^) $@',
				'',
				sprintf('%s : %s', bigbed, combined_bed),
				'	bedToBigBed $^ $(CHROM_PATH) $@',
				'')
			cat(pepr_rule, file = sprintf('projects/%s/makefile', project), sep='\n', append=T)

			# Track the number of pulldown compares
			pulldown_compares = c(pulldown_compares, sprintf('pulldown_compare_%s', i))

		} else if (platform == 'bisulfite') {

			cytfiles = paste(c(
				paste('$(DIR_BIS_BISMARK)/', gsub('.fastq.gz', '_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt', groupA$humanFile), sep=''),
				paste('$(DIR_BIS_BISMARK)/', gsub('.fastq.gz', '_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt', groupB$humanFile), sep='')), collapse=',')
			sampleids = paste(c(
				gsub('.fastq.gz','',groupA$humanFile),
				gsub('.fastq.gz','',groupB$humanFile)), collapse=',')
			treatment = paste(c(
				rep.int(groups[1],nrow(groupA)),
				rep.int(groups[2],nrow(groupB))), collapse=',')
			assembly = genome
			comparison = sprintf('%s_%s_%s', humanID, mark, platform)

			# make prerequisites for methylSig
			cytfiles_pre = paste(c(
				paste('$(DIR_BIS_BISMARK)/', gsub('.fastq.gz', '_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt', groupA$humanFile), sep=''),
				paste('$(DIR_BIS_BISMARK)/', gsub('.fastq.gz', '_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt', groupB$humanFile), sep='')), collapse=' ')

			# make targets
			msig_results = sprintf('$(DIR_BIS_MSIG)/%s_methylSig.txt', comparison)
			msig_tmp_results = sprintf('$(DIR_BIS_MSIG)/%s_methylSig_tmp.txt', comparison)
			msig_bigwig = sprintf('$(DIR_TRACK)/%s_methylSig.bw', comparison)

			# Write the methylSig variables for this comparison
			msig_vars = c(
				sprintf('BISULFITE_COMPARE_%s_PREREQS := %s %s', i, msig_results, msig_bigwig),
				sprintf('BISULFITE_COMPARE_%s_CYTFILES := %s', i, cytfiles),
				sprintf('BISULFITE_COMPARE_%s_SAMPLEIDS := %s', i, sampleids),
				sprintf('BISULFITE_COMPARE_%s_TREATMENT := %s', i, treatment),
				sprintf('BISULFITE_COMPARE_%s_COMPARISON := %s', i, comparison))
			cat(msig_vars, file = sprintf('projects/%s/variables.mk', project), sep='\n', append=T)

			# Write the methylSig rule for this comparison
			msig_rule = c(
				sprintf('.PHONY : bisulfite_compare_%s', i),
				sprintf('bisulfite_compare_%s : $(BISULFITE_COMPARE_%s_PREREQS)', i, i),
				'',
				sprintf('%s : %s', msig_results, cytfiles_pre),
				sprintf('	Rscript ../../scripts/process_bisulfite_comparison_run_methylSig.R --project $(PROJECT) --cytfiles $(BISULFITE_COMPARE_%s_CYTFILES) --sampleids $(BISULFITE_COMPARE_%s_SAMPLEIDS) --treatment $(BISULFITE_COMPARE_%s_TREATMENT) --assembly $(GENOME) --pipeline mint --comparison $(BISULFITE_COMPARE_%s_COMPARISON) $(OPTS_METHYLSIG)', i, i, i, i),
				'',
				sprintf('.INTERMEDIATE : %s', msig_tmp_results),
				sprintf('%s : %s', msig_tmp_results, msig_results),
				"	awk -v OFS='\\t' '$$5 < 0.05 {print $$1, $$2, $$3, $$7 }' $^ | sort -T . -k1,1 -k2,2n > $@",
				'',
				sprintf('%s : %s', msig_bigwig, msig_tmp_results),
				'	bedGraphToBigWig $^ $(CHROM_PATH) $@',
				'')
			cat(msig_rule, file = sprintf('projects/%s/makefile', project), sep='\n', append=T)

			# Track the number of bisulfite compares
			bisulfite_compares = c(bisulfite_compares, sprintf('bisulfite_compare_%s', i))
		}
	}

	# Make master rules for bisulfite_compare and pulldown_compare which
	# have the individual compares as their prerequisites
	if(length(bisulfite_compares) > 0) {
		master_bisulfite_compare_rule = c(
			'.PHONY : bisulfite_compare',
			sprintf('bisulfite_compare : bisulfite_align %s', paste(bisulfite_compares, collapse=' ')),
			'')
		cat(master_bisulfite_compare_rule, file = sprintf('projects/%s/makefile', project), sep='\n', append=T)
	}

	if(length(pulldown_compares) > 0) {
		master_pulldown_compare_rule = c(
			'.PHONY : pulldown_compare',
			sprintf('pulldown_compare : pulldown_align %s', paste(pulldown_compares, collapse=' ')),
			'')
		cat(master_pulldown_compare_rule, file = sprintf('projects/%s/makefile', project), sep='\n', append=T)
	}
}

########################################################
# Build PBS scripts

bisulfite_align_q = c(
	'#!/bin/bash',
	'#### Begin PBS preamble',
	'#PBS -N bis_align',
	'#PBS -l procs=10,mem=80gb,walltime=6:00:00',
	'#PBS -A sartor_lab',
	'#PBS -q first',
	'#PBS -M rcavalca@umich.edu',
	'#PBS -m abe',
	'#PBS -j oe',
	'#PBS -V',
	'#### End PBS preamble',
	'# Put your job commands after this line',
	sprintf('cd ~/latte/mint/projects/%s/',project),
	'make -j 2 bisulfite_align')
cat(bisulfite_align_q, file=sprintf('projects/%s/bisulfite_align.q', project), sep='\n')

pulldown_align_q = c(
	'#!/bin/bash',
	'#### Begin PBS preamble',
	'#PBS -N pull_align',
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
	'make -j 4 pulldown_align')
cat(pulldown_align_q, file=sprintf('projects/%s/pulldown_align.q', project), sep='\n')

pulldown_sample_q = c(
	'#!/bin/bash',
	'#### Begin PBS preamble',
	'#PBS -N pull_sample',
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
	'make -j 4 pulldown_sample')
cat(pulldown_sample_q, file=sprintf('projects/%s/pulldown_sample.q', project), sep='\n')

pulldown_compare_q = c(
	'#!/bin/bash',
	'#### Begin PBS preamble',
	'#PBS -N pull_compare',
	'#PBS -l procs=1,mem=32gb,walltime=6:00:00',
	'#PBS -A sartor_lab',
	'#PBS -q first',
	'#PBS -M rcavalca@umich.edu',
	'#PBS -m abe',
	'#PBS -j oe',
	'#PBS -V',
	'#### End PBS preamble',
	'# Put your job commands after this line',
	sprintf('cd ~/latte/mint/projects/%s/',project),
	'make -j pulldown_compare')
cat(pulldown_compare_q, file=sprintf('projects/%s/pulldown_compare.q', project), sep='\n')

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
	'make -j bisulfite_compare')
cat(bisulfite_compare_q, file=sprintf('projects/%s/bisulfite_compare.q', project), sep='\n')
