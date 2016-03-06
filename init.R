# Working directory is mint/

################################################################################
# Parse arguments
library(optparse)

option_list = list(
  make_option('--project', type='character'),
  make_option('--genome', type='character'),
  make_option('--datapath', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

project = opt$project
genome = opt$genome
datapath = opt$datapath

file_make = sprintf('projects/%s/makefile', project)
dir_hub = sprintf('projects/%s/%s_hub', project, project)
dir_track = sprintf('%s/%s', dir_hub, genome)

################################################################################
################################################################################
################################################################################

################################################################################
# READ project information
annots = read.table(file = sprintf('projects/%s_annotation.txt', project), header = T, sep = '\t', stringsAsFactors = F)

# Add the fullHumanID column which will be the prefix for all created files
annots$fullHumanID = NA
for(i in 1:nrow(annots)) {
	humanID = annots[i,'humanID']
	pull = annots[i,'pulldown']
	bis = annots[i,'bisulfite']
	mc = annots[i,'mc']
	hmc = annots[i,'hmc']
	input = annots[i,'input']

	if( pull == 1 ) {
	  platform = "pulldown"
	} else if ( bis == 1 ) {
	  platform = "bisulfite"
	} else {
		stop('Annotation Error: For each row, either the pulldown or the bisulfite column must be 1.')
	}

	if( mc == 1 && hmc == 1 ) {
	  mark = "mc_hmc"
	} else if ( mc == 1 && hmc == 0 ) {
	  mark = "mc"
	} else if ( mc == 0 && hmc == 1 ) {
	  mark = "hmc"
	} else {
		stop('Annotation Error: For each row, mc and/or hmc must be 1.')
	}

	if ( input == 1 ) {
		annots[i,'fullHumanID'] = sprintf("%s_%s_input_%s", humanID, mark, platform)
	} else {
		annots[i,'fullHumanID'] = sprintf("%s_%s_%s", humanID, mark, platform)
	}
}

# Split samples by bisulfite and pulldown
bisulfite = subset(annots, bisulfite == 1)
pulldown = subset(annots, pulldown == 1)

# Determine the sample
samples = subset(annots, !grepl('comparison', sampleID))
comparisons = subset(annots, grepl('comparison', sampleID))

# Split by sample and comparisons
bisulfite_samples = subset(bisulfite, !grepl('comparison', sampleID))
bisulfite_comparisons = subset(bisulfite, grepl('comparison', sampleID))

pulldown_samples = subset(pulldown, !grepl('comparison', sampleID))
pulldown_samples_noinput = subset(pulldown, !grepl('comparison', sampleID) & input == 0)
pulldown_comparisons = subset(pulldown, grepl('comparison', sampleID))

# Control variables
bool_bis_samp = nrow(bisulfite_samples) > 0
bool_pull_samp = nrow(pulldown_samples) > 0
bool_bis_comp = nrow(bisulfite_comparisons) > 0
bool_pull_comp = nrow(pulldown_comparisons) > 0

# XXX: Error for pure bisulfite setup. Will eventually support this.
if (bool_bis_samp && !bool_pull_samp) {
	stop('Error: Pure bisulfite experimental setups are not currently supported.')
}

# NOTE: ADD ERROR CHECKING
# 1. Are there input samples for pulldowns?
# 2. If there are comparisons, do the group numbers correspond to the groups assigned to the samples?

# NOTE: Have to think about how to handle inputs shared between mc and hmc pulldowns.
# One option is to create an mc and hmc symlink to the same input file, but then all
# analysis would be doubled for that input file...

################################################################################
################################################################################
################################################################################

################################################################################
# INITIALIZATION: Directory creation

# Always create these folders
setup_commands = c(
	sprintf('mkdir projects/%s', project),
	sprintf('mkdir projects/%s/data', project),
	sprintf('mkdir projects/%s/data/raw_fastqs', project),
	sprintf('mkdir projects/%s/summary', project),
	sprintf('mkdir projects/%s/summary/{figures,tables,reports}', project),
	sprintf('mkdir projects/%s/%s_hub', project, project),
	sprintf('mkdir projects/%s/%s_hub/%s', project, project, genome),
	sprintf('mkdir projects/%s/classifications', project),
	sprintf('mkdir projects/%s/classifications/{simple,sample}', project)
)

# Create folders for bisulfite samples if there are any
if(bool_bis_samp) {
	setup_commands = c(
		setup_commands,
		sprintf('mkdir projects/%s/bisulfite', project),
		sprintf('mkdir projects/%s/bisulfite/{raw_fastqs,raw_fastqcs,trim_fastqs,trim_fastqcs,bismark}', project)
	)
}
# Create folders for bisulfite comparisons if there are any
if(bool_bis_comp) {
	setup_commands = c(
		setup_commands,
		sprintf('mkdir projects/%s/bisulfite/methylsig_calls', project)
	)
}
# Create folders for pulldown samples if there are ny
if(bool_pull_samp) {
	setup_commands = c(
		setup_commands,
		sprintf('mkdir projects/%s/pulldown', project),
		sprintf('mkdir projects/%s/pulldown/{raw_fastqs,raw_fastqcs,trim_fastqs,trim_fastqcs,bowtie2_bams,pulldown_coverages,macs2_peaks}', project)
	)
}
# Create folders for pulldown comparisons if there are any
if(bool_pull_comp) {
	setup_commands = c(
		setup_commands,
		sprintf('mkdir projects/%s/pulldown/pepr_peaks', project)
	)
}
# Create folders for comparison classification if any comparisons are done
if(bool_bis_comp || bool_pull_comp) {
	setup_commands = c(
		setup_commands,
		sprintf('mkdir projects/%s/classifications/comparison', project)
	)
}

################################################################################
# INITIALIZATION: File copying into project folder

setup_commands = c(
	setup_commands,
	sprintf('cp template_makefile projects/%s/makefile', project),
	sprintf('cp template_config.mk projects/%s/config.mk', project),
	sprintf('cp narrowPeak.as projects/%s/', project),
	sprintf('cp projects/%s_annotation.txt projects/%s/data', project, project)
)

################################################################################
# INITIALIZATION: Run the setup_commands
for(command in setup_commands) {
	message(command)
	system(command)
}

################################################################################
# INITIALIZATION: Symlink files in datapath into projects/project/data/raw_fastqs

datafiles = list.files(datapath, full.names = TRUE)
for(file in datafiles) {
	command = sprintf('ln -s %s projects/%s/data/raw_fastqs/%s', file, project, basename(file))
	message(command)
	system(command)
}

################################################################################
# INITIALIZATION: Symlink files in projects/project/data/raw_fastqs to the
# structured file names in raw_fastqs folders in pulldown and bisulfite as needed

if(bool_bis_samp) {
	for(i in 1:nrow(bisulfite_samples)) {
		command = sprintf('ln -s %s/projects/%s/data/raw_fastqs/%s.fastq.gz %s/projects/%s/bisulfite/raw_fastqs/%s.fastq.gz',
			getwd(), project, bisulfite_samples[i,'sampleID'], getwd(), project, bisulfite_samples[i,'fullHumanID'])
		message(command)
		system(command)
	}
}

if(bool_pull_samp) {
	for(i in 1:nrow(pulldown_samples)) {
		command = sprintf('ln -s %s/projects/%s/data/raw_fastqs/%s.fastq.gz %s/projects/%s/pulldown/raw_fastqs/%s.fastq.gz',
			getwd(), project, pulldown_samples[i,'sampleID'], getwd(), project, pulldown_samples[i,'fullHumanID'])
		message(command)
		system(command)
	}
}

################################################################################
################################################################################
################################################################################

################################################################################
# Track hub files
# project_hub/genomes.txt
# project_hub/hub.txt
# project_hub/project_name_hub.html

# The trackDb.txt file is used no matter the situation
# Make sure we start fresh each time
hubtrackdbfile = sprintf('%s/trackDb.txt', dir_track, genome)
cat('', file=hubtrackdbfile)

# genomes.txt
  genomes = c(
    sprintf('genome %s', genome),
    sprintf('trackDb %s/trackDb.txt', genome))
  cat(genomes, file=sprintf('%s/genomes.txt', dir_hub), sep='\n')
# hub.txt
  hub = c(
    sprintf('hub %s_hub', project),
    sprintf('shortLabel %s', project),
    sprintf('longLabel %s', project),
    'genomesFile genomes.txt',
    'email rcavalca@umich.edu',
    sprintf('descriptionUrl %s_hub.html', project))
  cat(hub, file=sprintf('%s/hub.txt', dir_hub), sep='\n')
# project_name_hub.html
# TODO: Document what's in the hub!
  cat(' ', file=sprintf('%s/%s_hub.html',dir_hub, project))

priority = 1
# superTracks
if(bool_bis_comp || bool_pull_comp) {
	for(comparison in unique(comparisons$humanID)) {
		trackEntry = c(
			sprintf('track %s_group_comparison', comparison),
			'superTrack on show',
			sprintf('group %s_group_comparison', comparison),
			sprintf('shortLabel %s comparison', comparison),
			sprintf('longLabel %s comparison', comparison),
			sprintf('priority %s', priority),
			' ')
		cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
		priority = priority + 1
	}
}

if(bool_bis_samp || bool_pull_samp) {
	for(sample in sort(unique(samples$humanID))) {
		trackEntry = c(
			sprintf('track %s_sample', sample),
			'superTrack on show',
			sprintf('group %s_sample', sample),
			sprintf('shortLabel %s_sample', sample),
			sprintf('longLabel Sample-level tracks for %s', sample),
			sprintf('priority %s', priority),
			' ')
		cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
		priority = priority + 1
	}
}

################################################################################
################################################################################
################################################################################

################################################################################
# MAKEFILE: bisulfite_align rules

if(bool_bis_samp) {

make_var_bis_align_prefix = sprintf('
################################################################################
# Workflow for bisulfite_align
BISULFITE_ALIGN_PREFIXES := %s', paste(bisulfite_samples$fullHumanID, collapse=' '))

make_var_bis_align = 'BISULFITE_ALIGN_PREREQS := 	$(patsubst %,$(DIR_TRACK)/%_simple_classification.bb,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_CLASS_SIMPLE)/%_simple_classification.bed,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_TRACK)/%_trimmed.fq.gz_bismark_bt2.bw,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_annotatr.txt,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph.gz,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bam,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_TRIM_FASTQCS)/%_trimmed.fq_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz,$(BISULFITE_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip,$(BISULFITE_ALIGN_PREFIXES))'

# NOTE: This cannot be indented because they would mess up the makefile
make_rule_bis_align = '
.PHONY : bisulfite_align
bisulfite_align : $(BISULFITE_ALIGN_PREREQS)

# Rule for UCSC bigBed track of simple classifiation
$(DIR_TRACK)/%_bisulfite_simple_classification.bb : $(DIR_CLASS_SIMPLE)/%_bisulfite_simple_classification.bed
	bedToBigBed $< $(CHROM_PATH) $@

# Simple classification for percent methylation
$(DIR_CLASS_SIMPLE)/%_bisulfite_simple_classification.bed : $(DIR_BIS_BISMARK)/%_bisulfite_trimmed.fq.gz_bismark_bt2.bedGraph.gz
	Rscript ../../scripts/classify_simple.R --project $(PROJECT) --inFile $< --outFile $@

# Rules for UCSC bigWig track
$(DIR_TRACK)/%_trimmed.fq.gz_bismark_bt2.bw : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph
	bedGraphToBigWig $< $(CHROM_PATH) $@

.INTERMEDIATE : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph
$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph.gz
	gunzip -c $< | awk \'NR > 1 {print $$0}\' | sort -T . -k1,1 -k2,2n > $@

# Rule for methylSig input
$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_methylSig.txt : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt
	awk -f ../../scripts/to_methylSig.awk $< | sort -T . -k2,2 -k3,3n > $@

# Rule for annotatr input
$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report_for_annotatr.txt : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt
	awk -f ../../scripts/to_annotatr.awk $< | sort -T . -k1,1 -k2,2n > $@

# Rule for bismark methylation extractor
$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bedGraph.gz $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.CpG_report.txt : $(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bam
	cd $(DIR_BIS_BISMARK); \\
	bismark_methylation_extractor $(OPTS_EXTRACTOR) $(<F)

# Rule for bismark alignment
$(DIR_BIS_BISMARK)/%_trimmed.fq.gz_bismark_bt2.bam : $(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz $(DIR_BIS_TRIM_FASTQCS)/%_trimmed.fq_fastqc.zip
	bismark $(OPTS_BISMARK) --output_dir $(@D) --temp_dir $(@D) $<
	samtools sort $@ $(patsubst %.bam,%,$@)
	samtools index $@

# Rule for FastQC on trimmed
$(DIR_BIS_TRIM_FASTQCS)/%_trimmed.fq_fastqc.zip : $(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz
	fastqc $(OPTS_FASTQC) --outdir $(@D) $<

# Rule for trim_galore
$(DIR_BIS_TRIM_FASTQS)/%_trimmed.fq.gz : $(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip
	trim_galore $(OPTS_TRIMGALORE_BISULFITE) --output_dir $(@D) $(DIR_BIS_RAW_FASTQS)/$*.fastq.gz

# Rule for FastQC on raw
$(DIR_BIS_RAW_FASTQCS)/%_fastqc.zip :
	fastqc $(OPTS_FASTQC) --outdir $(@D) $(DIR_BIS_RAW_FASTQS)/$*.fastq.gz
'

cat(make_var_bis_align_prefix, file = file_make, sep = '\n', append = TRUE)
cat(make_var_bis_align, file = file_make, sep = '\n', append = TRUE)
cat(make_rule_bis_align, file = file_make, sep = '\n', append = TRUE)

#######################################
# PBS script
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

for(i in 1:nrow(bisulfite_samples)) {
	# trackDb.txt entry for Bismark methylation calls
	trackEntry = c(
	  sprintf('track %s_pct_meth', bisulfite_samples[i,'fullHumanID']),
	  sprintf('parent %s_sample', bisulfite_samples[i,'humanID']),
	  sprintf('bigDataUrl %s_trimmed.fq.gz_bismark_bt2.bw', bisulfite_samples[i,'fullHumanID']),
	  sprintf('shortLabel %s_pct_meth', bisulfite_samples[i,'fullHumanID']),
	  sprintf('longLabel %s_percent_methylation', bisulfite_samples[i,'fullHumanID']),
	  'visibility full',
	  'viewLimits 0:100',
	  'type bigWig',
	  'priority 1.4',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

	# trackDb.txt entry for bisulfite simple classification results
	trackEntry = c(
	  sprintf('track %s_simple_class', bisulfite_samples[i,'fullHumanID']),
	  sprintf('parent %s_sample', bisulfite_samples[i,'humanID']),
	  sprintf('bigDataUrl %s_simple_classification.bb', bisulfite_samples[i,'fullHumanID']),
	  sprintf('shortLabel %s_simp_class', bisulfite_samples[i,'fullHumanID']),
	  sprintf('longLabel %s_simple_classification', bisulfite_samples[i,'fullHumanID']),
	  'visibility pack',
	  'itemRgb on',
	  'type bigBed 9 .',
	  'priority 1.2',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
}

}

################################################################################

################################################################################
# MAKEFILE: pulldown_align rules

if(bool_pull_samp) {

make_var_pull_align_prefix = sprintf('
################################################################################
# Workflow for pulldown_align
PULLDOWN_ALIGN_PREFIXES := %s', paste(pulldown_samples$fullHumanID, collapse=' '))

make_var_pull_align = 'PULLDOWN_ALIGN_PREREQS :=  $(patsubst %,$(DIR_TRACK)/%_coverage.bw,$(PULLDOWN_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_PULL_COVERAGES)/%_coverage.bdg,$(PULLDOWN_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_PULL_COVERAGES)/%_merged_coverage.bdg,$(PULLDOWN_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_PULL_BOWTIE2)/%_trimmed.fq.gz_aligned.bam,$(PULLDOWN_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_PULL_TRIM_FASTQCS)/%_trimmed.fq_fastqc.zip,$(PULLDOWN_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_PULL_TRIM_FASTQS)/%_trimmed.fq.gz,$(PULLDOWN_ALIGN_PREFIXES)) \\
												$(patsubst %,$(DIR_PULL_RAW_FASTQCS)/%_fastqc.zip,$(PULLDOWN_ALIGN_PREFIXES))'

# NOTE: This cannot be indented because they would mess up the makefile
make_rule_pull_align = '
.PHONY : pulldown_align
pulldown_align : $(PULLDOWN_ALIGN_PREREQS)

# Rule for UCSC bigWig track
$(DIR_TRACK)/%_coverage.bw : $(DIR_PULL_COVERAGES)/%_coverage.bdg
	bedGraphToBigWig $< $(CHROM_PATH) $@

# Rule for coverage bedGraph
$(DIR_PULL_COVERAGES)/%_coverage.bdg : $(DIR_PULL_BOWTIE2)/%_trimmed.fq.gz_aligned.bam
	bedtools genomecov -bg -g $(CHROM_PATH) -ibam $< | sort -T . -k1,1 -k2,2n > $@

# Rule for merged coverage BED
# For use in signal BEDs downstream
$(DIR_PULL_COVERAGES)/%_merged_coverage.bdg : $(DIR_PULL_COVERAGES)/%_coverage.bdg
	bedtools merge -d 20 -i $< | sort -T . -k1,1 -k2,2n > $@

# Rule for bowtie2 alignment
$(DIR_PULL_BOWTIE2)/%_trimmed.fq.gz_aligned.bam : $(DIR_PULL_TRIM_FASTQS)/%_trimmed.fq.gz $(DIR_PULL_TRIM_FASTQCS)/%_trimmed.fq_fastqc.zip
	bowtie2 $(OPTS_BOWTIE2) $< | samtools view -bS - > $@
	samtools sort $@ $(patsubst %.bam,%,$@)
	samtools index $@

# Rule for FastQC on trimmed
$(DIR_PULL_TRIM_FASTQCS)/%_trimmed.fq_fastqc.zip : $(DIR_PULL_TRIM_FASTQS)/%_trimmed.fq.gz
	fastqc $(OPTS_FASTQC) --outdir $(@D) $<

# Rule for trim_galore
$(DIR_PULL_TRIM_FASTQS)/%_trimmed.fq.gz : $(DIR_PULL_RAW_FASTQCS)/%_fastqc.zip
	trim_galore $(OPTS_TRIMGALORE_PULLDOWN) --output_dir $(@D) $(DIR_PULL_RAW_FASTQS)/$*.fastq.gz

# Rule for FastQC on raw
$(DIR_PULL_RAW_FASTQCS)/%_fastqc.zip :
	fastqc $(OPTS_FASTQC) --outdir $(@D) $(DIR_PULL_RAW_FASTQS)/$*.fastq.gz
'

cat(make_var_pull_align_prefix, file = file_make, sep = '\n', append = TRUE)
cat(make_var_pull_align, file = file_make, sep = '\n', append = TRUE)
cat(make_rule_pull_align, file = file_make, sep = '\n', append = TRUE)

#######################################
# PBS script
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

for(i in 1:nrow(pulldown_samples)) {
	# trackDb.txt entry for chip/input pulldown coverages
	trackEntry = c(
	  sprintf('track %s_cov', pulldown_samples[i,'fullHumanID']),
	  sprintf('parent %s_sample', pulldown_samples[i,'humanID']),
	  sprintf('bigDataUrl %s_coverage.bw', pulldown_samples[i,'fullHumanID']),
	  sprintf('shortLabel %s_cov', pulldown_samples[i,'fullHumanID']),
	  sprintf('longLabel %s_coverage', pulldown_samples[i,'fullHumanID']),
	  'visibility full',
	  'autoScale on',
	  'type bigWig',
	  'priority 1.6',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

}

################################################################################

################################################################################
# MAKEFILE: pulldown_sample rules

make_var_pull_samp_prefix = sprintf('
################################################################################
# Workflow for pulldown_sample
PULLDOWN_SAMPLE_PREFIXES := %s', paste(pulldown_samples_noinput$fullHumanID, collapse=' '))


make_var_pull_samp = 'PULLDOWN_SAMPLE_PREREQS :=	$(patsubst %,$(DIR_TRACK)/%_simple_classification.bb,$(PULLDOWN_SAMPLE_PREFIXES)) \\
												$(patsubst %,$(DIR_CLASS_SIMPLE)/%_simple_classification.bed,$(PULLDOWN_SAMPLE_PREFIXES)) \\
												$(patsubst %,$(DIR_TRACK)/%_macs2_peaks.bb,$(PULLDOWN_SAMPLE_PREFIXES)) \\
												$(patsubst %,$(DIR_PULL_MACS)/%_macs2_peaks.narrowPeak,$(PULLDOWN_SAMPLE_PREFIXES))'

# NOTE: This cannot be indented because they would mess up the makefile
make_rule_pull_samp = '
.PHONY : pulldown_sample
pulldown_sample : pulldown_align $(PULLDOWN_SAMPLE_PREREQS)

# Rule for UCSC bigBed track of simple classifiation
$(DIR_TRACK)/%_pulldown_simple_classification.bb : $(DIR_CLASS_SIMPLE)/%_pulldown_simple_classification.bed
	bedToBigBed $< $(CHROM_PATH) $@

# Simple classification
$(DIR_CLASS_SIMPLE)/%_pulldown_simple_classification.bed : $(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak
	Rscript ../../scripts/classify_simple.R --project $(PROJECT) --inFile $< --outFile $@

# Rule for UCSC bigBed track
$(DIR_TRACK)/%_macs2_peaks.bb : $(DIR_PULL_MACS)/%_macs2_peaks_tmp.narrowPeak
	bedToBigBed -type=bed6+4 -as=narrowPeak.as $^ $(CHROM_PATH) $@

# Rule for macs2 peak fix
.INTERMEDIATE : $(DIR_PULL_MACS)/%_macs2_peaks_tmp.narrowPeak
$(DIR_PULL_MACS)/%_macs2_peaks_tmp.narrowPeak : $(DIR_PULL_MACS)/%_macs2_peaks.narrowPeak
	awk -f ../../scripts/fix_narrowPeak.awk $^ | sort -T . -k1,1 -k2,2n > $@

# Rule for macs2 peaks
$(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak : 	$(DIR_PULL_BOWTIE2)/%_pulldown_trimmed.fq.gz_aligned.bam \\
														$(DIR_PULL_BOWTIE2)/%_input_pulldown_trimmed.fq.gz_aligned.bam
	macs2 callpeak -t $(word 1, $^) -c $(word 2, $^) -f BAM -g hs --name $(patsubst %_peaks.narrowPeak,%,$(@F)) --outdir $(@D)
'

cat(make_var_pull_samp_prefix, file = file_make, sep = '\n', append = TRUE)
cat(make_var_pull_samp, file = file_make, sep = '\n', append = TRUE)
cat(make_rule_pull_samp, file = file_make, sep = '\n', append = TRUE)

#######################################
# PBS script
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

for(i in 1:nrow(pulldown_samples_noinput)) {
	# trackDb.txt entry for MACS2 output
	trackEntry = c(
	  sprintf('track %s_peaks', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('parent %s_sample', pulldown_samples_noinput[i,'humanID']),
	  sprintf('bigDataUrl %s_macs2_peaks.bb', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('shortLabel %s_peaks', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('longLabel %s_macs2_peaks', pulldown_samples_noinput[i,'fullHumanID']),
	  'visibility dense',
	  'type bigBed',
	  'priority 1.5',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)

	# trackDb.txt entry for pulldown simple classification output
	trackEntry = c(
	  sprintf('track %s_simple_class', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('parent %s_sample', pulldown_samples_noinput[i,'humanID']),
	  sprintf('bigDataUrl %s_simple_classification.bb', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('shortLabel %s_simp_class', pulldown_samples_noinput[i,'fullHumanID']),
	  sprintf('longLabel %s_simple_classification', pulldown_samples_noinput[i,'fullHumanID']),
	  'visibility pack',
	  'itemRgb on',
	  'type bigBed 9 .',
	  'priority 1.2',
	  ' ')
	cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
}



}

################################################################################

################################################################################
# MAKEFILE: sample_classification rules

if(bool_bis_samp || bool_pull_samp) {

make_rule_sample_class_bis_module = '
# Intermediates for the bisulfite piece
.PRECIOUS : $(DIR_BIS_BISMARK)/%_bisulfite_highmeth.txt $(DIR_BIS_BISMARK)/%_bisulfite_lowmeth.txt $(DIR_BIS_BISMARK)/%_bisulfite_nometh_signal.txt $(DIR_BIS_BISMARK)/%_bisulfite_nometh_nosignal.txt
$(DIR_BIS_BISMARK)/%_bisulfite_highmeth.txt $(DIR_BIS_BISMARK)/%_bisulfite_lowmeth.txt $(DIR_BIS_BISMARK)/%_bisulfite_nometh_signal.txt $(DIR_BIS_BISMARK)/%_bisulfite_nometh_nosignal.txt : $(DIR_BIS_BISMARK)/%_bisulfite_trimmed.fq.gz_bismark_bt2.CpG_report.txt
	awk -f ../../scripts/bisulfite_sample_module.awk $<
'

make_rule_sample_class_pull_module = '
# Intermediates for the pulldown piece
.PRECIOUS : $(DIR_PULL_MACS)/%_pulldown_peak.txt $(DIR_PULL_MACS)/%_pulldown_nopeak_signal.txt $(DIR_PULL_MACS)/%_pulldown_nopeak_nosignal.txt

$(DIR_PULL_MACS)/%_pulldown_peak.txt : $(DIR_PULL_MACS)/%_pulldown_macs2_peaks.narrowPeak
	awk -v OFS="\\t" \'{print $$1, $$2, $$3}\' $< \\
	| sort -T . -k1,1 -k2,2n \\
	> $@
$(DIR_PULL_MACS)/%_pulldown_nopeak_signal.txt : $(DIR_PULL_MACS)/%_pulldown_nopeak.txt $(DIR_PULL_MACS)/%_pulldown_signal.txt
	bedtools intersect -a $(word 1, $^) -b $(word 2, $^) \\
	| sort -T . -k1,1 -k2,2n \\
	> $@
$(DIR_PULL_MACS)/%_pulldown_nopeak_nosignal.txt : $(DIR_PULL_MACS)/%_pulldown_nopeak.txt $(DIR_PULL_MACS)/%_pulldown_nosignal.txt
	bedtools intersect -a $(word 1, $^) -b $(word 2, $^) \\
	| sort -T . -k1,1 -k2,2n \\
	> $@

.INTERMEDIATE : $(DIR_PULL_MACS)/%_pulldown_nopeak.txt
$(DIR_PULL_MACS)/%_pulldown_nopeak.txt : $(DIR_PULL_MACS)/%_pulldown_peak.txt
	bedtools complement -g <(sort -T . -k1,1 $(CHROM_PATH)) -i $< \\
	| sort -T . -k1,1 -k2,2n \\
	> $@

.INTERMEDIATE : $(DIR_PULL_MACS)/%_pulldown_signal.txt
$(DIR_PULL_MACS)/%_pulldown_signal.txt : $(DIR_PULL_COVERAGES)/%_input_pulldown_merged_coverage.bdg
	cp $< $@

.INTERMEDIATE : $(DIR_PULL_MACS)/%_pulldown_nosignal.txt
$(DIR_PULL_MACS)/%_pulldown_nosignal.txt : $(DIR_PULL_MACS)/%_pulldown_signal.txt
	bedtools complement -g <(sort -T . -k1,1 $(CHROM_PATH)) -i $< \\
	| sort -T . -k1,1 -k2,2n \\
	> $@
'

make_var_sample_class_prefix = sprintf('
################################################################################
# Workflow for sample_classification
SAMPLE_CLASS_PREFIXES := %s', paste(unique(samples$humanID), collapse=' '))
cat(make_var_sample_class_prefix, file = file_make, sep = '\n', append = TRUE)

# The sample class type depends on the type of samples present
if(bool_bis_samp && bool_pull_samp) {
	sample_class_type = 'hybrid_sample_classification'
	sample_class_target = '$(DIR_CLASS_SAMPLE)/%_sample_classification.bed : 	$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_highmeth.txt \\
														$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_lowmeth.txt \\
														$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_nometh_signal.txt \\
														$(DIR_BIS_BISMARK)/%_mc_hmc_bisulfite_nometh_nosignal.txt \\
														$(DIR_PULL_MACS)/%_hmc_pulldown_peak.txt \\
														$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_signal.txt \\
														$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_nosignal.txt'
	rule1 = make_rule_sample_class_bis_module
	rule2 = make_rule_sample_class_pull_module
	class_script = '../../scripts/classify_hybrid_sample.sh'
} else if (bool_bis_samp && !bool_pull_samp) {
	############################################################
	# NOTE: THIS IS NOT EXPLICITLY SUPPORTED RIGHT NOW
	############################################################
	sample_class_type = 'bisulfite_sample_classification'
	sample_class_target = '$(DIR_CLASS_SAMPLE)/%_sample_classification.bed : 	$(DIR_BIS_BISMARK)/%_mc_bisulfite_highmeth.txt \\
														$(DIR_BIS_BISMARK)/%_mc_bisulfite_lowmeth.txt \\
														$(DIR_BIS_BISMARK)/%_mc_bisulfite_nometh_signal.txt \\
														$(DIR_BIS_BISMARK)/%_mc_bisulfite_nometh_nosignal.txt \\
														$(DIR_BIS_BISMARK)/%_hmc_bisulfite_highmeth.txt \\
														$(DIR_BIS_BISMARK)/%_hmc_bisulfite_lowmeth.txt \\
														$(DIR_BIS_BISMARK)/%_hmc_bisulfite_nometh_signal.txt \\
														$(DIR_BIS_BISMARK)/%_hmc_bisulfite_nometh_nosignal.txt'
	rule1 = make_rule_sample_class_bis_module
	rule2 = ''
	class_script = '../../scripts/classify_bisulfite_sample.sh'
} else {
	sample_class_type = 'pulldown_sample_classification'
	sample_class_target = '$(DIR_CLASS_SAMPLE)/%_sample_classification.bed : 	$(DIR_PULL_MACS)/%_mc_pulldown_peak.txt \\
														$(DIR_PULL_MACS)/%_mc_pulldown_nopeak_signal.txt \\
														$(DIR_PULL_MACS)/%_mc_pulldown_nopeak_nosignal.txt \\
														$(DIR_PULL_MACS)/%_hmc_pulldown_peak.txt \\
														$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_signal.txt \\
														$(DIR_PULL_MACS)/%_hmc_pulldown_nopeak_nosignal.txt'
	rule1 = make_rule_sample_class_pull_module
	rule2 = ''
	class_script = '../../scripts/classify_pulldown_sample.sh'
}

make_rule_class_sample = sprintf('
# Master rule
.PHONY : sample_classification
sample_classification : %s

# Generated based on the experiment setup
.PHONY : %s
%s : 	$(patsubst %%,$(DIR_TRACK)/%%_sample_classification.bb,$(SAMPLE_CLASS_PREFIXES)) \\
		$(patsubst %%,$(DIR_CLASS_SAMPLE)/%%_sample_classification.bed,$(SAMPLE_CLASS_PREFIXES))

# Rule for sample classification bigBed
$(DIR_TRACK)/%%_sample_classification.bb : $(DIR_CLASS_SAMPLE)/%%_sample_classification.bed
	bedToBigBed $^ $(CHROM_PATH) $@

# NOTE: There is a known bug in make that incorrectly determines implicit intermediate
# files when they occur in a list of multiple targets and prerequisites.
# https://savannah.gnu.org/bugs/index.php?32042
# The easiest workaround is to make them precious and remove them

# Classification BED
.PRECIOUS : $(DIR_CLASS_SAMPLE)/%%_sample_classification.bed
%s
	bash %s $(CHROM_PATH) $@ $^
	rm -f $^

%s
%s',
	sample_class_type, sample_class_type, sample_class_type, sample_class_target, class_script, rule1, rule2)
cat(make_rule_class_sample, file = file_make, sep = '\n', append = TRUE)

#######################################
# PBS script
pulldown_sample_q = c(
	'#!/bin/bash',
	'#### Begin PBS preamble',
	'#PBS -N class_sample',
	'#PBS -l procs=4,mem=48gb,walltime=6:00:00',
	'#PBS -A sartor_lab',
	'#PBS -q first',
	'#PBS -M rcavalca@umich.edu',
	'#PBS -m abe',
	'#PBS -j oe',
	'#PBS -V',
	'#### End PBS preamble',
	'# Put your job commands after this line',
	sprintf('cd ~/latte/mint/projects/%s/',project),
	'make -j 4 sample_classification')
cat(pulldown_sample_q, file=sprintf('projects/%s/classify_sample.q', project), sep='\n')

for(sample in unique(samples$humanID)) {
# trackDb.txt entry for sample classification
trackEntry = c(
  sprintf('track %s_sample_classification', sample),
  sprintf('parent %s_sample', sample),
  sprintf('bigDataUrl %s_sample_classification.bb', sample),
  sprintf('shortLabel %s_sample_class', sample),
  sprintf('longLabel %s_sample_classification', sample),
  'visibility pack',
  'itemRgb on',
  'type bigBed 9 .',
  'priority 1.1',
  ' ')
cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
}

}

################################################################################

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
		mc = bisulfite_comparisons[i,'mc']
		hmc = bisulfite_comparisons[i,'hmc']
		group = bisulfite_comparisons[i,'group']
		fullHumanID = bisulfite_comparisons[i,'fullHumanID']

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

		groupA = subset(bisulfite_samples,
			grepl(groups[1], bisulfite_samples$group) &
			bisulfite_samples$pulldown == pull &
			bisulfite_samples$bisulfite == bis &
			bisulfite_samples$mc == mc &
			bisulfite_samples$hmc == hmc &
			bisulfite_samples$input == 0)
		groupB = subset(bisulfite_samples,
			grepl(groups[2], bisulfite_samples$group) &
			bisulfite_samples$pulldown == pull &
			bisulfite_samples$bisulfite == bis &
			bisulfite_samples$mc == mc &
			bisulfite_samples$hmc == hmc &
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
	cat(bisulfite_compare_q, file=sprintf('projects/%s/bisulfite_compare.q', project), sep='\n')
}

################################################################################

################################################################################
# MAKEFILE: pulldown_compare rules

if(bool_pull_comp) {
	# Keep track of the compares for the master make rules
	pulldown_compares = c()

	for(i in 1:nrow(pulldown_comparisons)) {
		# Establish row variables
		projectID = pulldown_comparisons[i,'projectID']
		sampleID = pulldown_comparisons[i,'sampleID']
		humanID = pulldown_comparisons[i,'humanID']
		pull = pulldown_comparisons[i,'pulldown']
		bis = pulldown_comparisons[i,'bisulfite']
		mc = pulldown_comparisons[i,'mc']
		hmc = pulldown_comparisons[i,'hmc']
		input = pulldown_comparisons[i,'input']
		group = pulldown_comparisons[i,'group']
		fullHumanID = pulldown_comparisons[i,'fullHumanID']

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

		groupA = subset(pulldown_samples,
			grepl(groups[1], pulldown_samples$group) &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc &
			pulldown_samples$hmc == hmc &
			pulldown_samples$input == 0)
		groupB = subset(pulldown_samples,
			grepl(groups[2], pulldown_samples$group) &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc &
			pulldown_samples$hmc == hmc &
			pulldown_samples$input == 0)

		inputGroupA = subset(pulldown_samples,
			grepl(groups[1], pulldown_samples$group) &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc &
			pulldown_samples$hmc == hmc &
			pulldown_samples$input == 1)
		inputGroupB = subset(pulldown_samples,
			grepl(groups[2], pulldown_samples$group) &
			pulldown_samples$pulldown == pull &
			pulldown_samples$bisulfite == bis &
			pulldown_samples$mc == mc &
			pulldown_samples$hmc == hmc &
			pulldown_samples$input == 1)

		########################################################################
		# Setup variables to put into the makefile

		# For the PePr call from projects/project/pepr_peaks/
		var_input1 = paste(sprintf('../bowtie2_bams/%s_trimmed.fq.gz_aligned.bam', inputGroupA$fullHumanID), sep='', collapse=',')
		var_input2 = paste(sprintf('../bowtie2_bams/%s_trimmed.fq.gz_aligned.bam', inputGroupB$fullHumanID), sep='', collapse=',')
		var_chip1 = paste(sprintf('../bowtie2_bams/%s_trimmed.fq.gz_aligned.bam', groupA$fullHumanID), sep='', collapse=',')
		var_chip2 = paste(sprintf('../bowtie2_bams/%s_trimmed.fq.gz_aligned.bam', groupB$fullHumanID), sep='', collapse=',')

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
		up_bed = sprintf('$(DIR_PULL_PEPR)/%s__PePr_up_peaks.bed', var_name)
		down_bed = sprintf('$(DIR_PULL_PEPR)/%s__PePr_down_peaks.bed', var_name)
		combined_bed = sprintf('$(DIR_PULL_PEPR)/%s_PePr_combined.bed', var_name)
		bigbed = sprintf('$(DIR_TRACK)/%s_PePr_peaks.bb', var_name)

		########################################################################
		# Variables for the makefile

		# Write the pulldown_compare variables for this comparison
		make_var_pull_compare = c(
			'################################################################################',
			'# Workflow for pulldown_compare',
			sprintf('PULLDOWN_COMPARE_%s_PREREQS := %s %s %s', i, up_bed, bigbed, input_signal),
			sprintf('PULLDOWN_COMPARE_%s_INPUT1 := %s', i, var_input1),
			sprintf('PULLDOWN_COMPARE_%s_INPUT2 := %s', i, var_input2),
			sprintf('PULLDOWN_COMPARE_%s_CHIP1 := %s', i, var_chip1),
			sprintf('PULLDOWN_COMPARE_%s_CHIP2 := %s', i, var_chip2),
			sprintf('PULLDOWN_COMPARE_%s_NAME := %s', i, var_name))
		cat(make_var_pull_compare, file = file_make, sep='\n', append=T)

		# Write the pulldown_compare rule for this comparison
		make_rule_pull_compare = c(
			sprintf('.PHONY : pulldown_compare_%s', i),
			sprintf('pulldown_compare_%s : $(PULLDOWN_COMPARE_%s_PREREQS)', i, i),
			'',
			sprintf('%s : %s %s %s %s', up_bed, var_input1_pre, var_input2_pre, var_chip1_pre, var_chip2_pre),
			'	cd pulldown/pepr_peaks; \\',
			sprintf('	$(PATH_TO_PEPR) --input1=$(PULLDOWN_COMPARE_%s_INPUT1) --input2=$(PULLDOWN_COMPARE_%s_INPUT2) --chip1=$(PULLDOWN_COMPARE_%s_CHIP1) --chip2=$(PULLDOWN_COMPARE_%s_CHIP2) --name=$(PULLDOWN_COMPARE_%s_NAME) $(OPTS_PEPR)', i, i, i, i, i),
			sprintf('%s : %s', down_bed, up_bed),
			'',
			sprintf('%s : %s %s', combined_bed, up_bed, down_bed),
			'	bash ../../scripts/combine_pepr.sh $(word 1,$^) $(word 2,$^) $@',
			'',
			sprintf('%s : %s %s', input_signal, var_merged_input1_pre, var_merged_input2_pre),
			'	cat $^ | sort -T . -k1,1 -k2,2n | bedtools merge -d 20 | sort -T . -k1,1 -k2,2n > $@',
			sprintf('%s : %s', bigbed, combined_bed),
			'	bedToBigBed $^ $(CHROM_PATH) $@',
			'')
		cat(make_rule_pull_compare, file = file_make, sep='\n', append=T)

		# Track the number of pulldown compares
		pulldown_compares = c(pulldown_compares, sprintf('pulldown_compare_%s', i))

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
		'')
	cat(make_rule_master_pull_compare, file = file_make, sep='\n', append=T)

	#######################################
	# PBS script
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

}

################################################################################

################################################################################
# MAKEFILE: compare_classification rules

if(bool_bis_comp || bool_pull_comp) {

make_rule_compare_class_bis_module = '
# Intermediates for the bisulfite piece
.PRECIOUS : $(DIR_BIS_MSIG)/%_bisulfite_DMup.txt $(DIR_BIS_MSIG)/%_bisulfite_DMdown.txt $(DIR_BIS_MSIG)/%_bisulfite_noDM_signal.txt $(DIR_BIS_MSIG)/%_bisulfite_noDM_nosignal.txt

$(DIR_BIS_MSIG)/%_bisulfite_DMup.txt : $(DIR_BIS_MSIG)/%_bisulfite_methylSig.txt
	awk -v OFS="\\t" \'NR > 1 && $$5 < 0.05 && $$7 > 0 { print $$1, $$2, $$3 }\' $< | sort -T . -k1,1 -k2,2n > $@

$(DIR_BIS_MSIG)/%_bisulfite_DMdown.txt : $(DIR_BIS_MSIG)/%_bisulfite_methylSig.txt
	awk -v OFS="\\t" \'NR > 1 && $$5 < 0.05 && $$7 < 0 { print $$1, $$2, $$3 }\' $< | sort -T . -k1,1 -k2,2n > $@

$(DIR_BIS_MSIG)/%_bisulfite_noDM_signal.txt : $(DIR_BIS_MSIG)/%_bisulfite_methylSig.txt
	awk -v OFS="\\t" \'NR > 1 && $$5 > 0.05 { print $$1, $$2, $$3 }\' $< | sort -T . -k1,1 -k2,2n > $@

$(DIR_BIS_MSIG)/%_bisulfite_noDM_nosignal.txt : $(DIR_BIS_MSIG)/%_bisulfite_methylSig.txt
	bedtools complement -i <(awk -v OFS="\\t" \'NR > 1 { print $$1, $$2, $$3 }\' $<) -g <(sort -T . -k1,1 $(CHROM_PATH)) | sort -T . -k1,1 -k2,2n > $@
'

make_rule_compare_class_pull_module = '
# Intermediates for the pulldown piece
.PRECIOUS : $(DIR_PULL_PEPR)/%_pulldown_DMup.txt $(DIR_PULL_PEPR)/%_pulldown_DMdown.txt $(DIR_PULL_PEPR)/%_pulldown_noDM_signal.txt $(DIR_PULL_PEPR)/%_pulldown_noDM_nosignal.txt

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

$(DIR_PULL_PEPR)/%_pulldown_DMup.txt : $(DIR_PULL_PEPR)/%_pulldown_tmp_up.txt $(DIR_PULL_PEPR)/%_pulldown_tmp_down.txt
	bedops --difference $^ > $@

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

$(DIR_PULL_PEPR)/%_pulldown_noDM_signal.txt : $(DIR_PULL_PEPR)/%_pulldown_tmp_disjoint_noDM.txt $(DIR_PULL_PEPR)/%_pulldown_tmp_signal.txt
	bedtools intersect -a $(word 1, $^) -b $(word 2, $^) | sort -T . -k1,1 -k2,2n > $@

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
	compare_class_type = 'hybrid_compare_classification'
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
	compare_class_type = 'bisulfite_compare_classification'
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
	compare_class_type = 'pulldown_compare_classification'
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
compare_classification : %s

# Generated based on the experiment setup
.PHONY : %s
%s : 	$(patsubst %%,$(DIR_TRACK)/%%_compare_classification.bb,$(COMPARE_CLASS_PREFIXES)) \\
		$(patsubst %%,$(DIR_CLASS_COMPARE)/%%_compare_classification.bed,$(COMPARE_CLASS_PREFIXES))

# Rule for compare classification bigBed
$(DIR_TRACK)/%%_compare_classification.bb : $(DIR_CLASS_COMPARE)/%%_compare_classification.bed
	bedToBigBed $^ $(CHROM_PATH) $@

# NOTE: There is a known bug in make that incorrectly determines implicit intermediate
# files when they occur in a list of multiple targets and prerequisites.
# https://savannah.gnu.org/bugs/index.php?32042
# The easiest workaround is to make them precious and remove them

# Classification BED
.PRECIOUS : $(DIR_CLASS_COMPARE)/%%_compare_classification.bed
%s
	bash %s $(CHROM_PATH) $@ $^
	rm -f $^

%s
%s',
	compare_class_type, compare_class_type, compare_class_type, compare_class_target, class_script, rule1, rule2)
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
cat(bisulfite_compare_q, file=sprintf('projects/%s/compare_classification.q', project), sep='\n')

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
