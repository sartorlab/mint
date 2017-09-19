# Configuration header
config_header = sprintf('# Configuration for mint pipeline analyses

# This makefile was generated using mint v0.4.1

################################################################################
# Project and experimental information

# Project name
PROJECT = %s
# Genome
GENOME = %s

################################################################################
# Genome paths
# Genomes, like as downloaded from iGenomes, will also contain bisulfite
# converted genomes for
GENOME_PATH := %s
# Location of bowtie2 indexes
BOWTIE2_GENOME_PATH := %s/genome
# Location of chromosome length file
CHROM_PATH := %s',
	project, genome, genomepath, genomepath, chrompath)
cat(config_header, file = 'config_header.mk', sep='\n')

################################################################################
# INITIALIZATION: Directory creation

# Always create these folders
setup_commands = c(
	sprintf('mkdir projects/%s', project),
	sprintf('mkdir projects/%s/tmp', project),
	sprintf('mkdir projects/%s/pbs_jobs', project),
	sprintf('mkdir projects/%s/data', project),
	sprintf('mkdir projects/%s/data/raw_fastqs', project),
	sprintf('mkdir projects/%s/summary', project),
	sprintf('mkdir projects/%s/summary/{figures,tables,reports}', project),
	sprintf('mkdir projects/%s/summary/reports/multiqc', project),
	sprintf('mkdir projects/%s/%s_hub', project, project),
	sprintf('mkdir projects/%s/%s_hub/%s', project, project, genome),
	sprintf('mkdir projects/%s/classifications', project),
	sprintf('mkdir projects/%s/classifications/{simple,sample}', project),
	sprintf('mkdir projects/%s/RData', project)
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
		sprintf('mkdir projects/%s/bisulfite/dss', project)
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
		sprintf('mkdir projects/%s/pulldown/csaw', project)
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
	sprintf('cat config_header.mk template_config.mk > projects/%s/config.mk', project, project),
	sprintf('rm config_header.mk', project),
	sprintf('cp narrowPeak.as projects/%s/', project),
	sprintf('cp projects/%s_samples.txt projects/%s/data', project, project),
	sprintf('cp projects/%s_comparisons.txt projects/%s/data', project, project)
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
