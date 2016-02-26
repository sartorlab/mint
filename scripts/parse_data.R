# Working directory is mint/

library(optparse)

# Parse arguments
option_list = list(
  make_option('--project', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

project = opt$project

########################################################
# Read annotations

annots = read.table(file=sprintf('projects/%s/data/%s_annotation.txt', project, project), header=T, sep='\t', stringsAsFactors=F)

# Split by sample and comparisons
sample_annots = subset(annots, !grepl('comparison', annots$sampleID))
comparison_annots = subset(annots, grepl('comparison', annots$sampleID))

########################################################
# Deal with samples

pullAlignFiles = c()
bisAlignFiles = c()
pullSampleFiles = c()
for(i in 1:nrow(sample_annots)) {
	# Establish row variables
	projectID = sample_annots[i,'projectID']
	sampleID = sample_annots[i,'sampleID']
	humanID = sample_annots[i,'humanID']
	pulldown = sample_annots[i,'pulldown']
	bisulfite = sample_annots[i,'bisulfite']
	mc = sample_annots[i,'mc']
	hmc = sample_annots[i,'hmc']
	input = sample_annots[i,'input']
	group = sample_annots[i,'group']

	# Sample file
	sampleFile = sprintf("%s/projects/%s/data/raw_fastqs/%s.fastq.gz", getwd(), project, sampleID)

	# Establish platform for row
	if( pulldown == 1 ) {
	  platform = "pulldown"
	} else if ( bisulfite == 1 ) {
	  platform = "bisulfite"
	}

	# Establish mark for row
	if( mc == 1 && hmc == 1 ) {
	  mark = "mc_hmc"
	} else if ( mc == 1 && hmc == 0 ) {
	  mark = "mc"
	} else if ( mc == 0 && hmc == 1 ) {
	  mark = "hmc"
	}

	# Make symlink on the basis of input
	if ( input == 1 ) {
	  readableFile = sprintf("projects/%s/%s/raw_fastqs/%s_%s_input_%s.fastq.gz", project, platform, humanID, mark, platform)
	} else {
	  readableFile = sprintf("projects/%s/%s/raw_fastqs/%s_%s_%s.fastq.gz", project, platform, humanID, mark, platform)
	}

	# symlink data/raw_fastq to pulldown/raw_fastq or bisulfite/raw_fastq
	command = sprintf('ln -s %s %s', sampleFile, readableFile)
	message(sprintf('Creating symlink: %s', command))
	system(command)

	# Add to the appropriate file group variable
	if ( pulldown == 1 ) {
	  if ( input == 1 ) {
		pullAlignFiles = c(pullAlignFiles, sprintf("%s_%s_input_%s.fastq.gz", humanID, mark, platform))
	  } else {
		pullAlignFiles = c(pullAlignFiles, sprintf("%s_%s_%s.fastq.gz", humanID, mark, platform))
		pullSampleFiles = c(pullSampleFiles, sprintf("%s_%s_%s.fastq.gz", humanID, mark, platform))
	  }
	} else if ( bisulfite == 1 ) {
	  bisAlignFiles = c(bisAlignFiles, sprintf("%s_%s_%s.fastq.gz", humanID, mark, platform))
	}
}

# Write the variables to projects/project/variables.mk
message(sprintf('Writing projects/%s/variables.mk ...', project))
vars = c(
	sprintf('BISULFITE_FASTQ_FILES := %s', paste(bisAlignFiles, collapse=' ')),
	sprintf('PULLDOWN_FASTQ_FILES := %s', paste(pullAlignFiles, collapse=' ')),
	sprintf('PULLDOWN_SAMPLE_FILES := %s', paste(pullSampleFiles, collapse=' '))
	)
cat(vars, file = sprintf('projects/%s/variables.mk', project), sep='\n')
