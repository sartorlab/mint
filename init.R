# Working directory is mint/

################################################################################
# Parse arguments
library(optparse)

option_list = list(
    make_option('--project', type='character',
        help='[Required] The name for the project. Corresponding to projects/[project]_annotation.txt file.'),
    make_option('--genome', type='character',
        help='[Required] The genome to be used for analyses.'),
    make_option('--genomepath', type='character',
        help='[Required] The absolute path to the reference genome.'),
    make_option('--chrompath', type='character',
        help='[Required] The absolute path to the chromosome length file.'),
    make_option('--datapath', type='character',
        help='[Required] The absolute path to the .fastq.gz sequencing data whose filenames match those in the sampleID column of the projects/[project]_annotation.txt file.')
)
option_parser = OptionParser(
    usage = "usage: Rscript %prog [options]",
    option_list = option_list,
    add_help_option = TRUE,
    prog = "init.R")
opt = parse_args(option_parser)

project = opt$project
genome = opt$genome
genomepath = opt$genomepath
chrompath = opt$chrompath
datapath = opt$datapath

file_make = sprintf('projects/%s/makefile', project)
file_config = sprintf('projects/%s/config.mk', project)
dir_hub = sprintf('projects/%s/%s_hub', project, project)
dir_track = sprintf('%s/%s', dir_hub, genome)

################################################################################
################################################################################
################################################################################

################################################################################
# READ project information
annots = read.table(file = sprintf('projects/%s_samples.txt', project), header = T, sep = '\t', stringsAsFactors = F)
comparisons = read.table(file = sprintf('projects/%s_comparisons.txt', project), header = T, sep = '\t', stringsAsFactors = F)

# Add the fullHumanID column to the annots table, which will be the prefix for all created files
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

# Add the fullHumanID column to the comparisons table, which will be the prefix for all created files
comparisons$fullHumanID = NA
for(i in 1:nrow(comparisons)) {
    comparison = comparisons[i,'comparison']
    pull = comparisons[i,'pulldown']
    bis = comparisons[i,'bisulfite']
    mc = comparisons[i,'mc']
    hmc = comparisons[i,'hmc']
    input = comparisons[i,'input']

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

    comparisons[i,'fullHumanID'] = sprintf("%s_%s_%s", comparison, mark, platform)
}

# Split samples by bisulfite and pulldown
bisulfite = subset(annots, bisulfite == 1)
pulldown = subset(annots, pulldown == 1)

# Determine the sample
samples = subset(annots, !grepl('comparison', sampleID))

# Split by sample and comparisons
bisulfite_samples = subset(bisulfite, !grepl('comparison', sampleID))
bisulfite_comparisons = subset(comparisons, bisulfite == 1)

pulldown_samples = subset(pulldown, !grepl('comparison', sampleID))
pulldown_samples_noinput = subset(pulldown, !grepl('comparison', sampleID) & input == 0)
pulldown_samples_input = subset(pulldown, !grepl('comparison', sampleID) & input == 1)
pulldown_comparisons = subset(comparisons, pulldown == 1)

# Control variables
bool_bis_samp = nrow(bisulfite_samples) > 0
bool_pull_samp = nrow(pulldown_samples) > 0
bool_bis_comp = nrow(bisulfite_comparisons) > 0
bool_pull_comp = nrow(pulldown_comparisons) > 0

################################################################################
# READ group names
# if(bool_bis_comp || bool_pull_comp) {
#     group_file = sprintf('projects/%s_comparisons.txt', project)
#     if(!file.exists(group_file)) {
#         stop(sprintf('You must include a tab-delimited file named %s_groups.txt in the projects/ folder. It should have two columns, named group and name, giving the group number and the name of the group.', project))
#     } else {
#         group_names = read.table(group_file, sep='\t', header = TRUE, quote = '', comment.char = '', stringsAsFactors= FALSE)
#     }
# }

# NOTE: ADD ERROR CHECKING
# 1. Are there input samples for pulldowns?
# 2. If there are comparisons, do the group numbers correspond to the groups assigned to the samples?

# NOTE: Have to think about how to handle inputs shared between mc and hmc pulldowns.
# One option is to create an mc and hmc symlink to the same input file, but then all
# analysis would be doubled for that input file. This is the easiest path right now.

################################################################################

source('scripts/sub_init_setup.R')

source('scripts/sub_init_hub.R')

source('scripts/sub_init_bis_align.R')

source('scripts/sub_init_pull_align.R')

source('scripts/sub_init_bis_sample.R')

source('scripts/sub_init_pull_sample.R')

source('scripts/sub_init_sample_class.R')

source('scripts/sub_init_bis_compare.R')

source('scripts/sub_init_pull_compare.R')

source('scripts/sub_init_compare_class.R')

source('scripts/sub_init_clean.R')
