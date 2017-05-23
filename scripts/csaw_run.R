library(optparse)

######
###### Options
######

option_list = list(
	make_option('--project', type='character', help='(Required) Name of the project the analysis is a part of'),
	make_option('--chipfiles', type='character', help='(Required) Comma-delimited list of BAM files representing aligned ChIP reads.'),
	make_option('--inputfiles', type='character', help='(Optional) Comma-delimited list of BAM files representing aligned input reads.'),
	make_option('--chipnames', type='character', help='(Required) Comma-delimited list of sample names to give the chipfiles.'),
	make_option('--fraglength', type='integer', help='(Required) An integer scalar containing the average length(s) of the sequenced fragmentsin each library.'),
	make_option('--winwidth', type='integer', help='(Required) An integer scalar specifying the width of the windows to tile the genome.'),
	make_option('--winspacing', type='integer', help='(Required) An integer scalar specifying how much the start of each window should be shifted to the left.'),
	make_option('--useinput', type='logical', help='(Required) TRUE/FALSE, specifying if inputfiles should be used to filter windows based on counts and fold change.'),
	make_option('--prior.count', type='numeric', help='(Required if useinput=TRUE) A numeric scalar specifying the minimum number of counts required in a window before filtering on fold change.'),
	make_option('--chipfold', type='numeric', help='(Required if useinput=TRUE) A numeric scalar specifying the fold change (log2(chipfold)) below which windows will be removed from analysis.'),
	make_option('--model', type='character', help='(Required) A formula (e.g. ~1+group) to be used to create a design matrix for the analysis.'),
	make_option('--groups', type='character', help='(Required) A comma-delimited list of the groups to which chipfiles belong. (e.g. if chipfiles has four files, 1,1,0,0)'),
	make_option('--contrast', type='character', help='(Required) A comma-delimited list as long as the number of components in model, specifying the contrast to test.'),
	make_option('--covariates', type='character', help='(Optional) Encoded covariates to use in the design matrix. E.g. "subject:1,2,1,2;cont:0.6,0.1,0.4,0.3". Use NA if none.'),
	make_option('--covIsNumeric', type='character', help='(Required if covariates!= NA) A comma-delimited list of 0s and 1s indicating if the covariates are continuous (1) or a factor (0).'),
	make_option('--mergewithin', type='numeric', help='(Required) A numeric scalar specifying the maximum distance between adjacent windows when combining windows after testing.'),
	make_option('--maxmerged', type='numeric', help='(Required) A numeric scalar specifying the maximum size of merged intervals.'),
	make_option('--FDRthreshold', type='numeric', help='(Required) A numeric scalar specifying the required FDR to be considered "significant" and to be returned in the outprefix_csaw_significant.txt table.'),
	make_option('--pvalthreshold', type='numeric', help='(Required) A numeric scalaer specifying the alternative p-value threshold for significance if no results FDR < 0.25.'),
	make_option('--interpretation', type='character', help='(Required) A comma-delimited list of how to interpret logFC < 0 (first entry) and logFC >= 0 (second entry).'),
	make_option('--quiet', type='logical', default=FALSE, help='TRUE/FALSE indicating whether progress messages should be printed.'),
	make_option('--outprefix', type='character', help='(Required) A string inicating the prefix of the output file names.')
)

opt = parse_args(OptionParser(option_list=option_list))

use_input = opt$useinput

chip_files = opt$chipfiles
if(use_input) {
	input_files = opt$inputfiles
} else {
	input_files = NA
}
chip_names = opt$chipnames

fragment_length = opt$fraglength
window_width = opt$winwidth
window_spacing = opt$winspacing

prior.count = opt$prior.count
chipfold = opt$chipfold

model = opt$model
groups = opt$groups
contrast = opt$contrast
covariates = opt$covariates
covIsNumeric = opt$covIsNumeric

mergewithin = opt$mergewithin
maxmerged = opt$maxmerged
FDRthreshold = opt$FDRthreshold
pvalthreshold = opt$pvalthreshold
interpretation = opt$interpretation

quiet = opt$quiet
prefix = opt$outprefix

library(csaw)
library(edgeR)

# # Things from opt that are hardcoded right now
# use_input = TRUE
#
# chip_files = 'pulldown/bowtie2_bams/IDH2mut_1_hmc_pulldown_trimmed.fq.gz_aligned.bam,pulldown/bowtie2_bams/IDH2mut_2_hmc_pulldown_trimmed.fq.gz_aligned.bam,pulldown/bowtie2_bams/NBM_1_hmc_pulldown_trimmed.fq.gz_aligned.bam,pulldown/bowtie2_bams/NBM_2_hmc_pulldown_trimmed.fq.gz_aligned.bam'
# if(use_input) {
# 	input_files = 'pulldown/bowtie2_bams/IDH2mut_1_hmc_input_pulldown_trimmed.fq.gz_aligned.bam,pulldown/bowtie2_bams/IDH2mut_2_hmc_input_pulldown_trimmed.fq.gz_aligned.bam,pulldown/bowtie2_bams/NBM_1_hmc_input_pulldown_trimmed.fq.gz_aligned.bam,pulldown/bowtie2_bams/NBM_2_hmc_input_pulldown_trimmed.fq.gz_aligned.bam'
# } else {
# 	input_files = NA
# }
# chip_names = 'IDH2mut_1,IDH2mut_2,NBM_1,NBM_2'
#
# fragment_length = 110
# window_width = 100
# window_spacing = 50
#
# prior.count = 5
# chipfold = 2
#
# model = '~ 1 + group'
# groups = '1,1,0,0'
# contrast = '0,1'
# covariates = 'NA'
# covIsNumeric = 0
#
# mergewithin = 500
# maxmerged = 2000
# FDRthreshold = 0.05
# interpretation = 'IDH2mut,NBM'
#
# quiet = TRUE
# prefix = 'IDH2mut_v_NBM_hmc_pulldown'

######
###### Parse arguments
######

chip_files = unlist(strsplit(chip_files, ','))
if(use_input) {
	input_files = unlist(strsplit(input_files, ','))
	all_files = c(chip_files, input_files)
} else {
	all_files = chip_files
}
chip_names = unlist(strsplit(chip_names, ','))
groups = unlist(strsplit(groups, ','))

if(!all(is.na(covIsNumeric))) {
	covIsNumeric = unlist(strsplit(covIsNumeric, ','))
}

if(covariates != 'NA') {
	covs = unlist(strsplit(covariates, ';'))

	covs_colnames = sapply(covs, function(cov){ unlist(strsplit(cov, ':'))[1] }, USE.NAMES = FALSE)
	covs_values = sapply(covs, function(cov){ unlist(strsplit(cov, ':'))[2] }, USE.NAMES = FALSE)

	cov_df = as.data.frame(sapply(covs_values, function(cov_value) {unlist(strsplit(cov_value, ',')) }), stringsAsFactors=F)
	colnames(cov_df) = covs_colnames

	for(i in seq_along(colnames(cov_df))) {
		if(covIsNumeric[i] == 1) {
			cov_df[,i] = as.numeric(cov_df[,i])
		} else {
			cov_df[,i] = as.factor(cov_df[,1])
		}
	}
	use_cov = TRUE
} else {
	cov_df = NA
	use_cov = FALSE
}

interpretation = unlist(strsplit(interpretation, ','))
contrast = as.integer(unlist(strsplit(contrast, ',')))

######
###### Read
######

# Read all in and create subsets
message('Reading in all data')
all = windowCounts(all_files, ext = fragment_length, width = window_width, spacing = window_spacing)
message('Splitting data')
chip  = all[, seq_along(chip_files)]
if(use_input) {
	input  = all[, seq(from = length(chip_files) + 1, to = length(all_files), by = 1)]
} else {
	input = NA
}

######
###### Window filtering with the inputs
######

if(use_input) {
	message('Determining filtering by input')
	all_binned = windowCounts(all_files, bin = TRUE, width = 10000)
	chip_binned = all_binned[, seq_along(chip_files)]
	input_binned = all_binned[, seq(from = length(chip_files) + 1, to = length(all_files), by = 1)]

	filter_stat = filterWindows(data = chip, background = input, type = 'control', prior.count = prior.count,
		norm.fac = list(chip_binned, input_binned))

	keep = (filter_stat$filter > log2(chipfold))
	chip = chip[keep,]
}

# Compute normalization
message('Computing offsetes')
chip_offsets = normOffsets(chip)

######
###### Setup data for differential binding test
######

message('Converting to DGEList')
y = asDGEList(chip, norm.factors = chip_offsets, group = groups)
rownames(y$samples) = chip_names
colnames(y$counts) = chip_names
if(use_cov) {
	y$samples = cbind(y$samples, cov_df)
}

######
###### Setup design matrix and contrasts
######

design = model.matrix(as.formula(model), data = y$samples)

######
###### Estimate dispersion
######

message('Estimating dispersion')
y = estimateDisp(y, design)

######
###### MDS and BCV plots
######

message('Plotting MDS plot')
mds_file = sprintf('pulldown/csaw/%s_QC_MDSplot.pdf', prefix)
pdf(file = mds_file, height = 6, width = 6)
	plotMDS(y, cex = 0.5)
dev.off()

message('Plotting BCV plot')
bcv_file = sprintf('pulldown/csaw/%s_QC_BCVplot.pdf', prefix)
pdf(file = bcv_file, height = 6, width = 6)
	plotBCV(y)
dev.off()

######
###### Do fit and test
######

message('Computing fit')
fit = glmQLFit(y, design, robust=TRUE)

qld_file = sprintf('pulldown/csaw/%s_QC_QLDispplot.pdf', prefix)
pdf(file = qld_file, height = 6, width = 6)
	plotQLDisp(fit)
dev.off()

message('Computing test results')
results = glmQLFTest(fit, contrast=contrast)

######
###### Window merging and multiple test correction
######

message('Merging windows')
merged = mergeWindows(rowRanges(chip), tol = mergewithin, max.width = maxmerged)
combined_tests = combineTests(merged$id, results$table)

# Return the average logFCs of the windows contained in each merged window
combined_tests$logFC = sapply(split(results$table$logFC, merged$id), mean)

# Return the direction based on the interpretation vector
combined_tests$direction = ifelse(combined_tests$logFC < 0, interpretation[1], interpretation[2])

# Assign a color to each of the interpretation vectors
combined_tests$color = ifelse(combined_tests$logFC < 0, '0,0,255', '102,102,255')

# Make a GRanges with mcols of combined_tests
combined_gr = merged$region
mcols(combined_gr) = combined_tests[,c('nWindows', 'logFC.up', 'logFC.down', 'logFC', 'direction', 'PValue', 'FDR', 'color')]

# Convert to a data.frame and create the different variants
combined_df = data.frame(combined_gr, stringsAsFactors=F)
colnames(combined_df) = c('chr','start','end','width','strand','nWindows','logFC.up','logFC.down','logFC','direction','PValue','FDR','color')
combined_df = combined_df[, c('chr','start','end','strand','nWindows','logFC.up','logFC.down','logFC','direction','PValue','FDR','color')]

# Significant regions
significant_df = subset(combined_df, FDR <= FDRthreshold)
increment = 0.05
while(nrow(significant_df) == 0 && FDRthreshold < 0.25) {
	FDRthreshold = FDRthreshold + increment
	significant_df = subset(combined_df, FDR <= FDRthreshold)
}

if(nrow(significant_df) == 0) {
	message(sprintf('No tests are FDR <= %s, using alternative p-value threshold.', FDRthreshold))
	significant_df = subset(combined_df, PValue <= pvalthreshold)

	if(nrow(significant_df) == 0) {
		stop(sprintf('No differential methylation at FDR <= %s or p-value <= %s!', FDRthreshold, pvalthreshold))
	}
} else {
	sprintf('Using FDR <= %s', FDRthreshold)
}


# For annotatr
annotatr_df = significant_df[, c('chr','start','end','direction','logFC','strand','PValue')]
annotatr_df$start = format(annotatr_df$start - 1, scientific = FALSE)
annotatr_df$end = format(annotatr_df$end, scientific = FALSE)

# For bigBed
bigbed_df = significant_df[, c('chr','start','end','direction','start','end','color')]
bigbed_df$start = format(bigbed_df$start - 1, scientific = FALSE)
bigbed_df$end = format(bigbed_df$end, scientific = FALSE)
bigbed_df$score = '1000'
bigbed_df$strand = '.'
bigbed_df = bigbed_df[, c('chr','start','end','direction','score','strand','start.1','end.1','color')]

######
###### Save results meeting FDR threshold
######

message('Saving all results')
all_results_file = sprintf('pulldown/csaw/%s_csaw_all.txt', prefix)
write.table(combined_df, file = all_results_file, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)

message('Saving significant results')
sig_results_file = sprintf('pulldown/csaw/%s_csaw_significant.txt', prefix)
write.table(significant_df, file = sig_results_file, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)

message('Saving significant results for annotatr')
annotatr_file = sprintf('pulldown/csaw/%s_csaw_for_annotatr.txt', prefix)
write.table(annotatr_df, file = annotatr_file, sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

message('Saving significant results for bigBed')
bigbed_file = sprintf('pulldown/csaw/%s_csaw_for_bigBed.bed', prefix)
write.table(bigbed_df, file = bigbed_file, sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
