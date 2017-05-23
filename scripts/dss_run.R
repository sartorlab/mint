library(optparse)

######
###### Options
######

option_list = list(
	make_option('--project', type='character', help='(Required) Name of the project the analysis is a part of.'),
	make_option('--genome', type='character', help='(Required) A supported genome from ??GenomeInfoDb::fetchExtendedChromInfoFromUCSC.'),
	make_option('--files', type='character', help='(Required) Comma-delimited list of input files containing columns chr, pos, N (number of reads), and X (number of methylated reads).'),
	make_option('--samplenames', type='character', help='(Required) Comma-delimited list of sample names to give the files.'),
	make_option('--destrand', type='logical', help='(Required) Combine CpGs that occur sequentially on opposite strands?'),
	make_option('--tilewidth', type='integer', help='(Required) An integer tile width to use to group CpGs into tiles. Use 0 for no tiling.'),
	make_option('--model', type='character', help='(Required) A formula (e.g. ~1+group) to be used to create a design matrix for the analysis.'),
	make_option('--groups', type='character', help='(Required) A comma-delimited list of the groups to which chipfiles belong. (e.g. if chipfiles has four files, 1,1,0,0)'),
	make_option('--contrast', type='character', help='(Required) A comma-delimited list as long as the number of components in model, specifying the contrast to test.'),
	make_option('--covariates', type='character', help='(Optional) Encoded covariates to use in the design matrix. E.g. "subject:1,2,1,2;cont:0.6,0.1,0.4,0.3". Use NA if none.'),
	make_option('--covIsNumeric', type='character', help='(Required if covariates!= NA) A comma-delimited list of 0s and 1s indicating if the covariates are continuous (1) or a factor (0).'),
	make_option('--methdiffthreshold', type='numeric', help='(Required) A numeric scalar between 0 and 100. The threshold above which (in combination with FDRthreshold) a CpG or region to be significant.'),
	make_option('--FDRthreshold', type='numeric', help='(Required) A numeric scalar specifying the required FDR to be considered "significant" and to be returned.'),
	make_option('--pvalthreshold', type='numeric', help='(Required) A numeric scalaer specifying the alternative p-value threshold for significance if no results FDR < 0.25.'),
	make_option('--interpretation', type='character', help='(Required) A comma-delimited list of how to interpret diff < 0 (first entry) and diff >= 0 (second entry).'),
	make_option('--quiet', type='logical', default=FALSE, help='TRUE/FALSE indicating whether progress messages should be printed.'),
	make_option('--outprefix', type='character', help='(Required) A string inicating the prefix of the output file names.')
)

opt = parse_args(OptionParser(option_list=option_list))

files = opt$files
sample_names = opt$samplenames
genome = opt$genome

destrand = opt$destrand
tilewidth = opt$tilewidth

model = opt$model
groups = opt$groups
contrast = opt$contrast
covariates = opt$covariates
covIsNumeric = opt$covIsNumeric

methdiffthreshold = opt$methdiffthreshold
FDRthreshold = opt$FDRthreshold
pvalthreshold = opt$pvalthreshold
interpretation = opt$interpretation

quiet = opt$quiet
prefix = opt$outprefix

library(readr)
library(DSS)
library(bsseq)

# project = 'test_hybrid_small'
# files = 'bisulfite/bismark/IDH2mut_1_mc_hmc_bisulfite_trimmed_bismark_bt2.CpG_report_for_dss.txt,bisulfite/bismark/IDH2mut_2_mc_hmc_bisulfite_trimmed_bismark_bt2.CpG_report_for_dss.txt,bisulfite/bismark/NBM_1_mc_hmc_bisulfite_trimmed_bismark_bt2.CpG_report_for_dss.txt,bisulfite/bismark/NBM_2_mc_hmc_bisulfite_trimmed_bismark_bt2.CpG_report_for_dss.txt'
# sample_names = 'IDH2mut_1,IDH2mut_2,NBM_1,NBM_2'
# genome = 'hg19'
#
# destrand = TRUE
# tilewidth = 50
#
# model = '~group'
# groups = '1,1,0,0'
# contrast = '0,1'
# covariates = 'NA'
# covIsNumeric = 0
#
# methdiffthreshold = 0.05
# FDRthreshold = 0.05
# interpretation = 'NBM,IDH2mut'
#
# quiet = FALSE
# prefix = 'IDH2mut_v_NBM_mc_hmc_bisulfite'

###
### Parse arguments
###

files = unlist(strsplit(files, ','))
sample_names = unlist(strsplit(sample_names, ','))

group = unlist(strsplit(groups, ','))

if(!all(is.na(covIsNumeric))) {
	covIsNumeric = unlist(strsplit(as.character(covIsNumeric), ','))
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


# Need to determine the coefficient of interest from the contrast
# NOTE: DSS is only able to handle testing if a coefficient is 0 because
# of the interface in DMLtest.multiFactor. This consequently limits csaw.
# NOTE: When we add DSS capabilities to methylSig we should generalize this.
coef = which(contrast == 1)

###
### Read data and construct combined bsseq object
###

data_list = as.list(seq_along(files))
for(i in seq_along(files)) {
	tmp = readr::read_tsv(file = files[i], col_types = 'cicii', col_names = c('chr','pos','strand','Cov','M'))
	tmp_bs = BSseq(chr = tmp$chr, pos = tmp$pos, M = as.matrix(tmp$M), Cov = as.matrix(tmp$Cov), sampleNames = sample_names[i])
	strand(tmp_bs) = tmp$strand
	data_list[[i]] = tmp_bs

	rm(tmp)
	rm(tmp_bs)
}
data = Reduce(combine, data_list)
rm(data_list)

###
### Destranding
###

if(destrand) {
	data = strandCollapse(data, shift = TRUE)
}

###
### Tiling
###

if(tilewidth != 0) {
	tiles = tileGenome(GenomeInfoDb::Seqinfo(genome = genome), tilewidth = tilewidth, cut.last.tile.in.chrom = TRUE)

	tiled_M = getCoverage(data, regions = tiles, what = "perRegionTotal", type = 'M')
	tiled_M[is.na(tiled_M)] = 0
	tiled_Cov = getCoverage(data, regions = tiles, what = "perRegionTotal", type = 'Cov')
	tiled_Cov[is.na(tiled_Cov)] = 0

	data = BSseq(gr = tiles, M = tiled_M, Cov = tiled_Cov, rmZeroCov = TRUE)

	rm(tiles)
	rm(tiled_M)
	rm(tiled_Cov)
}

###
### Recover difference of methylation averages between groups
###

sorted_groups = sort(unique(as.integer(group)), decreasing = TRUE)

group1_M = rowSums(getCoverage(data, type='M')[,which(group == sorted_groups[1])])
group0_M = rowSums(getCoverage(data, type='M')[,which(group == sorted_groups[2])])
group1_Cov = rowSums(getCoverage(data, type='Cov')[,which(group == sorted_groups[1])])
group0_Cov = rowSums(getCoverage(data, type='Cov')[,which(group == sorted_groups[2])])

mu1 = (group1_M / group1_Cov)*100
mu0 = (group0_M / group0_Cov)*100
methdiff = mu1 - mu0

rm(group1_M)
rm(group0_M)
rm(group1_Cov)
rm(group0_Cov)

###
### Construct design data.frame
###

design = data.frame(group = group)
if(use_cov) {
	design = cbind(design, cov_df)
}

###
### Do the fit and test
###

fit = DMLfit.multiFactor(data, design = design, formula = as.formula(model))

result = DMLtest.multiFactor(fit, coef=coef)

if(tilewidth != 0) {
	# If tiling, DSS ignores the fact that the ranges of the bsseq object are not CpGs
	result$start = result$pos
	result$end = result$pos + tilewidth - 1
} else {
	result$start = result$pos
	result$end = result$pos
}

result$mu1 = mu1
result$mu0 = mu0
result$methdiff = methdiff
result$direction = ifelse(methdiff < 0, interpretation[1], interpretation[2])
result = result[!is.na(result$stat),]

###
### Prepare for writing
###

# Significant
significant_df = subset(result, fdrs < FDRthreshold & abs(methdiff) > methdiffthreshold)
increment = 0.05
while(nrow(significant_df) == 0 && FDRthreshold < 0.25) {
	FDRthreshold = FDRthreshold + increment
	significant_df = subset(combined_df, fdrs < FDRthreshold)
}

if(nrow(significant_df) == 0) {
	message(sprintf('No tests are FDR <= %s, using alternative p-value threshold.', FDRthreshold))
	significant_df = subset(combined_df, pvals <= pvalthreshold)

	if(nrow(significant_df) == 0) {
		stop(sprintf('No differential methylation at FDR <= %s or p-value <= %s!', FDRthreshold, pvalthreshold))
	}
} else {
	sprintf('Using FDR <= %s', FDRthreshold)
}

# For annotatr
annotatr_df = significant_df
annotatr_df$start = format(annotatr_df$start - 1, scientific = FALSE)
annotatr_df$end = format(annotatr_df$end, scientific = FALSE)
annotatr_df$strand = '*'
annotatr_df = annotatr_df[,c('chr','start','end','direction','pvals','strand','methdiff','mu1','mu0')]

# For bedgraph
bedgraph_df = significant_df
bedgraph_df$start = format(bedgraph_df$start - 1, scientific = FALSE)
bedgraph_df$end = format(bedgraph_df$end, scientific = FALSE)
bedgraph_df = bedgraph_df[, c('chr','start','end','methdiff')]

######
###### Save results meeting FDR threshold
######

message('Saving all results')
all_results_file = sprintf('bisulfite/dss/%s_dss_all.txt', prefix)
write.table(result[, c('chr','start','end','stat','mu0','mu1','methdiff','direction','pvals','fdrs')],
	file = all_results_file, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)

message('Saving significant results')
sig_results_file = sprintf('bisulfite/dss/%s_dss_significant.txt', prefix)
write.table(significant_df[,c('chr','start','end','stat','mu0','mu1','methdiff','direction','pvals','fdrs')],
	file = sig_results_file, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)

message('Saving significant results for annotatr')
annotatr_file = sprintf('bisulfite/dss/%s_dss_for_annotatr.txt', prefix)
write.table(annotatr_df, file = annotatr_file, sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

message('Saving significant results for bedgraph')
bedgraph_file = sprintf('bisulfite/dss/%s_dss_for_bigWig.bedGraph', prefix)
write.table(bedgraph_df, file = bedgraph_file, sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
