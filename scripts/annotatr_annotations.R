library(annotatr)
library(ggplot2)
library(regioneR)
library(optparse)

################################################################################
# Deal with command line options

option_list = list(
	make_option('--file', type='character', help='[Required] Tab-delimited file with genomic locations and possibly associated data.'),
	make_option('--genome', type='character', help='[Required] The shortname for the genome used, e.g. hg19, mm9, rn4.'),
	make_option('--annot_type', type='character', help='[Required] One of bismark, simple_bisulfite_bismark, macs2, simple_pulldown_macs2, sample_class, dss, csaw, simple_pulldown_csaw, or compare_class. Indicates what type of data is being annotated.'),
	make_option('--group1', type='character', help='[Required] A character indicating the name of group1 or NULL.'),
	make_option('--group0', type='character', help='[Required] A character indicating the name of group0 or NULL.')
)
opt = parse_args(OptionParser(option_list=option_list))

file = opt$file
genome = opt$genome
annot_type = opt$annot_type

# If the group names are not NULL, use them
if(!is.null(opt$group1)) {
	group1 = opt$group1
	chip1 = group1
}
if(!is.null(opt$group0)) {
	group0 = opt$group0
	chip0 = group0
}

# Interpret mark
# This order matters very much!
if(grepl('mc_hmc', file)) {
	mark = 'mc_hmc'
} else if (grepl('hmc', file)) {
	mark = 'hmc'
} else if (grepl('mc', file)) {
	mark = 'mc'
}

################################################################################
# Deal with the annot_type
if(annot_type == 'bismark') {
	# bis_align
	# BARPLOT + NUMERICALS
	# Random for barplot

	# head ./bisulfite/bismark/IDH2mut_2_mc_hmc_bisulfite_trimmed_bismark_bt2.bismark.cov
	# chr21	9437517	9437517	86.241610738255	.   .	257
	# chr21	9437519	9437519	75.503355704698	.   .	225
	# chr21	9437531	9437531	58.9225589225589   .   .	175
	# chr21	9437533	9437533	86.5771812080537   .   .	258
	extraCols = c(coverage = 'integer')
	rename_score = 'perc_meth'
	rename_name = NULL
	keepCols = c('perc_meth','coverage')
	# col_names = c('chr','start','perc_meth','coverage')
	# col_types = 'ci-di-'
	# stranded = FALSE
	classifier_flag = FALSE
	log10p_flag = FALSE
	resolution = 'base'
	prefix = gsub('_bt2.bismark.cov','', basename(file))
	# No categories
} else if (annot_type == 'simple_bisulfite_bismark') {
	# bis_align
	# BARPLOT + CATEGORICALS
	# Random for barplot and categoricals

	# head ./classifications/simple/IDH2mut_1_mc_hmc_bisulfite_bismark_simple_classification.bed
	# chr1	249239097	249239098	mc_hmc_high	1000	.	249239097	249239098	102,0,102
	# chr1	249239100	249239101	mc_hmc_high	1000	.	249239100	249239101	102,0,102
	# chr1	249239104	249239105	mc_hmc_high	1000	.	249239104	249239105	102,0,102
	# chr1	249239117	249239118	mc_hmc_high	1000	.	249239117	249239118	102,0,102
	extraCols = NULL
	rename_score = NULL
	rename_name = 'class'
	keepCols = c('class')
	# col_names = c('chr','start','end','class')
	# col_types = 'ciic-----'
	# stranded = FALSE
	classifier_flag = TRUE
	log10p_flag = FALSE
	resolution = 'base'
	prefix = gsub('.bed','', basename(file))

	# Make the correct category ordering
	if(mark == 'mc_hmc') {
		cats_order = c('mc_hmc_high', 'mc_hmc_med', 'mc_hmc_low', 'mc_hmc_none')
	} else if (mark == 'hmc') {
		cats_order = c('hmc_high', 'hmc_med', 'hmc_low', 'hmc_none')
	} else if (mark == 'mc') {
		cats_order = c('mc_high', 'mc_med', 'mc_low', 'mc_none')
	}

	# Variables for plot_categorical usage
	x_str = 'class'
	x_desc = 'Bisulfite Simple Classification'

	# Variable for summarize_categorical usage
	by_vec = c('class','annot.type')
} else if (annot_type == 'macs2') {
	# pulldown_sample
	# BARPLOT + NUMERICALS
	# Random for barplot

	# head ./pulldown/macs2_peaks/IDH2mut_2_hmc_pulldown_macs2_peaks.narrowPeak
	# chr21	9909565	9910066	IDH2mut_2_hmc_pulldown_macs2_peak_1	126	.	3.04226	16.46495	12.64028	213
	# chr21	15471800	15472211	IDH2mut_2_hmc_pulldown_macs2_peak_2	146	.	7.33008	18.53385	14.66105	186
	# chr21	15855465	15855886	IDH2mut_2_hmc_pulldown_macs2_peak_3	77	.	4.05129	11.47467	7.78926	178
	# chr21	15865175	15865647	IDH2mut_2_hmc_pulldown_macs2_peak_4	322	.	6.33401	36.40632	32.21310	240
	extraCols = c(fold = 'numeric', pval = 'numeric', qval = 'numeric', summit = 'numeric')
	rename_score = NULL
	rename_name = NULL
	keepCols = c('fold','pval')
	# col_names = c('chr','start','end','name','fold','pval')
	# col_types = 'ciic--dd--'
	# stranded = FALSE
	classifier_flag = FALSE
	log10p_flag = FALSE
	resolution = 'region'
	prefix = gsub('_peaks.narrowPeak','', basename(file))
	# No categories
} else if (annot_type == 'simple_pulldown_macs2') {
	# pulldown_sample
	# BARPLOT + CATEGORICALS
	# Random for barplot and categoricals

	# head ./classifications/simple/IDH2mut_1_hmc_pulldown_macs2_simple_classification.bed
	# chr21	9944160	9944341	hmc_low	1000	.	9944160	9944341	102,102,255
	# chr21	9999539	9999701	hmc_low	1000	.	9999539	9999701	102,102,255
	# chr21	10132638	10133209	hmc_low	1000	.	10132638	10133209	102,102,255
	# chr21	10135234	10135504	hmc_low	1000	.	10135234	10135504	102,102,255
	extraCols = NULL
	rename_name = 'class'
	rename_score = NULL
	keepCols = c('class')
	# col_names = c('chr','start','end','class')
	# col_types = 'ciic-----'
	# stranded = FALSE
	classifier_flag = TRUE
	log10p_flag = FALSE
	resolution = 'region'
	prefix = gsub('.bed','', basename(file))

	# Make the correct category ordering
	if (mark == 'hmc') {
		cats_order = c('hmc_high', 'hmc_med', 'hmc_low')
	} else if (mark == 'mc') {
		cats_order = c('mc_high', 'mc_med', 'mc_low')
	}

	# Variables for plot_categorical usage
	x_str = 'class'
	x_desc = 'macs2 Simple Classification'

	# Variable for summarize_categorical usage
	by_vec = c('class','annot.type')
} else if (annot_type == 'sample_class') {
	# sample_classification
	# BARPLOT + CATEGORICALS
	# Random for barplot and categoricals
	# NOTE: In the hybrid case there are so many regions that the All bar could
	# serve as the backgr

	# head ./classifications/sample/NBM_1_sample_classification.bed
	# chr1	10468	10472	unclassifiable	1000	.	10468	10472	192,192,192
	# chr1	10483	10485	unclassifiable	1000	.	10483	10485	192,192,192
	# chr1	10488	10490	unclassifiable	1000	.	10488	10490	192,192,192
	# chr1	10492	10494	unclassifiable	1000	.	10492	10494	192,192,192
	extraCols = NULL
	rename_name = 'class'
	rename_score = NULL
	keepCols = c('class')
	# col_names = c('chr','start','end','class')
	# col_types = 'ciic-----'
	# stranded = FALSE
	classifier_flag = TRUE
	log10p_flag = FALSE
	resolution = 'region'
	prefix = gsub('.bed','', basename(file))

	# The category ordering for this annot_type is computed below, after the file
	# is read because information from the file determines hybrid or pulldown.

	# Variables for plot_categorical usage
	x_str = 'class'
	x_desc = 'Sample Classification'

	# Variable for summarize_categorical usage
	by_vec = c('class','annot.type')
} else if (annot_type == 'dss') {
	# bis_compare
	# BARPLOT + NUMERICALS + CATEGORICALS
	# Random for barplot and categoricals
	# NOTE: Randoms will lack the CpG bias of dss regions...

	# head ./bisulfite/dss_calls/IDH2mut_v_NBM_mc_hmc_bisulfite_DMR_dss_for_annotatr.txt
	# chr21	9437472	9437521   IDH2mut 0.00022857   .   28.675   81.424	52.748
	# chr21	9437522	9437571   NBM	 0.02085531   .   22.738   85.875	63.136
	extraCols = c(meth_diff = 'numeric', mu1 = 'numeric', mu0 = 'numeric')
	rename_score = 'pval'
	rename_name = 'DM_status'
	keepCols = c('pval','DM_status','meth_diff','mu1','mu0')
	# col_names = c('chr','start','end','pval','DM_status','meth_diff','mu1','mu0')
	# col_types = 'ciidcddd'
	# stranded = FALSE
	classifier_flag = FALSE
	log10p_flag = TRUE
	resolution = 'region'
	prefix = gsub('_for_annotatr.txt','', basename(file))

	# Categories are the group1 and group0 command line parameters
	cats_order = c(group1, group0)

	# Variables for plot_categorical usage
	x_str = 'DM_status'
	x_desc = sprintf('Differential %s', mark)

	# Variable for summarize_categorical usage
	by_vec = c('DM_status','annot.type')
} else if (annot_type == 'csaw') {
	# pulldown_compare
	# BARPLOT + NUMERICALS + CATEGORICALS
	# Random for barplot and categoricals

	# head pulldown/csaw/IDH2mut_v_NBM_hmc_pulldown_csaw_for_annotatr.txt
	# chr21	46387680	46388640	IDH2mut	2.80016862356	*	2.28100636773e-165
	# chr21	36207120	36208800	IDH2mut	3.36810346155	*	2.74839325284e-143
	# chr21	39810080	39810880	IDH2mut	4.23550374893	*	1.25690543952e-142
	# chr21	46902720	46904640	IDH2mut	5.45917702717	*	1.27951073684e-106
	extraCols = c(pval = 'numeric')
	rename_score = 'fold'
	rename_name = 'group'
	keepCols = c('pval','group','fold')
	# col_names = c('chr','start','end','group','fold','pval')
	# col_types = 'ciicd-d'
	# stranded = FALSE
	classifier_flag = FALSE
	log10p_flag = TRUE
	resolution = 'region'
	prefix = gsub('_for_annotatr.txt','', basename(file))

	# Categories are the group1 and group0 command line parameters
	cats_order = c(group1, group0)

	# Variables for plot_categorical usage
	x_str = 'group'
	x_desc = sprintf('Differential %s', mark)

	# Variable for summarize_categorical usage
	by_vec = c('group','annot.type')
} else if (annot_type == 'simple_pulldown_csaw') {
	# pulldown_compare
	# BARPLOT + CATEGORICALS
	# Random for barplot and categoricals

	# head classifications/simple/HPVpos_v_HPVneg_hmc_pulldown_csaw_simple_classification.bed
	# chr1	565250	565350	diff_HPVpos_hmc_weak	1000	.	565250	565350	102,102,255
	# chr1	565700	565800	diff_HPVpos_hmc_weak	1000	.	565700	565800	102,102,255
	# chr1	565900	566000	diff_HPVneg_hmc_weak	1000	.	565900	566000	102,102,255
	# chr1	565950	566100	diff_HPVpos_hmc_strong	1000	.	565950	566100	0,0,102
	# chr1	567350	567450	diff_HPVneg_hmc_strong	1000	.	567350	567450	0,0,102
	extraCols = NULL
	rename_name = 'class'
	rename_score = NULL
	keepCols = c('class')
	# col_names = c('chr','start','end','class')
	# col_types = 'ciic-----'
	# stranded = FALSE
	classifier_flag = TRUE
	log10p_flag = FALSE
	resolution = 'region'
	prefix = gsub('.bed','', basename(file))

	# Make the correct category ordering
	if (mark == 'hmc') {
		cats_order = c(
			paste('diff', chip1, c('hmc_strong','hmc_mod','hmc_weak'), sep='_'),
			paste('diff', chip0, c('hmc_strong','hmc_mod','hmc_weak'), sep='_'))
	} else if (mark == 'mc') {
		cats_order = c(
			paste('diff', chip1, c('mc_strong','mc_mod','mc_weak'), sep='_'),
			paste('diff', chip0, c('mc_strong','mc_mod','mc_weak'), sep='_'))
	}

	# Variables for plot_categorical usage
	x_str = 'class'
	x_desc = 'csaw Simple Classification'

	# Variable for summarize_categorical usage
	by_vec = c('class','annot.type')
} else if (annot_type == 'compare_class') {
	# compare_classification
	# BARPLOT + CATEGORICALS
	# NOTE: There are so many regions that the All bar serves as the background

	# head ./classifications/comparison/IDH2mut_v_NBM_compare_classification.bed
	# chr1	0	10000	unclassifiable	1000	.	0	10000	192,192,192
	# chr1	10000	10162	no_DM	1000	.	10000	10162	0,0,0
	# chr1	10162	10185	unclassifiable	1000	.	10162	10185	192,192,192
	# chr1	10185	10276	no_DM	1000	.	10185	10276	0,0,0
	extraCols = NULL
	rename_name = 'class'
	rename_score = NULL
	keepCols = c('class')
	# col_names = c('chr','start','end','class')
	# col_types = 'ciic-----'
	# stranded = FALSE
	classifier_flag = TRUE
	log10p_flag = FALSE
	resolution = 'region'
	prefix = gsub('.bed','', basename(file))

	# Make the correct category ordering
	cats_order = c(
		'hyper_mc_hyper_hmc',
		'hypo_mc_hypo_hmc',
		'hyper_mc_hypo_hmc',
		'hypo_mc_hyper_hmc',
		'hyper_mc',
		'hyper_hmc',
		'hypo_mc',
		'hypo_hmc',
		'no_DM')

	# Variables for plot_categorical usage
	x_str = 'class'
	x_desc = 'Compare Classification'

	# Variable for summarize_categorical usage
	by_vec = c('class','annot.type')
} else {
	stop('annot_type is invalid. Must be one of bismark, simple_bisulfite_bismark, macs2, simple_pulldown_macs2, sample_class, dss, csaw, simple_pulldown_csaw, or compare_class.')
}

################################################################################
# Read the data and do some filtering / transformations as necessary
if(all(is.null(rename_name), is.null(rename_score))) {
	# Both rename_name and rename_score are NULL
	regions = read_regions(
    	con = file,
    	genome = genome,
        format = 'bed',
    	extraCols = extraCols)
} else if (all(!is.null(rename_name), !is.null(rename_score))) {
	# Neither rename_name nor rename_score are NULL
    regions = read_regions(
    	con = file,
    	genome = genome,
        format = 'bed',
    	extraCols = extraCols,
		rename_name = rename_name,
		rename_score = rename_score)
} else if (!is.null(rename_name)) {
	# Just rename_name is NULL
    regions = read_regions(
    	con = file,
    	genome = genome,
        format = 'bed',
    	extraCols = extraCols,
		rename_name = rename_name)
} else {
	# Just rename_score is NULL
	regions = read_regions(
    	con = file,
    	genome = genome,
        format = 'bed',
    	extraCols = extraCols,
		rename_score = rename_score)
}

# Use the keepCols
# NOTE: Due to some weird issue in GenomicRanges, doing this column subsetting
# causes the columns to be renamed, so we need to rename them. Fortunately,
# the subsetting will give the order of keepCols so we can just rename them.
mcols(regions) = mcols(regions)[, keepCols]
colnames(mcols(regions)) = keepCols

# Fix resolution='base' to make start = end
# This is a new and fun addition to rtracklayer::import()
if(resolution == 'base') {
	start(regions) = end(regions)
}

# Filter out unclassifiable classification
if(classifier_flag) {
	regions = regions[regions$class != 'unclassifiable']
}

# Convert pval column to -log10(pval)
if(log10p_flag) {
	regions$pval = -log10(regions$pval)
}

################################################################################
# Deal with random regions depending on resolution
# NOTE: Do not do random regions for dss at the moment, uses too much RAM
if(resolution == 'region' && annot_type != 'sample_class' && annot_type != 'compare_class' && annot_type != 'dss') {
	regions_rnd = NULL
} else if (resolution == 'base' && genome == 'hg19') {
	regions_rnd = NULL
} else {
	regions_rnd = NULL
}

################################################################################
# Pick annotations and order them
# Created variables: annot_cpg_order, annot_gene_order, annot_all_order
if(genome %in% c('hg19','hg38','mm9','mm10')) {
	annot_cpg_order = paste(genome, c(
		'cpg_islands',
		'cpg_shores',
		'cpg_shelves',
		'cpg_inter'), sep='_')
	annot_gene_order = paste(genome, c(
		'genes_1to5kb',
		'genes_promoters',
		'genes_5UTRs',
		'genes_exons',
		'genes_introns',
		'genes_exonintronboundaries',
		'genes_intronexonboundaries',
		'genes_3UTRs',
		'genes_intergenic'), sep='_')
	annot_all_order = c(
		annot_cpg_order,
		annot_gene_order)

	if(genome == 'hg19' || genome == 'mm9') {
		annot_gene_order = c(annot_gene_order, paste(genome, 'enhancers_fantom', sep='_'))
		annot_all_order = c(annot_all_order, paste(genome, 'enhancers_fantom', sep='_'))
	}

	if(genome == 'hg19' || genome == 'hg38' || genome == 'mm10') {
		annot_gene_order = c(annot_gene_order, paste(genome, 'lncrna_gencode', sep='_'))
		annot_all_order = c(annot_all_order, paste(genome, 'lncrna_gencode', sep='_'))
	}
} else {
	# TODO: Add support for insertion of custom annotation path via command line
	stop('Only hg19, hg38, mm9, and mm10 annotations are supported for this part of mint at the moment.')
}

################################################################################
# Build the annotations
annotations = build_annotations(genome = genome, annotations = annot_all_order)

################################################################################
# Order categories for sample classification case (others are above)
# Must have read in the file to deduce the classes
# Created variables: cats_order
if (annot_type == 'sample_class') {
	cats = unique(regions$class)
	if('mc_or_hmc' %in% cats) {
		# hybrid sample
		cats_order = c(
			'mc',
			'mc_low',
			'hmc',
			'mc_or_hmc',
			'mc_or_hmc_low',
			'no_meth')
	} else if ('mc_and_hmc' %in% cats) {
		# pulldown sample
		cats_order = c(
			'mc_and_hmc',
			'mc',
			'hmc',
			'no_meth')
	}
}

################################################################################
# For strict macs2 and csaw output, create dplyr::tbl_df on which to create
# overall distributions of fold change
if(annot_type == 'macs2') {
	regions_tbl = dplyr::tbl_df(data.frame(
		'seqnames' = seqnames(regions),
		'start' = start(regions),
		'end' = end(regions),
		'fold' = mcols(regions)$fold,
		'pval' = mcols(regions)$pval,
		stringsAsFactors=F))
} else if (annot_type == 'csaw') {
	regions_tbl = dplyr::tbl_df(data.frame(
		'seqnames' = seqnames(regions),
		'start' = start(regions),
		'end' = end(regions),
		'group' = mcols(regions)$group,
		'fold' = mcols(regions)$fold,
		'pval' = mcols(regions)$pval,
		stringsAsFactors=F))
} else {
	regions_tbl = NULL
}

################################################################################
# Do the annotation of the regions

annotated_regions = annotate_regions(regions = regions, annotations = annotations)

# If regions_rnd has GRanges class, then annotate, otherwise you read in the CpG annotation table
if(!is.null(regions_rnd)) {
	annotated_regions_rnd = annotate_regions(regions = regions_rnd, annotations = annotations)
} else {
	annotated_regions_rnd = NULL
}

################################################################################
# Summarize annotations
if(!is.null(annotated_regions_rnd)) {
	# Use randomized regions
	count_annots = summarize_annotations(annotated_regions = annotated_regions, annotated_random = annotated_regions_rnd)
} else {
	# Do not use randomized regions
	count_annots = summarize_annotations(annotated_regions = annotated_regions)
}

# Write it
readr::write_tsv(x = count_annots, path = sprintf('summary/tables/%s_annotation_counts.txt', prefix))

################################################################################
# Plots

#############################################################
# Barplot of # of regions per annotation
# NOTE: This is the same for all the types
if(!is.null(annotated_regions_rnd)) {
	# Use randomized regions
	file_png = sprintf('summary/figures/%s_counts.png', prefix)
	plot_counts = plot_annotation(
		annotated_regions = annotated_regions,
		annotated_random = annotated_regions_rnd,
		annotation_order = annot_all_order,
		plot_title = sprintf('%s regions per annotation', prefix),
		x_label = 'Annotations',
		y_label = '# Regions')
	ggsave(filename = file_png, plot = plot_counts, width = 8, height = 6)
} else {
	# Use CpG annotations as a post annotatr add on
	file_png = sprintf('summary/figures/%s_counts.png', prefix)
	plot_counts = plot_annotation(
		annotated_regions = annotated_regions,
		annotation_order = annot_all_order,
		plot_title = sprintf('%s regions per annotation', prefix),
		x_label = 'Annotations',
		y_label = '# Regions')
	ggsave(filename = file_png, plot = plot_counts, width = 8, height = 6)
}

#############################################################
# plot_numerical specific to bismark
if(annot_type == 'bismark') {
	##############################
	# CpG read coverage (limited to 500)
	# NOTE: Will receive a warning like:
	# Warning messages:
	# 1: Removed 578 rows containing non-finite values (stat_bin).
	# 2: Removed 6936 rows containing non-finite values (stat_bin).
	# When there are coverages exceeding the x limits
	file_png = sprintf('summary/figures/%s_coverage.png', prefix)
	plot_coverage = plot_numerical(
		annotated_regions = annotated_regions,
		x = 'coverage',
		facet = 'annot.type',
		facet_order = annot_all_order,
		bin_width = 10,
		plot_title = sprintf('%s coverage over annotations', prefix),
		x_label = 'Coverage',
		legend_facet_label = 'Coverage in annotation',
		legend_cum_label = 'Overall coverage')
	plot_coverage = plot_coverage + xlim(0,500)
	ggsave(filename = file_png, plot = plot_coverage, width = 8, height = 8)

	##############################
	# Percent methylation over annotations
	file_png = sprintf('summary/figures/%s_percmeth.png', prefix)
	plot_percmeth = plot_numerical(
		annotated_regions = annotated_regions,
		x = 'perc_meth',
		facet = 'annot.type',
		facet_order = annot_all_order,
		bin_width = 5,
		plot_title = sprintf('%s perc. meth. over annotations', prefix),
		x_label = 'Percent Methylation',
		legend_facet_label = 'Percent methylation in annotation',
		legend_cum_label = 'Overall percent methylation')
	ggsave(filename = file_png, plot = plot_percmeth, width = 8, height = 8)
}

#############################################################
# plot_numerical specific to dss
if(annot_type == 'dss') {
	##############################
	# Annotation counts by categorical
	sum_cat = summarize_categorical(annotated_regions = annotated_regions, by = by_vec)
	readr::write_tsv(x = sum_cat, path = sprintf('summary/tables/%s_annotation_counts_by_category.txt', prefix))

	##############################
	# Methylation difference over annotations
	file_png = sprintf('summary/figures/%s_methdiff.png', prefix)
	plot_methdiff = plot_numerical(
		annotated_regions = annotated_regions,
		x = 'meth_diff',
		facet = 'annot.type',
		facet_order = annot_all_order,
		bin_width = 5,
		plot_title = sprintf('%s meth. diff. over annotations', prefix),
		x_label = sprintf('Methylation Difference (%s - %s)', group1, group0),
		legend_facet_label = 'Methylation difference in annotation',
		legend_cum_label = 'Overall methylation difference')
	ggplot2::ggsave(filename = file_png, plot = plot_methdiff, width = 8, height = 8)

	##############################
	# Volcano plots over annotations
	file_png = sprintf('summary/figures/%s_volcano.png', prefix)
	plot_volcano = plot_numerical(
		annotated_regions = annotated_regions,
		x = 'meth_diff',
		y = 'pval',
		facet = 'annot.type',
		facet_order = annot_all_order,
		plot_title = sprintf('%s meth. diff. vs -log10(pval)', prefix),
		x_label = sprintf('Methylation Difference (%s - %s)', group1, group0),
		y_label = '-log10(pval)')
	ggplot2::ggsave(filename = file_png, plot = plot_volcano, width = 8, height = 8)
}

#############################################################
# Overall fold-change, peak width, and volcano plots for macs2 & csaw

if(!is.null(regions_tbl)) {
	##############################
	# Peak widths
	widths = data.frame(
		region = 1:length(regions),
		width = end(regions) - start(regions), stringsAsFactors=F)

	file_png = sprintf('summary/figures/%s_peak_widths.png', prefix)
	plot_region_widths = ggplot(data = widths, aes(x = width, y=..density..)) +
		scale_x_log10() + geom_histogram(stat = 'bin', fill = NA, color='gray', bins=30) +
		theme_bw() +
		xlab('Peak Widths (log10 scale)') +
		ggtitle(sprintf('%s peak widths', prefix))
	ggsave(filename = file_png, plot = plot_region_widths, width = 6, height = 6)

	##############################
	# Fold changes
	file_png = sprintf('summary/figures/%s_foldchg_overall.png', prefix)
	plot_foldchg = ggplot(data = regions_tbl, aes(x = fold, y=..density..)) +
		geom_histogram(stat = 'bin', fill = NA, color='gray', bins=30) +
		theme_bw() +
		xlab('Fold Change') +
		ggtitle(sprintf('%s Overall Fold Change', prefix))
	ggsave(filename = file_png, plot = plot_foldchg, width = 6, height = 6)

	##############################
	# Volcano plots
	file_png = sprintf('summary/figures/%s_volcano_overall.png', prefix)
	plot_volcano = ggplot(data = regions_tbl, aes(x = fold, y = pval)) +
		geom_point(alpha = 1/8, size = 1) +
		theme_bw() +
		xlab('Fold Change') +
		ylab('-log10(pval)')
		ggtitle(sprintf('%s Overall Fold Change', prefix))
	ggsave(filename = file_png, plot = plot_volcano, width = 6, height = 6)
}

#############################################################
# plot_numerical specific to macs2
if(annot_type == 'macs2') {
	##############################
	# Fold change over annotations
	file_png = sprintf('summary/figures/%s_foldchg_annots.png', prefix)
	plot_foldchg = plot_numerical(
		annotated_regions = annotated_regions,
		x = 'fold',
		facet = 'annot.type',
		facet_order = annot_all_order,
		bin_width = 5,
		plot_title = sprintf('%s fold change over annotations', prefix),
		x_label = 'Fold Change',
		legend_facet_label = 'Fold change in annotation',
		legend_cum_label = 'Overall fold change')
	ggsave(filename = file_png, plot = plot_foldchg, width = 8, height = 8)

	##############################
	# Volcano plots over annotations
	file_png = sprintf('summary/figures/%s_volcano_annots.png', prefix)
	plot_volcano = plot_numerical(
		annotated_regions = annotated_regions,
		x = 'fold',
		y = 'pval',
		facet = 'annot.type',
		facet_order = annot_all_order,
		plot_title = sprintf('%s fold change vs -log10(pval)', prefix),
		x_label = 'Fold change IP / input',
		y_label = '-log10(pval)')
	ggsave(filename = file_png, plot = plot_volcano, width = 8, height = 8)
}

#############################################################
# plot_numerical specific to csaw
if(annot_type == 'csaw') {
	# ##############################
	# # Fold change with facet over group (chip1/chip0)
	# file_png = sprintf('summary/figures/%s_foldchg_groups.png', prefix)
	# plot_foldchg = plot_numerical(
	# 	annotated_regions = regions_tbl,
	# 	x = 'fold',
	# 	facet = 'group',
	# 	facet_order = cats_order,
	# 	bin_width = 5,
	# 	plot_title = sprintf('%s fold change over annotations', prefix),
	# 	x_label = 'Fold Change')
	# ggsave(filename = file_png, plot = plot_foldchg, width = 12, height = 6)
	#
	# ##############################
	# # Volcano plots with facet over group (chip1/chip0)
	# file_png = sprintf('summary/figures/%s_volcano_groups.png', prefix)
	# plot_volcano = plot_numerical(
	# 	annotated_regions = regions_tbl,
	# 	x = 'fold',
	# 	y = 'pval',
	# 	facet = 'group',
	# 	facet_order = cats_order,
	# 	plot_title = sprintf('%s fold change vs -log10(pval)', prefix),
	# 	x_label = 'Fold change',
	# 	y_label = '-log10(pval)')
	# ggsave(filename = file_png, plot = plot_volcano, width = 12, height = 6)

	##############################
	# Fold change in chip1 with facet over annots
	file_png = sprintf('summary/figures/%s_foldchg_%s_annots.png', prefix, chip1)
	plot_foldchg = plot_numerical(
		annotated_regions = subset(annotated_regions, group == chip1),
		x = 'fold',
		facet = 'annot.type',
		facet_order = annot_all_order,
		bin_width = 5,
		plot_title = sprintf('%s %s fold change over annotations', prefix, chip1),
		x_label = sprintf('%s fold change', chip1),
		legend_facet_label = 'Fold change in annotation',
		legend_cum_label = 'Overall fold change')
	ggsave(filename = file_png, plot = plot_foldchg, width = 8, height = 8)

	##############################
	# Volcano in chip1 with facet over annots
	file_png = sprintf('summary/figures/%s_volcano_%s_annots.png', prefix, chip1)
	plot_volcano = plot_numerical(
		annotated_regions = subset(annotated_regions, group == chip1),
		x = 'fold',
		y = 'pval',
		facet = 'annot.type',
		facet_order = annot_all_order,
		plot_title = sprintf('%s %s fold change vs -log10(pval)', prefix, chip1),
		x_label = sprintf('%s fold change', chip1),
		y_label = '-log10(pval)')
	ggsave(filename = file_png, plot = plot_volcano, width = 8, height = 8)

	##############################
	# Fold change in chip0 with facet over annots
	file_png = sprintf('summary/figures/%s_foldchg_%s_annots.png', prefix, chip0)
	plot_foldchg = plot_numerical(
		annotated_regions = subset(annotated_regions, group == chip0),
		x = 'fold',
		facet = 'annot.type',
		facet_order = annot_all_order,
		bin_width = 5,
		plot_title = sprintf('%s %s fold change over annotations', prefix, chip0),
		x_label = sprintf('%s fold change', chip0),
		legend_facet_label = 'Fold change in annotation',
		legend_cum_label = 'Overall fold change')
	ggsave(filename = file_png, plot = plot_foldchg, width = 8, height = 8)

	##############################
	# Volcano in chip0 with facet over annots
	file_png = sprintf('summary/figures/%s_volcano_%s_annots.png', prefix, chip0)
	plot_volcano = plot_numerical(
		annotated_regions = subset(annotated_regions, group == chip0),
		x = 'fold',
		y = 'pval',
		facet = 'annot.type',
		facet_order = annot_all_order,
		plot_title = sprintf('%s %s fold change vs -log10(pval)', prefix, chip0),
		x_label = sprintf('%s fold change', chip0),
		y_label = '-log10(pval)')
	ggsave(filename = file_png, plot = plot_volcano, width = 8, height = 8)
}

#############################################################
# plot_categorical specific to csaw, dss, simple, sample, and compare classifications
# Use x_str and x_desc variables determined near the top of this file when determining annot_type
if(annot_type != 'bismark' && annot_type != 'macs2') {
	##############################
	# Annotation counts by categorical
	sum_cat = summarize_categorical(annotated_regions = annotated_regions, by = by_vec)
	readr::write_tsv(x = sum_cat, path = sprintf('summary/tables/%s_annotation_counts_by_category.txt', prefix))

	if(!is.null(annotated_regions_rnd)) {
		##############################
		# Regions split by category and stacked by CpG annotations (count)
		file_png = sprintf('summary/figures/%s_cat_count_cpgs.png', prefix)
		plot_cat_count_cpgs = plot_categorical(
			annotated_regions = annotated_regions,
			annotated_random = annotated_regions_rnd,
			x=x_str, fill='annot.type',
			x_order = cats_order, fill_order = annot_cpg_order, position='stack',
			plot_title = prefix,
			legend_title = 'CpG annots.',
			x_label = x_desc,
			y_label = 'Count')
		ggsave(filename = file_png, plot = plot_cat_count_cpgs, width = 8, height = 8)

		##############################
		# Regions split by category and stacked by knownGene annotations (count)
		file_png = sprintf('summary/figures/%s_cat_count_genes.png', prefix)
		plot_cat_count_genes = plot_categorical(
			annotated_regions = annotated_regions,
			annotated_random = annotated_regions_rnd,
			x=x_str, fill='annot.type',
			x_order = cats_order, fill_order = annot_gene_order, position='stack',
			plot_title = prefix,
			legend_title = 'knownGene annots.',
			x_label = x_desc,
			y_label = 'Count')
		ggsave(filename = file_png, plot = plot_cat_count_genes, width = 8, height = 8)

		##############################
		# Regions split by category and filled by CpG annotations (prop)
		file_png = sprintf('summary/figures/%s_cat_prop_cpgs.png', prefix)
		plot_cat_prop_cpgs = plot_categorical(
			annotated_regions = annotated_regions,
			annotated_random = annotated_regions_rnd,
			x=x_str, fill='annot.type',
			x_order = cats_order, fill_order = annot_cpg_order, position='fill',
			plot_title = prefix,
			legend_title = 'CpG annots.',
			x_label = x_desc,
			y_label = 'Proportion')
		ggsave(filename = file_png, plot = plot_cat_prop_cpgs, width = 8, height = 8)

		##############################
		# Regions split by category and filled by knownGene annotations (prop)
		file_png = sprintf('summary/figures/%s_cat_prop_genes.png', prefix)
		plot_cat_prop_genes = plot_categorical(
			annotated_regions = annotated_regions,
			annotated_random = annotated_regions_rnd,
			x=x_str, fill='annot.type',
			x_order = cats_order, fill_order = annot_gene_order, position='fill',
			plot_title = prefix,
			legend_title = 'knownGene annots.',
			x_label = x_desc,
			y_label = 'Proportion')
		ggsave(filename = file_png, plot = plot_cat_prop_genes, width = 8, height = 8)
	} else {
		##############################
		# Regions split by category and stacked by CpG annotations (count)
		file_png = sprintf('summary/figures/%s_cat_count_cpgs.png', prefix)
		plot_cat_count_cpgs = plot_categorical(
			annotated_regions = annotated_regions,
			x=x_str, fill='annot.type',
			x_order = cats_order, fill_order = annot_cpg_order, position='stack',
			plot_title = prefix,
			legend_title = 'CpG annots.',
			x_label = x_desc,
			y_label = 'Count')
		ggsave(filename = file_png, plot = plot_cat_count_cpgs, width = 8, height = 8)

		##############################
		# Regions split by category and stacked by knownGene annotations (count)
		file_png = sprintf('summary/figures/%s_cat_count_genes.png', prefix)
		plot_cat_count_genes = plot_categorical(
			annotated_regions = annotated_regions,
			x=x_str, fill='annot.type',
			x_order = cats_order, fill_order = annot_gene_order, position='stack',
			plot_title = prefix,
			legend_title = 'knownGene annots.',
			x_label = x_desc,
			y_label = 'Count')
		ggsave(filename = file_png, plot = plot_cat_count_genes, width = 8, height = 8)

		##############################
		# Regions split by category and filled by CpG annotations (prop)
		file_png = sprintf('summary/figures/%s_cat_prop_cpgs.png', prefix)
		plot_cat_prop_cpgs = plot_categorical(
			annotated_regions = annotated_regions,
			x=x_str, fill='annot.type',
			x_order = cats_order, fill_order = annot_cpg_order, position='fill',
			plot_title = prefix,
			legend_title = 'CpG annots.',
			x_label = x_desc,
			y_label = 'Proportion')
		ggsave(filename = file_png, plot = plot_cat_prop_cpgs, width = 8, height = 8)

		##############################
		# Regions split by category and filled by knownGene annotations (prop)
		file_png = sprintf('summary/figures/%s_cat_prop_genes.png', prefix)
		plot_cat_prop_genes = plot_categorical(
			annotated_regions = annotated_regions,
			x=x_str, fill='annot.type',
			x_order = cats_order, fill_order = annot_gene_order, position='fill',
			plot_title = prefix,
			legend_title = 'knownGene annots.',
			x_label = x_desc,
			y_label = 'Proportion')
		ggsave(filename = file_png, plot = plot_cat_prop_genes, width = 8, height = 8)
	}
}

################################################################################
# Record .RData object of the session
save.image(file = sprintf('RData/%s_annotatr_analysis.RData', prefix))
