library(annotatr)
library(readr)
library(ggplot2)
library(optparse)
library(GenomicRanges)

option_list = list(
  make_option('--file', type='character'),
  make_option('--genome', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

file = opt$file
genome = opt$genome

prefix = gsub('_for_annotatr.txt','', basename(file))
if(grepl('simple', prefix)) {
	class_type = 'simple'
	display_type = 'Simple'
} else if (grepl('sample', prefix)) {
	class_type = 'sample'
	display_type = 'Sample'
} else if (grepl('compare', prefix)) {
	class_type = 'compare'
	display_type = 'Compare'
} else if (grepl('PePr', prefix)) {
	class_type = 'PePr'
	display_type = 'PePr'
} else if (grepl('macs2', prefix)) {
	class_type = 'macs2'
	display_type = 'macs2'
} else {
	stop('Invalid classification type. Must be one of simple, sample, or compare.')
}

###############################################################
# Read
if(class_type == 'macs2' || class_type == 'PePr') {
	r = read_bed(
		file = file,
		col.names = c('chr','start','end','name','fold','strand','pvalue'),
		genome = genome,
		stranded = FALSE,
		use.score = TRUE)

	# Remove unclassifiable regions
	r = r[r$name != 'unclassifiable']
} else {
	r = read_bed(
		file = file,
		col.names = FALSE,
		genome = genome,
		stranded = FALSE,
		use.score = FALSE)

	# Remove unclassifiable regions
	r = r[r$name != 'unclassifiable']
}

###############################################################
# Make random regions

r_rand = randomize_regions(
	regions = r,
	genome = 'hg19',)

###############################################################
# Pick annotations
if(genome %in% c('hg19','hg38','mm9','mm10')) {
	a = c(
		sprintf('%s_cpgs', genome),
		sprintf('%s_basicgenes', genome),
		sprintf('%s_knownGenes_intergenic', genome))
} else {
	# TODO: Add support for insertion of custom annotation path via command line
	stop('Only hg19, hg38, mm9, and mm10 annotations are supported for this part of mint at the moment.')
}

if(genome == 'hg19') {
	a = c(a, 'hg19_enhancers_fantom')
}

###############################################################
# Order annotations and categories
if(genome %in% c('hg19','hg38','mm9','mm10')) {
	if(genome == 'hg19') {
		a_cpg_order = paste(genome, c(
			'cpg_islands',
			'cpg_shores',
			'cpg_shelves',
			'cpg_inter'), sep='_')
		a_gene_order = paste(genome, c(
			'enhancers_fantom',
			'knownGenes_1to5kb',
			'knownGenes_promoters',
			'knownGenes_5UTRs',
			'knownGenes_exons',
			'knownGenes_introns',
			'knownGenes_3UTRs',
			'knownGenes_intergenic'), sep='_')
		a_all_order = c(
			a_cpg_order,
			a_gene_order)
	} else {
		a_cpg_order = paste(genome, c(
			'cpg_islands',
			'cpg_shores',
			'cpg_shelves',
			'cpg_inter'), sep='_')
		a_gene_order = paste(genome, c(
			'knownGenes_1to5kb',
			'knownGenes_promoters',
			'knownGenes_5UTRs',
			'knownGenes_exons',
			'knownGenes_introns',
			'knownGenes_3UTRs',
			'knownGenes_intergenic'), sep='_')
		a_all_order = c(
			a_cpg_order,
			a_gene_order)
	}
}

if(class_type == 'simple') {
	cats = unique(r$name)
	cats_order = c(
		grep('high', cats, value=TRUE),
		grep('med', cats, value=TRUE),
		grep('low', cats, value=TRUE),
		grep('none', cats, value=TRUE))
} else if (class_type == 'sample') {
	cats = unique(r$name)
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
} else if (class_type == 'compare') {
	cats_order = c(
		'hyper_mc_hyper_hmc',
		'hypo_mc_hypo_hmc',
		'hyper_mc_hypo_hmc',
		'hypo_mc_hyper_hmc',
		'hyper_mc',
		'hyper_hmc',
		'hypo-mc',
		'hypo_hmc',
		'no_DM')
} else if (class_type == 'PePr') {
	cats_order = c(
		'chip1',
		'chip2')
}

###############################################################
# Annotate regions
if(class_type == 'macs2' || class_type == 'PePr') {
	ar = annotate_regions(
		regions = r,
		annotations = a,
		ignore.strand = TRUE,
		use.score = TRUE)

	# Get just the plain tbl with no annotations for overall distributions
	r_tbl = dplyr::tbl_df(data.frame(
		'chrom' = seqnames(r),
		'start' = start(r),
		'end' = end(r),
		'name' = mcols(r)$name,
		'fold' = mcols(r)$fold,
		'strand' = strand(r),
		'pvalue' = mcols(r)$pvalue,
		stringsAsFactors=F))
} else {
	ar = annotate_regions(
		regions = r,
		annotations = a,
		ignore.strand = TRUE,
		use.score = FALSE)
}

# Write it
readr::write_tsv(x = ar, path = sprintf('summary/tables/%s_annotations.txt', prefix))

###############################################################
# Annotate random regions

ar_rand = annotate_regions(
	regions = r_rand,
	annotations = a,
	ignore.strand = TRUE,
	use.score= FALSE)

###############################################################
# Summarize annotations
count_annots = summarize_annotations(annotated_regions = ar, annotated_random = ar_rand)

# Write it
readr::write_tsv(x = count_annots, path = sprintf('summary/tables/%s_annotation_counts.txt', prefix))

###############################################################
# Visualizations

# Regular barplot of regions in annotations
counts_png = sprintf('summary/figures/%s_counts.png', prefix)
plot_counts = plot_annotation(
	annotated_regions = ar,
	annotated_random = ar_rand,
	annotation_order = a_all_order,
	plot_title = sprintf('%s regions per annotation', prefix),
	x_label = 'Annotations',
	y_label = '# Regions')
ggplot2::ggsave(filename = counts_png, plot = plot_counts, width = 8, height = 8)

# Other classifications are too big for this
# if(class_type == 'simple' || class_type == 'PePr') {
# 	# Heatmap of regions in pairs of annotations
# 	cocounts_png = sprintf('summary/figures/%s_cocounts.png', prefix)
# 	plot_cocounts = plot_coannotations(
# 		annotated_regions = ar,
# 		annotation_order = a_all_order,
# 		plot_title = sprintf('%s regions in pairs of annotations', prefix),
# 		axes_label = 'Annotations')
# 	ggplot2::ggsave(filename = cocounts_png, plot = plot_cocounts, width = 8, height = 8)
# }

# Histogram of region (peak) widths for *pulldown_simple* and PePr inputs
if( (class_type == 'simple' && grepl('pulldown', prefix)) || class_type == 'PePr' || class_type == 'macs2' ) {
	widths = data.frame(
		region = 1:length(r),
		width = end(r) - start(r), stringsAsFactors=F)

	widths_png = sprintf('summary/figures/%s_peak_widths.png', prefix)
	plot_region_widths = ggplot(data = widths, aes(x = width, y=..density..)) +
		scale_x_log10() + geom_histogram(stat = 'bin', fill = NA, color='gray', bins=30) +
		theme_bw() +
		xlab('Peak Widths (log10 scale)') +
		ggtitle(sprintf('%s peak widths', prefix))
	ggplot2::ggsave(filename = widths_png, plot = plot_region_widths, width = 6, height = 6)
}

# Fold change and volcano plots for macs2 narrowPeak
if( class_type == 'macs2') {
	foldchg_png = sprintf('summary/figures/%s_foldchg_overall.png', prefix)
	plot_foldchg = ggplot(data = r_tbl, aes(x = fold, y=..density..)) +
		geom_histogram(stat = 'bin', fill = NA, color='gray', bins=30) +
		theme_bw() +
		xlab('Fold Change') +
		ggtitle(sprintf('%s Overall Fold Change', prefix))
	ggplot2::ggsave(filename = foldchg_png, plot = plot_foldchg, width = 6, height = 6)

	volcano_png = sprintf('summary/figures/%s_volcano_overall.png', prefix)
	plot_volcano = ggplot(data = r_tbl, aes(x = fold, y = pvalue)) +
		geom_point(alpha = 1/8, size = 1) +
		theme_bw() +
		xlab('Fold Change IP / Input') +
		ylab('-log10(pval)')
		ggtitle(sprintf('%s Overall Fold Change', prefix))
	ggplot2::ggsave(filename = volcano_png, plot = plot_volcano, width = 6, height = 6)

	foldchg_png = sprintf('summary/figures/%s_foldchg_annots.png', prefix)
	plot_foldchg = plot_numerical(
		tbl = ar,
		x = 'fold',
		facet = 'annot_type',
		facet_order = a_all_order,
		bin_width = 5,
		plot_title = sprintf('%s fold change over annotations', prefix),
		x_label = 'Fold Change')
	ggplot2::ggsave(filename = foldchg_png, plot = plot_foldchg, width = 8, height = 8)

	volcano_png = sprintf('summary/figures/%s_volcano_annots.png', prefix)
	plot_volcano = plot_numerical(
		tbl = ar,
		x = 'fold',
		y = 'pvalue',
		facet = 'annot_type',
		facet_order = a_all_order,
		plot_title = sprintf('%s fold change vs -log10(pval)', prefix),
		x_label = 'Fold change IP / input',
		y_label = '-log10(pval)')
	ggplot2::ggsave(filename = volcano_png, plot = plot_volcano, width = 8, height = 8)
}

# Fold change and volcano plots for PePr output
if(class_type == 'PePr') {
	###############################################################
	# Fold change and volcano with facet over name (chip1/chip2)
	foldchg_png = sprintf('summary/figures/%s_foldchg_overall.png', prefix)
	plot_foldchg = plot_numerical(
		tbl = r_tbl,
		x = 'fold',
		facet = 'name',
		facet_order = c('chip1','chip2'),
		bin_width = 5,
		plot_title = sprintf('%s fold change over annotations', prefix),
		x_label = 'Fold Change')
	ggplot2::ggsave(filename = foldchg_png, plot = plot_foldchg, width = 8, height = 8)

	volcano_png = sprintf('summary/figures/%s_volcano_overall.png', prefix)
	plot_volcano = plot_numerical(
		tbl = r_tbl,
		x = 'fold',
		y = 'pvalue',
		facet = 'name',
		facet_order = c('chip1','chip2'),
		plot_title = sprintf('%s fold change vs -log10(pval)', prefix),
		x_label = 'Fold change',
		y_label = '-log10(pval)')
	ggplot2::ggsave(filename = volcano_png, plot = plot_volcano, width = 8, height = 8)

	###############################################################
	# Fold change and volcano in chip1 with facet over annots
	foldchg_png = sprintf('summary/figures/%s_foldchg_chip1_annots.png', prefix)
	plot_foldchg = plot_numerical(
		tbl = subset(ar, name == 'chip1'),
		x = 'fold',
		facet = 'annot_type',
		facet_order = a_all_order,
		bin_width = 5,
		plot_title = sprintf('%s chip1 fold change over annotations', prefix),
		x_label = 'chip1 fold Change')
	ggplot2::ggsave(filename = foldchg_png, plot = plot_foldchg, width = 8, height = 8)

	volcano_png = sprintf('summary/figures/%s_volcano_chip1_annots.png', prefix)
	plot_volcano = plot_numerical(
		tbl = subset(ar, name == 'chip1'),
		x = 'fold',
		y = 'pvalue',
		facet = 'annot_type',
		facet_order = a_all_order,
		plot_title = sprintf('%s chip1 fold change vs -log10(pval)', prefix),
		x_label = 'chip1 fold change',
		y_label = '-log10(pval)')
	ggplot2::ggsave(filename = volcano_png, plot = plot_volcano, width = 8, height = 8)

	###############################################################
	# Fold change and volcano in chip2 with facet over annots
	foldchg_png = sprintf('summary/figures/%s_foldchg_chip2_annots.png', prefix)
	plot_foldchg = plot_numerical(
		tbl = subset(ar, name == 'chip2'),
		x = 'fold',
		facet = 'annot_type',
		facet_order = a_all_order,
		bin_width = 5,
		plot_title = sprintf('%s chip2 fold change over annotations', prefix),
		x_label = 'Fold Change')
	ggplot2::ggsave(filename = foldchg_png, plot = plot_foldchg, width = 8, height = 8)

	volcano_png = sprintf('summary/figures/%s_volcano_chip2_annots.png', prefix)
	plot_volcano = plot_numerical(
		tbl = subset(ar, name == 'chip2'),
		x = 'fold',
		y = 'pvalue',
		facet = 'annot_type',
		facet_order = a_all_order,
		plot_title = sprintf('%s chip2 fold change vs -log10(pval)', prefix),
		x_label = 'chip2 fold change',
		y_label = '-log10(pval)')
	ggplot2::ggsave(filename = volcano_png, plot = plot_volcano, width = 8, height = 8)
}

if( class_type != 'macs2') {
	# Regions split by category and stacked by CpG annotations (count)
	cat_count_cpgs_png = sprintf('summary/figures/%s_cat_count_cpgs.png', prefix)
	plot_cat_count_cpgs = plot_categorical(
	  annotated_regions = ar,
	  annotated_random = ar_rand,
	  x='name', fill='annot_type',
	  x_order = cats_order, fill_order = a_cpg_order, position='stack',
	  plot_title = sprintf('%s classification by Annotation', display_type),
	  legend_title = 'Annotations',
	  x_label = sprintf('%s classification', display_type),
	  y_label = 'Count')
	ggplot2::ggsave(filename = cat_count_cpgs_png, plot = plot_cat_count_cpgs, width = 8, height = 8)

	# Regions split by category and stacked by knownGene annotations (count)
	cat_count_genes_png = sprintf('summary/figures/%s_cat_count_genes.png', prefix)
	plot_cat_count_genes = plot_categorical(
	  annotated_regions = ar,
	  annotated_random = ar_rand,
	  x='name', fill='annot_type',
	  x_order = cats_order, fill_order = a_gene_order, position='stack',
	  plot_title = sprintf('%s classification by Annotation', display_type),
	  legend_title = 'Annotations',
	  x_label = sprintf('%s classification', display_type),
	  y_label = 'Count')
	ggplot2::ggsave(filename = cat_count_genes_png, plot = plot_cat_count_genes, width = 8, height = 8)

	# Regions split by category and filled by CpG annotations (prop)
	cat_prop_cpgs_png = sprintf('summary/figures/%s_cat_prop_cpgs.png', prefix)
	plot_cat_prop_cpgs = plot_categorical(
	  annotated_regions = ar,
	  annotated_random = ar_rand,
	  x='name', fill='annot_type',
	  x_order = cats_order, fill_order = a_cpg_order, position='fill',
	  plot_title = sprintf('%s classification by Annotation', display_type),
	  legend_title = 'Annotations',
	  x_label = sprintf('%s classification', display_type),
	  y_label = 'Proportion')
	ggplot2::ggsave(filename = cat_prop_cpgs_png, plot = plot_cat_prop_cpgs, width = 8, height = 8)

	# Regions split by category and filled by knownGene annotations (prop)
	cat_prop_genes_png = sprintf('summary/figures/%s_cat_prop_genes.png', prefix)
	plot_cat_prop_genes = plot_categorical(
	  annotated_regions = ar,
	  annotated_random = ar_rand,
	  x='name', fill='annot_type',
	  x_order = cats_order, fill_order = a_gene_order, position='fill',
	  plot_title = sprintf('%s classification by Annotation', display_type),
	  legend_title = 'Annotations',
	  x_label = sprintf('%s classification', display_type),
	  y_label = 'Proportion')
	ggplot2::ggsave(filename = cat_prop_genes_png, plot = plot_cat_prop_genes, width = 8, height = 8)
}

save.image(file = sprintf('RData/%s_annotatr_analysis.RData', prefix))
