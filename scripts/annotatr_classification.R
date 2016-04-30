library(annotatr)
library(readr)
library(ggplot2)
library(optparse)

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
} else {
	stop('Invalid classification type. Must be one of simple, sample, or compare.')
}

###############################################################
# Read
r = read_bed(
	file = file,
	col.names = FALSE,
	genome = genome,
	stranded = FALSE,
	use.score = FALSE)

# Remove unclassifiable regions
r = r[r$name != 'unclassifiable']

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
		'hyper',
		'hypo')
}

###############################################################
# Annotate regions
ar = annotate_regions(
	regions = r,
	annotations = a,
	ignore.strand = TRUE,
	use.score = FALSE)

# Write it
readr::write_tsv(x = ar, path = sprintf('summary/tables/%s_annotations.txt', prefix))

###############################################################
# Summarize annotations
count_annots = summarize_annotations(ar)

# Write it
readr::write_tsv(x = count_annots, path = sprintf('summary/tables/%s_annotation_counts.txt', prefix))

###############################################################
# Visualizations

# Regular barplot of regions in annotations
counts_png = sprintf('summary/figures/%s_counts.png', prefix)
plot_counts = plot_annotation(
	annotated_regions = ar,
	annotation_order = a_all_order,
	plot_title = sprintf('%s regions per annotation', prefix),
	x_label = 'Annotations',
	y_label = '# Regions')
ggplot2::ggsave(filename = counts_png, plot = plot_counts, width = 8, height = 8)

# Other classifications are too big for this
if(class_type == 'simple' || class_type == 'PePr') {
	# Heatmap of regions in pairs of annotations
	cocounts_png = sprintf('summary/figures/%s_cocounts.png', prefix)
	plot_cocounts = plot_coannotations(
		annotated_regions = ar,
		annotation_order = a_all_order,
		plot_title = sprintf('%s regions in pairs of annotations', prefix),
		axes_label = 'Annotations')
	ggplot2::ggsave(filename = cocounts_png, plot = plot_cocounts, width = 8, height = 8)
}

# Regions split by category and stacked by CpG annotations (count)
cat_count_cpgs_png = sprintf('summary/figures/%s_cat_count_cpgs.png', prefix)
plot_cat_count_cpgs = plot_categorical(
  annotated_regions = ar, x='name', fill='annot_type',
  x_order = cats_order, fill_order = a_cpg_order, position='stack',
  plot_title = sprintf('%s classification by Annotation', display_type),
  legend_title = 'Annotations',
  x_label = sprintf('%s classification', display_type),
  y_label = 'Count')
ggplot2::ggsave(filename = cat_count_cpgs_png, plot = plot_cat_count_cpgs, width = 8, height = 8)

# Regions split by category and stacked by knownGene annotations (count)
cat_count_genes_png = sprintf('summary/figures/%s_cat_count_genes.png', prefix)
plot_cat_count_genes = plot_categorical(
  annotated_regions = ar, x='name', fill='annot_type',
  x_order = cats_order, fill_order = a_gene_order, position='stack',
  plot_title = sprintf('%s classification by Annotation', display_type),
  legend_title = 'Annotations',
  x_label = sprintf('%s classification', display_type),
  y_label = 'Count')
ggplot2::ggsave(filename = cat_count_genes_png, plot = plot_cat_count_genes, width = 8, height = 8)

# Regions split by category and filled by CpG annotations (prop)
cat_prop_cpgs_png = sprintf('summary/figures/%s_cat_prop_cpgs.png', prefix)
plot_cat_prop_cpgs = plot_categorical(
  annotated_regions = ar, x='name', fill='annot_type',
  x_order = cats_order, fill_order = a_cpg_order, position='fill',
  plot_title = sprintf('%s classification by Annotation', display_type),
  legend_title = 'Annotations',
  x_label = sprintf('%s classification', display_type),
  y_label = 'Proportion')
ggplot2::ggsave(filename = cat_prop_cpgs_png, plot = plot_cat_prop_cpgs, width = 8, height = 8)

# Regions split by category and filled by knownGene annotations (prop)
cat_prop_genes_png = sprintf('summary/figures/%s_cat_prop_genes.png', prefix)
plot_cat_prop_genes = plot_categorical(
  annotated_regions = ar, x='name', fill='annot_type',
  x_order = cats_order, fill_order = a_gene_order, position='fill',
  plot_title = sprintf('%s classification by Annotation', display_type),
  legend_title = 'Annotations',
  x_label = sprintf('%s classification', display_type),
  y_label = 'Proportion')
ggplot2::ggsave(filename = cat_prop_genes_png, plot = plot_cat_prop_genes, width = 8, height = 8)
