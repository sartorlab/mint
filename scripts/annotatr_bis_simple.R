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
sample = gsub('_simple_classification_for_annotatr.txt','', basename(file))
suffix = 'simple_class'

###############################################################
# Read
r = read_bed(
	file = file,
	col.names = FALSE,
	genome = genome,
	stranded = FALSE,
	use.score = FALSE)

###############################################################
# Pick annotations
if(genome %in% c('hg19','hg38','mm9','mm10')) {
	a = c(
		sprintf('%s_cpgs', genome),
		sprintf('%s_basicgenes', genome))
} else {
	# TODO: Add support for insertion of custom annotation path via command line
	stop('Only hg19, hg38, mm9, and mm10 annotations are supported for this part of mint at the moment.')
}

if(genome == 'hg19') {
	a = c(a, 'hg19_enhancers_fantom')
}

###############################################################
# Order annotations
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
			'knownGenes_3UTRs'), sep='_')
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
			'knownGenes_3UTRs'), sep='_')
		a_all_order = c(
			a_cpg_order,
			a_gene_order)
	}
}

cats = unique(r$name)
cat_order = c(
	grep('high', cats, value=TRUE),
	grep('med', cats, value=TRUE),
	grep('low', cats, value=TRUE),
	grep('none', cats, value=TRUE))

###############################################################
# Annotate regions
ar = annotate_regions(
	regions = r,
	annotations = a,
	ignore.strand = TRUE,
	use.score = FALSE)

# Write it
readr::write_tsv(x = ar, path = sprintf('summary/tables/%s_%s_annotations.txt', sample, suffix))

###############################################################
# Summarize annotations
count_annots = summarize_annotations(ar)

# Write it
readr::write_tsv(x = count_annots, path = sprintf('summary/tables/%s_%s_annotation_counts.txt', sample, suffix))

###############################################################
# Visualizations

counts_png = sprintf('summary/figures/%s_%s_counts.png', sample, suffix)
plot_counts = visualize_annotation(
	annotated_regions = ar,
	annotation_order = a_all_order,
	plot_title = sprintf('%s simple class. regions per annotation', sample),
	x_label = 'Annotations',
	y_label = '# CpGs')
ggplot2::ggsave(filename = counts_png, plot = plot_counts, width = 8, height = 8)

cocounts_png = sprintf('summary/figures/%s_%s_cocounts.png', sample, suffix)
plot_cocounts = visualize_coannotations(
	annotated_regions = ar,
	annotation_order = a_all_order,
	plot_title = sprintf('%s simple class. regions in pairs of annotations', sample),
	axes_label = 'Annotations')
ggplot2::ggsave(filename = cocounts_png, plot = plot_cocounts, width = 8, height = 8)

cat_prop_cpgs_png = sprintf('summary/figures/%s_%s_cat_prop_cpgs.png', sample, suffix)
plot_cat_prop_cpgs = visualize_categorical(
  annotated_regions = ar, x='name', fill='annot_type',
  x_order = cat_order, fill_order = a_cpg_order, position='fill',
  plot_title = 'Simple Classification by Annotation',
  legend_title = 'Annotations',
  x_label = 'Simple Classification',
  y_label = 'Proportion')
ggplot2::ggsave(filename = cat_prop_cpgs_png, plot = plot_cat_prop_cpgs, width = 8, height = 8)

cat_prop_genes_png = sprintf('summary/figures/%s_%s_cat_prop_genes.png', sample, suffix)
plot_cat_prop_genes = visualize_categorical(
  annotated_regions = ar, x='name', fill='annot_type',
  x_order = cat_order, fill_order = a_gene_order, position='fill',
  plot_title = 'Simple Classification by Annotation',
  legend_title = 'Annotations',
  x_label = 'Simple Classification',
  y_label = 'Proportion')
ggplot2::ggsave(filename = cat_prop_genes_png, plot = plot_cat_prop_genes, width = 8, height = 8)
