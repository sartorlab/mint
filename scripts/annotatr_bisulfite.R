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
if(grepl('_trimmed_bismark_bt2.CpG_report_for_annotatr.txt', file)) {
	sample = gsub('_trimmed_bismark_bt2.CpG_report_for_annotatr.txt','', basename(file))
	suffix = 'bismark'
} else if (grepl('_methylSig_for_annotatr.txt', file)) {
	sample = gsub('_methylSig_for_annotatr.txt','', basename(file))
	suffix = 'methylSig'
}


###############################################################
# Read
if(suffix == 'bismark') {
	column_names = c('chr','start','end','name','coverage','strand','perc_meth')
	strand = TRUE
} else if (suffix == 'methylSig') {
	column_names = c('chr','start','end','DM_status','pvalue','strand','meth_diff','mu1','mu0')
	strand = FALSE
}

r = read_bed(
	file = file,
	col.names = column_names,
	genome = genome,
	stranded = strand,
	use.score = TRUE)

# sign-log the p-values for volcano plots
if(suffix == 'methylSig') {
	r$pvalue = -log10(r$pvalue)
}

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

cats_order = c('hyper','hypo','noDM')

###############################################################
# Annotate regions
ar = annotate_regions(
	regions = r,
	annotations = a,
	ignore.strand = TRUE,
	use.score = TRUE)

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
plot_counts = plot_annotation(
	annotated_regions = ar,
	annotation_order = a_all_order,
	plot_title = sprintf('%s CpGs per annotation', sample),
	x_label = 'Annotations',
	y_label = '# CpGs')
ggplot2::ggsave(filename = counts_png, plot = plot_counts, width = 8, height = 8)

cocounts_png = sprintf('summary/figures/%s_%s_cocounts.png', sample, suffix)
plot_cocounts = plot_coannotations(
	annotated_regions = ar,
	annotation_order = a_all_order,
	plot_title = sprintf('%s CpGs in pairs of annotations', sample),
	axes_label = 'Annotations')
ggplot2::ggsave(filename = cocounts_png, plot = plot_cocounts, width = 8, height = 8)

if(suffix == 'bismark') {
	coverage_png = sprintf('summary/figures/%s_%s_coverage.png', sample, suffix)
	plot_coverage = plot_numerical(
		tbl = ar,
		x = 'coverage',
		facet = 'annot_type',
		facet_order = a_all_order,
		bin_width = 10,
		plot_title = sprintf('%s coverage over annotations', sample),
		x_label = 'Coverage')
	ggplot2::ggsave(filename = coverage_png, plot = plot_coverage, width = 8, height = 8)

	percmeth_png = sprintf('summary/figures/%s_%s_percmeth.png', sample, suffix)
	plot_percmeth = plot_numerical(
		tbl = ar,
		x = 'perc_meth',
		facet = 'annot_type',
		facet_order = a_all_order,
		bin_width = 5,
		plot_title = sprintf('%s percent meth. over annotations', sample),
		x_label = 'Percent Methylation')
	ggplot2::ggsave(filename = percmeth_png, plot = plot_percmeth, width = 8, height = 8)
}

if(suffix == 'methylSig') {
	methdiff_png = sprintf('summary/figures/%s_%s_methdiff.png', sample, suffix)
	plot_methdiff = plot_numerical(
		tbl = ar,
		x = 'meth_diff',
		facet = 'annot_type',
		facet_order = a_all_order,
		bin_width = 5,
		plot_title = sprintf('%s methylation difference over annotations', sample),
		x_label = 'Methylation Difference (Group1 - Group0)')
	ggplot2::ggsave(filename = methdiff_png, plot = plot_methdiff, width = 8, height = 8)

	volcano_png = sprintf('summary/figures/%s_%s_volcano.png', sample, suffix)
	plot_volcano = plot_numerical(
		tbl = ar,
		x = 'meth_diff',
		y = 'pvalue',
		facet = 'annot_type',
		facet_order = a_all_order,
		plot_title = sprintf('%s meth. diff. vs -log10(pval)', sample),
		x_label = 'Methylation Difference (Group1 - Group0)',
		y_label = '-log10(pvalue)')
	ggplot2::ggsave(filename = volcano_png, plot = plot_volcano, width = 8, height = 8)

	# Regions split by category and stacked by CpG annotations (count)
	cat_count_cpgs_png = sprintf('summary/figures/%s_%s_DMstatus_count_cpgs.png', sample, suffix)
	plot_cat_count_cpgs = plot_categorical(
	  annotated_regions = ar, x='DM_status', fill='annot_type',
	  x_order = cats_order, fill_order = a_cpg_order, position='stack',
	  plot_title = sprintf('%s DM status by Annotation', sample),
	  legend_title = 'Annotations',
	  x_label = sprintf('%s DM status', sample),
	  y_label = 'Count')
	ggplot2::ggsave(filename = cat_count_cpgs_png, plot = plot_cat_count_cpgs, width = 8, height = 8)

	# Regions split by category and stacked by knownGene annotations (count)
	cat_count_genes_png = sprintf('summary/figures/%s_%s_DMstatus_count_genes.png', sample, suffix)
	plot_cat_count_genes = plot_categorical(
	  annotated_regions = ar, x='DM_status', fill='annot_type',
	  x_order = cats_order, fill_order = a_gene_order, position='stack',
	  plot_title = sprintf('%s DM status by Annotation', sample),
	  legend_title = 'Annotations',
	  x_label = sprintf('%s DM status', sample),
	  y_label = 'Count')
	ggplot2::ggsave(filename = cat_count_genes_png, plot = plot_cat_count_genes, width = 8, height = 8)

	# Regions split by category and filled by CpG annotations (prop)
	cat_prop_cpgs_png = sprintf('summary/figures/%s_%s_DMstatus_prop_cpgs.png', sample, suffix)
	plot_cat_prop_cpgs = plot_categorical(
	  annotated_regions = ar, x='DM_status', fill='annot_type',
	  x_order = cats_order, fill_order = a_cpg_order, position='fill',
	  plot_title = sprintf('%s DM status by Annotation', sample),
	  legend_title = 'Annotations',
	  x_label = sprintf('%s DM status', sample),
	  y_label = 'Proportion')
	ggplot2::ggsave(filename = cat_prop_cpgs_png, plot = plot_cat_prop_cpgs, width = 8, height = 8)

	# Regions split by category and filled by knownGene annotations (prop)
	cat_prop_genes_png = sprintf('summary/figures/%s_%s_DMstatus_prop_genes.png', sample, suffix)
	plot_cat_prop_genes = plot_categorical(
	  annotated_regions = ar, x='DM_status', fill='annot_type',
	  x_order = cats_order, fill_order = a_gene_order, position='fill',
	  plot_title = sprintf('%s DM status by Annotation', sample),
	  legend_title = 'Annotations',
	  x_label = sprintf('%s DM status', sample),
	  y_label = 'Proportion')
	ggplot2::ggsave(filename = cat_prop_genes_png, plot = plot_cat_prop_genes, width = 8, height = 8)

}
