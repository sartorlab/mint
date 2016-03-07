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
sample = gsub('_trimmed.fq.gz_bismark_bt2.CpG_report_for_annotatr.txt','', basename(file))
suffix = 'bismark'

###############################################################
# Read
r = read_bed(
	file = file,
	col.names = c('chr','start','end','name','coverage','strand','perc_meth'),
	genome = genome,
	stranded = TRUE,
	use.score = TRUE)

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
plot_counts = visualize_annotation(
	annotated_regions = ar,
	annotation_order = a_all_order,
	plot_title = sprintf('%s CpGs per annotation', sample),
	x_label = 'Annotations',
	y_label = '# CpGs')
ggplot2::ggsave(filename = counts_png, plot = plot_counts, width = 8, height = 8)

cocounts_png = sprintf('summary/figures/%s_%s_cocounts.png', sample, suffix)
plot_cocounts = visualize_coannotations(
	annotated_regions = ar,
	annotation_order = a_all_order,
	plot_title = sprintf('%s CpGs in pairs of annotations', sample),
	axes_label = 'Annotations')
ggplot2::ggsave(filename = cocounts_png, plot = plot_cocounts, width = 8, height = 8)

coverage_png = sprintf('summary/figures/%s_%s_coverage.png', sample, suffix)
plot_coverage = visualize_numerical(
	tbl = ar,
	x = 'coverage',
	facet = 'annot_type',
	facet_order = a_all_order,
	bin_width = 10,
	plot_title = sprintf('%s coverage over annotations', sample),
	x_label = 'Coverage')
ggplot2::ggsave(filename = coverage_png, plot = plot_coverage, width = 8, height = 8)

percmeth_png = sprintf('summary/figures/%s_%s_percmeth.png', sample, suffix)
plot_percmeth = visualize_numerical(
	tbl = ar,
	x = 'perc_meth',
	facet = 'annot_type',
	facet_order = a_all_order,
	bin_width = 5,
	plot_title = sprintf('%s percent meth. over annotations', sample),
	x_label = 'Percent Methylation')
ggplot2::ggsave(filename = percmeth_png, plot = plot_percmeth, width = 8, height = 8)
