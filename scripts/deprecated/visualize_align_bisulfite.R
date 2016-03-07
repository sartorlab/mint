library(annotatr)
library(optparse)

# Parse arguments
  option_list = list(
      make_option('--project', type='character'),
      make_option('--humanID', type='character')
  )
  opt = parse_args(OptionParser(option_list=option_list))

  project = opt$project
  humanID = opt$humanID

# Set working directory and collect relevant files
  setwd(sprintf('~/latte/mint/%s/', project))

  summary_tables_dir = './analysis/summary/tables'

  bis_cov = sprintf('./analysis/bismark_extractor_calls/%s_trimmed.fq.gz_bismark_bt2.bismark.cov.gz', humanID)
  annot_cov = sprintf('./analysis/bismark_extractor_calls/%s_trimmed.fq.gz_bismark_bt2.CpG_report_for_annotatr.txt.gz', humanID)

################################################################################
# % methylation and coverage

  # Read
  bis_regions = read_bed(
    annot_cov,
    col.names = c('chr','start','end','name','coverage','strand','perc_meth'),
    genome = 'hg19',
    stranded = TRUE,
    use.score = TRUE)

  # Annotate
  annotations = c('hg19_basicgenes','hg19_cpgs','hg19_enhancers_fantom')
  bis_annotated = annotate_regions(
    regions = bis_regions,
    annotations = annotations,
    ignore.strand = TRUE,
    use.score = TRUE)

  # Setup order vectors
  all_order = c(
    'hg19_enhancers_fantom',
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter',
    'hg19_knownGenes_1to5kb',
    'hg19_knownGenes_promoters',
    'hg19_knownGenes_5UTRs',
    'hg19_knownGenes_exons',
    'hg19_knownGenes_introns',
    'hg19_knownGenes_3UTRs')
  cpgs_order = c(
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_cpg_shelves',
    'hg19_cpg_inter')
  genes_order = c(
    'hg19_knownGenes_1to5kb',
    'hg19_knownGenes_promoters',
    'hg19_knownGenes_5UTRs',
    'hg19_knownGenes_exons',
    'hg19_knownGenes_introns',
    'hg19_knownGenes_3UTRs')

  # Visualize
    # Single annotation counts
    bis_vis_annot_counts = visualize_annotation(
      annotated_regions = bis_annotated,
      annotation_order = all_order,
      plot_title = sprintf('Num. %s CpGs in Annotations', humanID),
      x_label = 'Annotations',
      y_label = '# CpGs')

    # Co-occurring annotation counts
    bis_vis_coannot_counts = visualize_coannotations(
      annotated_regions = bis_annotated,
      annotation_order = all_order,
      plot_title = sprintf('Num. %s CpGs in Co-occurring Annotations', humanID),
      axes_label = 'Annotations')

    # Coverage across Annotations
    bis_vis_cov = visualize_numerical(
      tbl = bis_annotated,
      x = 'coverage',
      facet = 'annot_type',
      facet_order = all_order,
      bin_width = 30,
      plot_title = sprintf('%s CpG Coverage in Annotations', humanID),
      x_label = 'Coverage',
      y_label = 'Density')

    # Coverage across Annotations
    bis_vis_meth = visualize_numerical(
      tbl = bis_annotated,
      x = 'perc_meth',
      facet = 'annot_type',
      facet_order = all_order,
      bin_width = 5,
      plot_title = sprintf('%s CpG Perc. Meth. in Annotations', humanID),
      x_label = 'Percent Methylation',
      y_label = 'Density')
