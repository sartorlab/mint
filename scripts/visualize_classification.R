library(ggplot2)
library(optparse)

# Parse arguments
option_list = list(
    make_option('--project', type='character'),
    make_option('--classAnnot', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

project = opt$project
classAnnot = opt$classAnnot

# Set working directory and collect relevant files
  setwd(sprintf('~/latte/mint/%s/', project))

  baseClassAnnot = basename(classAnnot)
  baseClassAnnot = gsub('_annot.bed', '', baseClassAnnot)

  annots = list(
    cpgs = c('cpg_islands','cpg_shores','cpg_shelves','cpg_inter'),
    prom_enh = c('promoters','enhancers')
  )

################################################################################
# Distribution of classifications overall and in CpG annotations
# Read-in takes a while, but we can't really avoid that

  data = read.table(classAnnot,
    header=F, sep='\t', comment.char='', quote='', stringsAsFactors=F,
    col.names=c('chr','start','end','class','annot_chr','annot_start','annot_end','annot_type'),
    colClasses=c('NULL','NULL','NULL','character','NULL','NULL','NULL','character'))
  data$annot_type[grepl('promoters',data$annot_type)] = 'promoters'

  data_list = lapply(unique(data$class), function(c) {
    table(subset(data, class == c))
  })
  names(data_list) = sapply(data_list, row.names)

  df = Reduce(rbind, lapply(names(data_list), function(n){
    data.frame(
      'type' = n,
      'annot' = colnames(data_list[[n]]),
      'prop'= as.numeric(data_list[[n]]),
      stringsAsFactors=F)
  }))

  df_sub = subset(df, !(type %in% c('unclassifiable','noDM')))

  lapply(names(annots), function(a){
    message(sprintf('Plotting prevalence of classes across %s annotations in %s', a, baseClassAnnot))

    df_annot = subset(df_sub, annot %in% annots[[a]])

    png_all_counts = sprintf('analysis/summary/figures/%s_%s_class_all_counts_barplot.png', baseClassAnnot, a)
    plot_all_counts =
      ggplot(df_annot, aes(x=type, y=prop)) +
      geom_bar(stat='identity') +
      ggtitle('Counts of Overall Classifications') +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    ggsave(filename=png_all_counts, plot=plot_all_counts, width=6, height=4, dpi=300)

    png_annot_counts = sprintf('analysis/summary/figures/%s_%s_class_annot_counts_barplot.png', baseClassAnnot, a)
    plot_annot_counts =
      ggplot(df_annot, aes(x=type, y=prop, fill=annot)) +
      geom_bar(stat='identity') +
      ggtitle('Counts of Classifications in annot.') +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    ggsave(filename=png_annot_counts, plot=plot_annot_counts, width=6, height=4, dpi=300)

    png_annot_props = sprintf('analysis/summary/figures/%s_%s_class_annot_props_barplot.png', baseClassAnnot, a)
    plot_annot_props =
      ggplot(df_annot, aes(x=type, y=prop, fill=annot)) +
      geom_bar(stat='identity', position='fill') + ggtitle('Prop. of Classifications in annot.') +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    ggsave(filename=png_annot_props, plot=plot_annot_props, width=6, height=4, dpi=300)

  })
