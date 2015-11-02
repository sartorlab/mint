library(ggplot2)
library(optparse)

# Parse arguments
option_list = list(
    make_option('--project', type='character'),
    make_option('--comparison', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

project = opt$project
comparison = opt$comparison

# Set working directory and collect relevant files
  setwd(sprintf('../%s/', project))

################################################################################
# Distribution of PePr peak widths

  data = read.table(sprintf('analysis/summary/tables/%s_PePr_peak_widths.bed', comparison),
    header=F, sep='\t', comment.char='', quote='', stringsAsFactors=F,
    col.names=c('chr','start','end','direction','width'))

  data_list = list(
    DM = data$width,
    pos_DM = subset(data, direction == 'up')$width,
    neg_DM = subset(data, direction == 'down')$width
  )

  df = Reduce(rbind, lapply(names(data_list), function(n){
    data.frame(
      'type' = n,
      'width'= data_list[[n]],
      stringsAsFactors=F)
  }))

  png_widths = sprintf('analysis/summary/figures/%s_PePr_peak_widths.png', comparison)
  plot_widths =
    ggplot(df, aes(log10(width))) +
    geom_histogram(aes(y=..density..)) +
    facet_grid(. ~ type) +
    ggtitle('PePr Peak Widths') +
    theme_bw()
  ggsave(filename = png_widths, plot = plot_widths, width=6, height=4, dpi=300)

################################################################################
# Distribution of PePr DM regions in CpG annotations

  annots = list(
    cpgs = c('cpg_islands','cpg_shores','cpg_shelves','cpg_inter'),
    prom_enh = c('promoters','enhancers')
  )

  data = read.table(sprintf('analysis/summary/tables/%s_PePr_annot.bed', comparison),
    header=F, sep='\t', comment.char='', quote='', stringsAsFactors=F,
    col.names=c('chr','start','end','direction','annot_chr','annot_start','annot_end','annot_type'))
  data$annot_type[grepl('promoters',data$annot_type)] = 'promoters'

  data_list = list(
    DM = table(data$annot_type),
    pos_DM = table(subset(data, direction == 'up')$annot_type),
    neg_DM = table(subset(data, direction == 'down')$annot_type)
  )

  df = Reduce(rbind, lapply(names(data_list), function(n){
    data.frame(
      'type' = n,
      'annot' = names(data_list[[n]]),
      'prop'= as.numeric(data_list[[n]]),
      stringsAsFactors=F)
  }))

  lapply(names(annots), function(a){
    message(sprintf('Plotting prevalence of PePr peaks across %s annotations in %s', a, comparison))
    df_annot = subset(df, annot %in% annots[[a]])

    if (length(annots[[a]]) == 2) {
      width = 6
      height = 4
    } else if (length(annots[[a]]) > 2) {
      width = 8
      height = 4
    } else {
      width = 4
      height = 4
    }

    # scale_fill_manual(values=c('inter'="#00008B", 'shelf'="#0000FF", 'shore'="#F4A460", 'island'="#003300"))

    png_counts = sprintf('analysis/summary/figures/%s_PePr_%s_annot_barplot_counts.png', comparison, a)
    plot_counts =
      ggplot(df_annot, aes(x=type, y=prop, fill=annot)) +
      geom_bar(stat='identity') +
      ggtitle(sprintf('Counts of PePr peaks in %s (%s)', a, comparison))
    ggsave(filename = png_counts, plot = plot_counts, width=width, height=height, dpi=300)

    png_props = sprintf('analysis/summary/figures/%s_PePr_%s_annot_barplot_props.png', comparison, a)
    plot_props =
      ggplot(df_annot, aes(x=type, y=prop, fill=annot)) +
      geom_bar(stat='identity', position='fill') +
      ggtitle(sprintf('Prop. of PePr peaks in %s (%s)', a, comparison))
    ggsave(filename = png_props, plot = plot_props, width=width, height=height, dpi=300)

  })
