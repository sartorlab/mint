library(ggplot2)
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
  setwd(sprintf('../%s/', project))

################################################################################
# Distribution of PePr peak widths

  data = read.table(sprintf('analysis/summary/tables/%s_macs2_peak_widths.bed', humanID),
    header=F, sep='\t', comment.char='', quote='', stringsAsFactors=F,
    col.names=c('chr','start','end','width'))

  png_widths = sprintf('analysis/summary/figures/%s_macs2_peak_widths.png', humanID)
  plot_widths =
    ggplot(data, aes(log10(width))) +
    geom_histogram(aes(y=..density..)) +
    ggtitle('MACS2 Peak Widths') +
    theme_bw()
  ggsave(filename = png_widths, plot = plot_widths, width=6, height=4, dpi=300)

################################################################################
# Distribution of macs2 regions in CpG annotations

  annots = list(
    cpgs = c('cpg_islands','cpg_shores','cpg_shelves','cpg_inter'),
    prom_enh = c('promoters','enhancers')
  )

  data = read.table(sprintf('analysis/summary/tables/%s_macs2_annot.bed',humanID),
    header=F, sep='\t', comment.char='', quote='', stringsAsFactors=F,
    colClasses=c('NULL','NULL','NULL','NULL','NULL','NULL','numeric','numeric','numeric','NULL','NULL','NULL','NULL','character'),
    col.names=c('chr','start','end','name','score','strand','signal','pval','qval','peak','annot_chr','annot_start','annot_end','annot_type'))
  data$annot_type[grepl('promoters',data$annot_type)] = 'promoters'

  lapply(names(annots), function(a){
    message('Plotting macs2 peaks in annotations and signalValue in %s annotations for ', a, humanID)
    df = subset(data, annot_type %in% annots[[a]])

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

    png_macs2_counts = sprintf('analysis/summary/figures/%s_macs2_%s_annot_barplot_counts.png', humanID, a)
    plot_macs2_count =
      ggplot(df, aes(annot_type)) +
      geom_bar() +
      ggtitle(sprintf('Counts of macs2 regions in annot. (%s)', humanID)) +
      theme_bw()
    ggsave(filename=png_macs2_counts, plot=plot_macs2_count, width=width, height=height, dpi=300)

    png_macs2_signal = sprintf('analysis/summary/figures/%s_macs2_%s_annot_signalValue.png', humanID, a)
    plot_macs2_signal =
      ggplot(df, aes(log2(signal))) +
      geom_histogram(aes(y=..density..)) +
      facet_grid(. ~ annot_type) +
      ggtitle(sprintf('signalValue in annot. (%s)', humanID)) +
      theme_bw()
    ggsave(filename=png_macs2_signal, plot=plot_macs2_signal, width=width, height=height, dpi=300)
  })
