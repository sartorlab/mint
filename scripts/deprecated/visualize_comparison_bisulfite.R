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
  setwd(sprintf('~/latte/mint/%s/', project))

  summary_files = list.files('analysis/summary/tables', full.names=T)
  summary_files = grep(comparison, summary_files, value=T)
  summary_files = grep('methylSig', summary_files, value=T)

  files = list(
    cpg = grep('cpg', summary_files, value=T),
    windows_5kb = grep('windows', summary_files, value=T),
    prom_enh = grep('(promoters)|(enhancers)', summary_files, value=T))

################################################################################
# % methylation change across all DMRs (methylSig)
  message('Plotting methylation change across all sites and all DMs.')
  data = read.table(sprintf('analysis/methylsig_calls/%s.txt', comparison),
    header=T, sep='\t', quote='', comment.char='', stringsAsFactors=F)
  sig = subset(data, pvalue < 0.05)

  DM = data.frame(
    'annotation'='all_tested',
    'meth_diff'=data$meth.diff,
    'pvalue'=data$pvalue,
    stringsAsFactors=F)
  DM_sig = data.frame(
    'annotation'='DM',
    'meth_diff'=sig$meth.diff,
    'pvalue'=sig$pvalue,
    stringsAsFactors=F)

  DM_groups = data.frame(
    'annotation'='all_tested',
    'meth1'=data$mu1,
    'meth0'=data$mu0,
    stringsAsFactors=F)
  DM_groups_sig = data.frame(
    'annotation'='DM',
    'meth1'=sig$mu1,
    'meth0'=sig$mu0,
    stringsAsFactors=F)

  df = rbind(DM, DM_sig)
  df_groups = rbind(DM_groups, DM_groups_sig)

  png_volcano = sprintf('analysis/summary/figures/%s_methylSig_diff_meth_pval_volcano.png', comparison)
  plot_volcano =
    ggplot(DM, aes(meth_diff,-log10(pvalue))) +
    geom_point(alpha=1/8, size=1) +
    geom_abline(intercept=-log10(0.05), slope=0, color='red') +
    ggtitle('% Methylation Difference vs -log10(P-value)') + theme(plot.title = element_text(size=8)) +
    theme_bw()
  ggsave(filename = png_volcano, plot = plot_volcano, width = 5, height = 5, dpi = 300)

  png_dm_groups = sprintf('analysis/summary/figures/%s_methylSig_diff_meth_groups_heatmap.png', comparison)
  plot_dm_groups =
    ggplot(df_groups, aes(meth1, meth0)) +
    geom_point(alpha=1/8, size=1) +
    geom_abline(intercept=0, slope=1, color='red') +
    facet_grid(. ~ annotation) +
    ggtitle('% Methylation of Groups (DMR-Matched)')
  ggsave(filename = png_dm_groups, plot = plot_dm_groups, width = 8, height = 4, dpi = 300)

################################################################################
# Average % methylation change across annotations (methylSig)
# These cannot be changed to volcano plots because we're taking the average
# over the CpG annotations of the methylation differences.

  # lapply(names(files), function(annot) {
  #   message(sprintf('Plotting average difference in methylation across %s annotations in %s', annot, comparison))
  #
  #   files_diff = files[[annot]]
  #
  #   annotations = gsub(sprintf('analysis/summary/tables/%s_methylSig_avg_diff_meth_in_DM_',comparison),'',files_diff)
  #   annotations = gsub('.bed','',annotations)
  #
  #   data_diff = lapply(files_diff, read.table, sep='\t', header=F, quote='', comment.char='', stringsAsFactors=F)
  #   names(data_diff) = annotations
  #
  #   if (length(data_diff) == 2) {
  #     width = 6
  #     height = 4
  #   } else if (length(data_diff) > 2) {
  #     width = 8
  #     height = 4
  #   } else {
  #     width = 4
  #     height = 4
  #   }
  #
  #   df_diff = Reduce(rbind,
  #     lapply(names(data_diff), function(n){
  #       data.frame(
  #         'annotation'=n,
  #         'diff_meth'=data_diff[[n]][,5],
  #         stringsAsFactors=F)}))
  #
  #   png_diff = sprintf('analysis/summary/figures/%s_methylSig_%s_diff_meth_in_DM.png', comparison, annot)
  #   plot_diff =
  #     ggplot(df_diff, aes(diff_meth)) +
  #     geom_histogram(binwidth=5, aes(y=..density..)) +
  #     facet_grid(. ~ annotation) +
  #     ggtitle('% Methylation Difference Across Annotations') +
  #     theme_bw()
  #   ggsave(filename = png_diff, plot = plot_diff, width = width, height = height, dpi = 300)
  #
  # })

################################################################################
# Distribution of methylSig DM regions in annotations

  annots = list(
    cpgs = c('cpg_islands','cpg_shores','cpg_shelves','cpg_inter'),
    prom_enh = c('promoters','enhancers')
  )

  data = read.table(sprintf('analysis/summary/tables/%s_methylSig_annot.bed', comparison),
    header=F, sep='\t', comment.char='', quote='', stringsAsFactors=F,
    col.names=c('chr','start','end','pval','qval','diff_meth','mu1',
      'mu0','annot_chr','annot_start','annot_end','annot_type'))
  data$annot_type[grepl('promoters',data$annot_type)] = 'promoters'

  data_list = list(
    all_tested = table(data$annot_type),
    DM = table(subset(data, pval < 0.05)$annot_type),
    pos_DM = table(subset(data, pval < 0.05 & diff_meth > 0)$annot_type),
    neg_DM = table(subset(data, pval < 0.05 & diff_meth < 0)$annot_type)
  )

  df = Reduce(rbind, lapply(names(data_list), function(n){
    data.frame(
      'type' = n,
      'annot' = names(data_list[[n]]),
      'prop'= as.numeric(data_list[[n]]),
      stringsAsFactors=F)
  }))

  lapply(names(annots), function(a){
    message(sprintf('Plotting prevalence of methylSig DM regions across %s annotations in %s', a, comparison))
    df_annot = subset(df, annot %in% annots[[a]])
    data_annot = subset(data, annot_type %in% annots[[a]])

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

    png_volcano = sprintf('analysis/summary/figures/%s_methylSig_%s_annot_volcanos.png', comparison, a)
    plot_volcano =
      ggplot(data_annot, aes(diff_meth,-log10(pval))) +
      geom_point(alpha=1/8, size=1) +
      geom_abline(intercept=-log10(0.05), slope=0, color='red') +
      facet_grid(. ~ annot_type) +
      ggtitle(sprintf('Pct. Methylation Difference vs -log10(P-value) in %s',a)) +
      theme(plot.title = element_text(size=8)) +
      theme_bw()
    ggsave(filename = png_volcano, plot = plot_volcano, width = width, height = height, dpi = 300)

    png_counts = sprintf('analysis/summary/figures/%s_methylSig_%s_annot_barplot_counts.png', comparison, a)
    plot_counts =
      ggplot(df_annot, aes(x=type, y=prop, fill=annot)) +
      geom_bar(stat='identity') +
      ggtitle(sprintf('Counts of mSig regions in %s (%s)', a, comparison))
    ggsave(filename = png_counts, plot = plot_counts, width=width, height=height, dpi=300)

    png_props = sprintf('analysis/summary/figures/%s_methylSig_%s_annot_barplot_props.png', comparison, a)
    plot_props =
      ggplot(df_annot, aes(x=type, y=prop, fill=annot)) +
      geom_bar(stat='identity', position='fill') +
      ggtitle(sprintf('Prop. of mSig regions in %s (%s)', a, comparison))
    ggsave(filename = png_props, plot = plot_props, width=width, height=height, dpi=300)

  })
