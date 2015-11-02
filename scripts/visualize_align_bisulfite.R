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
  setwd(sprintf('~/latte/mint/%s/', project))

  summary_files = list.files('analysis/summary/tables', full.names=T)
  summary_files = grep(humanID, summary_files, value=T)

  bis_cov = sprintf('analysis/bismark_extractor_calls/%s_trim.fastq.gz_bismark.bismark.cov.gz', humanID)

  files = list(
    cpg = grep('cpg', summary_files, value=T),
    windows_5kb = grep('windows', summary_files, value=T),
    prom_enh = grep('(promoters)|(enhancers)', summary_files, value=T))

################################################################################
# % methylation and coverage across all sites in sample
# Cannot include this in the next block because we get data from one file instead of two

  message(sprintf('Plotting percent methylation and coverage across all sites in %s', humanID))
  data = read.table(bis_cov,
    header=F, sep='\t', quote='', comment.char='', stringsAsFactors=F,
    col.names=c('chrom','start','end','perc_meth','numC','numT'))
  data$coverage = data$numC + data$numT

  # % methylation
  png_pm = sprintf('analysis/summary/figures/%s_bisulfite_all_sites_perc_meth.png', humanID)
  plot_pm =
    ggplot(data, aes(perc_meth)) +
    geom_histogram(binwidth=5, aes(y=..density..)) +
    ggtitle('% Methylation Across All Sites') +
    theme_bw()
  ggsave(filename = png_pm, plot = plot_pm, width=6, height=6, dpi=300)

  # coverage
  png_cov = sprintf('analysis/summary/figures/%s_bisulfite_all_sites_cov.png', humanID)
  plot_cov =
    ggplot(data, aes(log10(coverage))) +
    geom_histogram(aes(y=..density..)) +
    ggtitle('Coverage Across All Sites') +
    theme_bw()
  ggsave(filename = png_cov, plot = plot_cov, width=6, height=6, dpi=300)

################################################################################
# Average % methylation and coverage across CpG annotations and promoter/enhancers

  lapply(names(files), function(annot){
    message(sprintf('Plotting average percent methylation and average coverage across %s annotations in %s', annot, humanID))

    files_pm = grep('avg_methylation', files[[annot]], value=T)
    files_cov = grep('avg_coverage', files[[annot]], value=T)

    annotations = gsub(sprintf('analysis/summary/tables/%s_avg_methylation_',humanID),'',files_pm)
    annotations = gsub('.bed','',annotations)

    data_pm = lapply(files_pm, read.table, sep='\t', header=F, quote='', comment.char='', stringsAsFactors=F)
    names(data_pm) = annotations

    data_cov = lapply(files_cov, read.table, sep='\t', header=F, quote='', comment.char='', stringsAsFactors=F)
    names(data_cov) = annotations

    df_pm = Reduce(rbind,
      lapply(names(data_pm), function(n){
        data.frame(
          'annotation'=n,
          'avg_meth'=data_pm[[n]][,5],
          stringsAsFactors=F)}))

    df_cov = Reduce(rbind,
      lapply(names(data_cov), function(n){
        data.frame(
          'annotation'=n,
          'avg_coverage'=log10(data_cov[[n]][,5]),
          stringsAsFactors=F)}))

    if (length(data_pm) == 2) {
      width = 6
      height = 4
    } else if (length(data_pm) > 2) {
      width = 8
      height = 4
    } else {
      width = 4
      height = 4
    }

    png_pm = sprintf('analysis/summary/figures/%s_bisulfite_%s_avg_perc_meth.png', humanID, annot)
    plot_pm =
      ggplot(df_pm, aes(avg_meth)) +
      geom_histogram(binwidth=5, aes(y=..density..)) +
      facet_grid(. ~ annotation) +
      ggtitle('Avg. % Meth. Across Annotations') +
      theme_bw()
    ggsave(filename = png_pm, plot = plot_pm, width = width, height = height, dpi = 300)

    png_cov = sprintf('analysis/summary/figures/%s_bisulfite_%s_avg_cov.png', humanID, annot)
    plot_cov =
      ggplot(df_cov, aes(avg_coverage)) +
      geom_histogram(aes(y=..density..)) +
      facet_grid(. ~ annotation) +
      ggtitle('Avg. Cov. Across Annotations') +
      theme_bw()
    ggsave(filename = png_cov, plot = plot_cov, width = width, height = height, dpi = 300)

  })
