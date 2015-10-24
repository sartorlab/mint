library(ggplot2)

setwd('~/Desktop')

################################################################################
# % methylation and coverage across all sites in sample

  data = read.table('IDH2mut_1_mc_hmc_trim.fastq.gz_bismark.bismark.cov.gz', header=F, sep='\t', quote='', comment.char='', stringsAsFactors=F, col.names=c('chrom','start','end','perc_meth','numC','numT'))
  data$coverage = log10(data$numC + data$numT)

  # % methylation
  per_meth_plot = ggplot(data, aes(perc_meth)) + geom_histogram(binwidth=5, aes(y=..density..)) + theme_bw()
  ggsave(filename='IDH2mut_1_mc_hmc_sample_per_meth.png', plot = per_meth_plot, width=6, height=6, dpi=300)

  # coverage
  coverage_plot = ggplot(data, aes(coverage)) + geom_histogram(aes(y=..density..)) + theme_bw()
  ggsave(filename='IDH2mut_1_mc_hmc_sample_coverage.png', plot = coverage_plot, width=6, height=6, dpi=300)

################################################################################
# Average % methylation across annotations

  files = list.files(pattern='avg_meth')

  annotations = gsub('IDH2mut_1_mc_hmc_avg_meth_','',files)
  annotations = gsub('.bed','',annotations)

  data = lapply(files, read.table, sep='\t', header=F, quote='', comment.char='', stringsAsFactors=F)
  names(data) = annotations

  df = Reduce(rbind, lapply(names(data), function(n){data.frame('annotation'=n,'avg_meth'=data[[n]][,4],stringsAsFactors=F)}))

  avg_meth_plot = ggplot(df, aes(avg_meth)) + geom_histogram(binwidth=5, aes(y=..density..)) + facet_grid(. ~ annotation) + theme_bw()

  ggsave(filename='IDH2mut_1_mc_hmc_avg_meth_annotations.png', plot = avg_meth_plot, width = 8, height = 4, dpi = 300)

################################################################################
# Average coverage across annotations

  files = list.files(pattern='avg_coverage')

  annotations = gsub('IDH2mut_1_mc_hmc_avg_coverage_','',files)
  annotations = gsub('.bed','',annotations)

  data = lapply(files, read.table, sep='\t', header=F, quote='', comment.char='', stringsAsFactors=F)
  names(data) = annotations

  df = Reduce(rbind, lapply(names(data), function(n){data.frame('annotation'=n,'avg_coverage'=log10(data[[n]][,4]),stringsAsFactors=F)}))

  avg_meth_plot = ggplot(df, aes(avg_coverage)) + geom_histogram(aes(y=..density..)) + facet_grid(. ~ annotation) + theme_bw()

  ggsave(filename='IDH2mut_1_mc_hmc_avg_coverage_annotations.png', plot = avg_meth_plot, width = 8, height = 4, dpi = 300)
