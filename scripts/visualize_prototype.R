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
# Average % methylation across annotations (bismark.cov.gz)

  files = list.files(pattern='avg_meth')

  annotations = gsub('IDH2mut_1_mc_hmc_avg_meth_','',files)
  annotations = gsub('.bed','',annotations)

  data = lapply(files, read.table, sep='\t', header=F, quote='', comment.char='', stringsAsFactors=F)
  names(data) = annotations

  df = Reduce(rbind, lapply(names(data), function(n){data.frame('annotation'=n,'avg_meth'=data[[n]][,4],stringsAsFactors=F)}))

  avg_meth_plot = ggplot(df, aes(avg_meth)) + geom_histogram(binwidth=5, aes(y=..density..)) + facet_grid(. ~ annotation) + theme_bw()

  ggsave(filename='IDH2mut_1_mc_hmc_avg_meth_annotations.png', plot = avg_meth_plot, width = 8, height = 4, dpi = 300)

################################################################################
# Average coverage across annotations (bismark.cov.gz)

  files = list.files(pattern='avg_coverage')

  annotations = gsub('IDH2mut_1_mc_hmc_avg_coverage_','',files)
  annotations = gsub('.bed','',annotations)

  data = lapply(files, read.table, sep='\t', header=F, quote='', comment.char='', stringsAsFactors=F)
  names(data) = annotations

  df = Reduce(rbind, lapply(names(data), function(n){data.frame('annotation'=n,'avg_coverage'=log10(data[[n]][,4]),stringsAsFactors=F)}))

  avg_meth_plot = ggplot(df, aes(avg_coverage)) + geom_histogram(aes(y=..density..)) + facet_grid(. ~ annotation) + theme_bw()

  ggsave(filename='IDH2mut_1_mc_hmc_avg_coverage_annotations.png', plot = avg_meth_plot, width = 8, height = 4, dpi = 300)

################################################################################
# Average % methylation change across all DMRs (methylSig)

  data = read.table('IDH2mut_v_NBM.txt', header=T, sep='\t', quote='', comment.char='', stringsAsFactors=F)
  sig = subset(data, pvalue < 0.05)

  avg_DM = data.frame('annotation'='all_tested', 'meth_diff'=data$meth.diff, stringsAsFactors=F)
  avg_DM_sig = data.frame('annotation'='significantly_DM', 'meth_diff'=sig$meth.diff, stringsAsFactors=F)

  region_DM = data.frame('annotation'='all_tested', 'meth1'=data$mu1, 'meth0'=data$mu0, stringsAsFactors=F)
  region_DM_sig = data.frame('annotation'='significantly_DM', 'meth1'=sig$mu1, 'meth0'=sig$mu0, stringsAsFactors=F)

  df = rbind(avg_DM, avg_DM_sig)
  region_df = rbind(region_DM, region_DM_sig)

  diff_meth_plot = ggplot(df, aes(meth_diff)) + geom_histogram(binwidth=5, aes(y=..density..)) + facet_grid(. ~ annotation) + theme_bw()
  point_plot = ggplot(region_df, aes(meth1, meth0)) + geom_bin2d(binwidth=c(1,1)) + facet_grid(. ~ annotation)

  ggsave(filename='IDH2mut_v_NBM_mc_hmc_diff_meth.png', plot = diff_meth_plot, width = 8, height = 4, dpi = 300)
  ggsave(filename='IDH2mut_v_NBM_mc_hmc_diff_meth_point_heatmap.png', plot = point_plot, width = 8, height = 4, dpi = 300)

################################################################################
# Average % methylation change across annotations (methylSig)

  files = list.files(pattern='avg_meth_chg')

  annotations = gsub('IDH2mut_v_NBM_mc_hmc_avg_meth_chg_in_DM_','',files)
  annotations = gsub('.bed','',annotations)

  data = lapply(files, read.table, sep='\t', header=F, quote='', comment.char='', stringsAsFactors=F)
  names(data) = annotations

  df = Reduce(rbind, lapply(names(data), function(n){data.frame('annotation'=n, 'avg_diff_meth_in_DM'=data[[n]][,4], stringsAsFactors=F)}))

  avg_diff_meth_in_DM_plot = ggplot(df, aes(avg_diff_meth_in_DM)) + geom_histogram(binwidth=5, aes(y=..density..)) + facet_grid(. ~ annotation) + theme_bw()

  ggsave(filename='IDH2mut_v_NBM_mc_hmc_avg_diff_meth_in_DM_annotations.png', plot = avg_diff_meth_in_DM_plot, width = 8, height = 4, dpi = 300)
