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

  summary_files = list.files('analysis/summary/tables', full.names=T)
  summary_files = grep(humanID, summary_files, value=T)

  files = list(
    cpg = grep('cpg', summary_files, value=T),
    windows_5kb = grep('windows', summary_files, value=T),
    prom_enh = grep('(promoters)|(enhancers)', summary_files, value=T))

################################################################################
# Average % methylation and coverage across CpG annotations and promoter/enhancers

  lapply(names(files), function(annot){
    message(sprintf('Plotting average average coverage across %s annotations in %s', annot, humanID))

    files_cov = grep('avg_coverage', files[[annot]], value=T)

    annotations = gsub(sprintf('analysis/summary/tables/%s_avg_coverage_',humanID),'',files_cov)
    annotations = gsub('.bed','',annotations)

    data_cov = lapply(files_cov, read.table, sep='\t', header=F, quote='', comment.char='', stringsAsFactors=F)
    names(data_cov) = annotations

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

    png_cov = sprintf('analysis/summary/figures/%s_pulldown_%s_avg_cov.png', humanID, annot)
    plot_cov =
      ggplot(df_cov, aes(avg_coverage)) +
      geom_histogram(aes(y=..density..)) +
      facet_grid(. ~ annotation) +
      ggtitle('Avg. Cov. Across Annotations') +
      theme_bw()
    ggsave(filename = png_cov, plot = plot_cov, width = width, height = height, dpi = 300)

  })
