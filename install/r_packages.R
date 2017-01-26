install.packages(c('devtools','optparse','readr','dplyr','ggplot2'), repos='http://cran.rstudio.com')
source('http://bioconductor.org/biocLite.R')
biocLite(c('GenomeInfoDb','GenomicRanges','IRanges','regioneR'))

devtools::install_github('rcavalcante/annotatr@v1.1.7')
devtools::install_github('sartorlab/methylSig@v0.4.4')
