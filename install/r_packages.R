install.packages(c('devtools','optparse'), repos='http://cran.rstudio.com')

devtools::install_github('sartorlab/methylSig@v0.4.4')

source("http://bioconductor.org/biocLite.R")
biocLite(c("annotatr","csaw"))
