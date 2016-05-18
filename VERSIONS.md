### Versions

This file gives the current version of software used in the pipeline as of May 15, 2016.

* bedtools v2.25.0
* bedops v2.4.14
* bismark v0.16.1
* bowtie2 v2.2.4
* macs2 v2.1.0.20140616
* PePr v1.1.3
* R v3.2.5
	* annotatr v0.7.1
	* methylSig v0.4.3
* samtools v0.1.19
* trim_galore v0.4.1
* cutadapt v1.9.1
* FastQC v0.11.5

The following R command will install the required R packages (assuming fresh R install):

```{r}
install.packages(c('readr','optparse','ggplot2','dplyr','devtools'), repos='http://cran.rstudio.com')
source("http://bioconductor.org/biocLite.R")
biocLite(c("BiocStyle","GenomeInfoDb","IRanges","GenomicRanges"))
devtools::install_github('sartorlab/methylSig')
devtools::install_github('rcavalcante/annotatr')
```
