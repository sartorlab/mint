### Versions

This file gives the current version of software used in the pipeline as of May 15, 2016.

* [`bedtools` v2.25.0](https://github.com/arq5x/bedtools2/releases/tag/v2.25.0)
* [`bedops` v2.4.14](https://github.com/bedops/bedops/releases/tag/v2.4.14)
* [`bismark` v0.16.1](https://github.com/FelixKrueger/Bismark/releases/tag/0.16.1)
* [`bowtie2` v2.2.4](https://github.com/BenLangmead/bowtie2/releases/tag/v2.2.4)
* [`cutadapt` v1.9.1](https://pypi.python.org/pypi/cutadapt/1.9.1)
* [`FastQC` v0.11.5](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5_source.zip)
* [`macs2` v2.1.0.20140616](https://pypi.python.org/pypi/MACS2/2.1.0.20140616)
* [`multiqc` v0.6.0](https://github.com/ewels/MultiQC/releases/tag/v0.6)
* [`PePr` v1.1.9](https://github.com/shawnzhangyx/PePr/releases/tag/1.1.9)
* [`R` >= v3.2.5](https://cran.r-project.org)
	* [`annotatr` v0.7.3](https://github.com/rcavalcante/annotatr/releases/tag/v0.7.3)
	* `devtools`
	* `dplyr`
	* `ggplot2`
	* [`methylSig` v0.4.3](https://github.com/sartorlab/methylSig/releases/tag/v0.4.3)
	* `optparse`
	* `readr`
* [`samtools` v0.1.19](https://github.com/samtools/samtools/releases/tag/0.1.19)
* [`trim_galore` v0.4.1](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip)

The following R command will install the required R packages (assuming fresh R install):

```{r}
# Install CRAN packages
install.packages(c('devtools','optparse','readr','dplyr','ggplot2'), repos='http://cran.rstudio.com')

# Install Bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocStyle","GenomeInfoDb","IRanges","GenomicRanges"))

# Install GitHub packages
install_github('rcavalcante/annotatr@v0.7.3')
install_github('sartorlab/methylSig@v0.4.3')
```
