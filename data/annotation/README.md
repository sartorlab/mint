# Annotation Data
This file describes the data in `mint/data/annotation`.

## andersson_permissive_enhancers.bed
These are the ~43,000 permissive enhancers determined by [Andersson, et al](http://www.nature.com/nature/journal/v507/n7493/abs/nature12787.html). Data were downloaded directly from the pre-defined enhancer sets of the [Transcribed Enhancer Atlas](http://enhancer.binf.ku.dk/presets/). The direct URL for the data is [http://enhancer.binf.ku.dk/presets/robust_enhancers.bed](http://enhancer.binf.ku.dk/presets/robust_enhancers.bed).

## chromInfo_hg19.txt
These are chromosome lengths for hg19 taken directly from the UCSC Genome Browser. Datestamp 04-27-2009.

## cpg_islands_hg19_ucsc.txt
These are the CpG islands taken from the UCSC Table Browser. Group: Regulation, Track: CpG Islands. Header removed. Datestamp 10-24-2015.

```
#bin	chrom	chromStart	chromEnd	name	length	cpgNum	gcNum	perCpg	perGc	obsExp
```

## cpg_islands_hg19_ucsc.bed
Taken from `cpg_islands_hg19_ucsc.txt` but modified with the awk one-liner:

```{bash}
awk -v OFS='\t' '{print $2, $3, $4}' cpg_islands_hg19_ucsc.txt | sort -T . -k1,1 -k2,2n > cpg_islands_hg19_ucsc.bed
```

## cpg_shores_hg19_ucsc.bed
These are the CpG shores, where a shore is defined as 1kb up and downstream of the CpG island start and end. The following one-liner was used to generate the shores. Datestamp 10-24-2015.

```{bash}
awk -v OFS='\t' 'NR > 1 {print $2, $3, $4}' ~/latte/mint/data/annotation/cpg_islands_hg19_ucsc.txt | bedtools flank -b 1000 -g ~/latte/mint/data/annotation/chromInfo_hg19.txt | sort -T . -k1,1 -k2,2n | bedtools merge > cpg_shores_hg19_ucsc.bed
```

## ldef_5kb_hg19_reduced.bed
These are the 5kb locus definitions from `chipenrich.data` run through the `mint/data/scripts/reduce_ldef.R` script. The result is that each locus should correspond to one line, except in the cases where a locus may be interrupted by another on the opposite strand. They will be used as the promoter regions.

To get original 5kb locus definitions from `chipenrich.data`:
```{r}
library(chipenrich.data)
data('locusdef.hg19.5kb')
write.table(locusdef.hg19.5kb@dframe, file='path/to/ldef_5kb_hg19.bed', sep='\t', col.names=F, row.names=F, quote=F)
```
And to get the reduced locus definitions:
```{bash}
Rscript reduce_ldef.R --inpath path/to/ldef_5kb_hg19.bed --outpath path/to/ldef_5kb_hg19_reduced.bed
```
