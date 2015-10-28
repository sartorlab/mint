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
These are the CpG shores, where a shore is defined as 2kb up and downstream of the CpG island start and end. The following one-liner was used to generate the shores. Datestamp 10-28-2015.

```{bash}
bedtools subtract -a <(bedtools flank -b 2000 -i cpg_islands_hg19_ucsc.bed -g chromInfo_hg19.txt | sort -T . -k1,1 -k2,2n | bedtools merge) -b cpg_islands_hg19_ucsc.bed > cpg_shores_hg19_ucsc.bed
```

CpG islands and CpG shores are disjoint. The following intersection returns no lines:

```{bash}
bedtools intersect -a cpg_islands_hg19_ucsc.bed -b cpg_shores_hg19_ucsc.bed -wa -wb | head -10
```

## cpg_shelves_hg19_ucsc.bed
These are the CpG shelves, where a shelf is defined as the region between 2kb and 4kb away from a CpG island boundary. Datestamp 10-28-2015.

```{bash}
bedtools subtract -a <(bedtools flank -b 2000 -i cpg_shores_hg19_ucsc.bed -g chromInfo_hg19.txt | sort -T . -k1,1 -k2,2n | bedtools merge) -b <(cat cpg_islands_hg19_ucsc.bed cpg_shores_hg19_ucsc.bed | sort -T . -k1,1 -k2,2n) > cpg_shelves_hg19_ucsc.bed
```

CpG islands and CpG shores are disjoint from CpG shelves. The following intersections returns no lines:

```{bash}
bedtools intersect -a cpg_islands_hg19_ucsc.bed -b cpg_shelves_hg19_ucsc.bed -wa -wb | head -10
bedtools intersect -a cpg_shelves_hg19_ucsc.bed -b cpg_shores_hg19_ucsc.bed -wa -wb | head -10
```

## cpg_inter_hg19_ucsc.bed

These are regions that are not CpG islands, shores, nor shelves. **NOTE**: `bedtools complement` seems to complain about the third column of the chromosome sizes from UCSC which other bedtools subroutines don't complain about. The sort order of the chromosomes also seems to matter in this case, but not others. Datestamp 10-28-2015.

```{bash}
bedtools complement -i <(cat cpg_islands_hg19_ucsc.bed cpg_shores_hg19_ucsc.bed cpg_shelves_hg19_ucsc.bed | sort -T . -k1,1 -k2,2n) -g <(awk '{print $1 "\t" $2}' chromInfo_hg19.txt | sort -T . -k1,1 -k2,2n) > cpg_inter_hg19_ucsc.bed
```

CpG islands, shores, and shelves should be disjoint with the Inter CpG file.

```{bash}
bedtools intersect -a cpg_inter_hg19_ucsc.bed -b cpg_islands_hg19_ucsc.bed | head -10
bedtools intersect -a cpg_inter_hg19_ucsc.bed -b cpg_shores_hg19_ucsc.bed | head -10
bedtools intersect -a cpg_inter_hg19_ucsc.bed -b cpg_shelves_hg19_ucsc.bed | head -10
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

## windows_5kb_hg19.bed
The hg19 genome tiled into non-overlapping 5kb windows using the following [bedtools](http://bedtools.readthedocs.org) command:
```{bash}
bedtools makewindows -g chromInfo_hg19.txt -w 5000 > windows_5kb_hg19.bed
```
