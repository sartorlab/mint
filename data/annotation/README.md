# Annotation Data
This file describes the data in `mint/data/annotation`.

## andersson_permissive_enhancers.bed
These are the ~43,000 permissive enhancers determined by [Andersson, et al](http://www.nature.com/nature/journal/v507/n7493/abs/nature12787.html). Data were downloaded directly from the pre-defined enhancer sets of the [Transcribed Enhancer Atlas](http://enhancer.binf.ku.dk/presets/). The direct URL for the data is [http://enhancer.binf.ku.dk/presets/robust_enhancers.bed](http://enhancer.binf.ku.dk/presets/robust_enhancers.bed).

## andersson_permissive_enhancers_less_5kb_loci.bed
Some of the Andersson enhancers actually occur within the 5kb locus definitions from `chipenrich`. Consequently we need to do the following:

```{bash}
bedtools subtract -a andersson_permissive_enhancers.bed -b annot_ldefs_5kb_hg19.bed | awk -v OFS='\t' '{print $1, $2, $3, "enhancers"}' > annot_enhancers_hg19.bed
```

## chromInfo_hg19.txt
These are chromosome lengths for hg19 taken directly from the UCSC Genome Browser. Datestamp 04-27-2009.

## cpg_islands_hg19_ucsc.txt
These are the CpG islands taken from the UCSC Table Browser. Group: Regulation, Track: CpG Islands. Header removed. Datestamp 10-24-2015.

```
#bin	chrom	chromStart	chromEnd	name	length	cpgNum	gcNum	perCpg	perGc	obsExp
```

## annot_cpg_islands_hg19.bed
Taken from `cpg_islands_hg19_ucsc.txt` but modified with the awk one-liner:

```{bash}
awk -v OFS='\t' '{print $2, $3, $4}' cpg_islands_hg19_ucsc.txt | sort -T . -k1,1 -k2,2n | awk -v OFS='\t' '{print $1, $2, $3, "cpg_islands"}' > annot_cpg_islands_hg19.bed
```

## annot_cpg_shores_hg19.bed
These are the CpG shores, where a shore is defined as 2kb up and downstream of the CpG island start and end. The following one-liner was used to generate the shores. Datestamp 10-28-2015.

```{bash}
bedtools subtract -a <(bedtools flank -b 2000 -i annot_cpg_islands_hg19.bed -g chromInfo_hg19.txt | sort -T . -k1,1 -k2,2n | bedtools merge) -b annot_cpg_islands_hg19.bed | awk -v OFS='\t' '{print $1, $2, $3, "cpg_shores"}' > annot_cpg_shores_hg19.bed
```

CpG islands and CpG shores are disjoint. The following intersection returns no lines:

```{bash}
bedtools intersect -a annot_cpg_islands_hg19.bed -b annot_cpg_shores_hg19.bed -wa -wb | head -10
```

## annot_cpg_shelves_hg19.bed
These are the CpG shelves, where a shelf is defined as the region between 2kb and 4kb away from a CpG island boundary. Datestamp 10-28-2015.

```{bash}
bedtools subtract -a <(bedtools flank -b 2000 -i annot_cpg_shores_hg19.bed -g chromInfo_hg19.txt | sort -T . -k1,1 -k2,2n | bedtools merge) -b <(cat annot_cpg_islands_hg19.bed annot_cpg_shores_hg19.bed | sort -T . -k1,1 -k2,2n) | awk -v OFS='\t' '{print $1, $2, $3, "cpg_shelves"}' > annot_cpg_shelves_hg19.bed
```

CpG islands and CpG shores are disjoint from CpG shelves. The following intersections returns no lines:

```{bash}
bedtools intersect -a annot_cpg_islands_hg19.bed -b annot_cpg_shelves_hg19.bed -wa -wb | head -10
bedtools intersect -a annot_cpg_shelves_hg19.bed -b annot_cpg_shores_hg19.bed -wa -wb | head -10
```

## annot_cpg_inter_hg19.bed

These are regions that are not CpG islands, shores, nor shelves. **NOTE**: `bedtools complement` seems to complain about the third column of the chromosome sizes from UCSC which other bedtools subroutines don't complain about. The sort order of the chromosomes also seems to matter in this case, but not others. **NOTE**: Some lines in the result have start 0 and end 0 which prevents the proper use of the file downstream. Hence the last awk pipe. Datestamp 10-28-2015.

```{bash}
bedtools complement -i <(cat annot_cpg_islands_hg19.bed annot_cpg_shores_hg19.bed annot_cpg_shelves_hg19.bed | sort -T . -k1,1 -k2,2n | bedtools merge) -g <(awk '{print $1 "\t" $2}' chromInfo_hg19.txt | sort -T . -k1,1 -k2,2n) | awk -v OFS='\t' '$2 != $3 {print $1, $2, $3, "cpg_inter"}' > annot_cpg_inter_hg19.bed
```

CpG islands, shores, and shelves should be disjoint with the Inter CpG file.

```{bash}
bedtools intersect -a annot_cpg_inter_hg19.bed -b annot_cpg_islands_hg19.bed | head -10
bedtools intersect -a annot_cpg_inter_hg19.bed -b annot_cpg_shores_hg19.bed | head -10
bedtools intersect -a annot_cpg_inter_hg19.bed -b annot_cpg_shelves_hg19.bed | head -10
```

## ldef_5kb_hg19_reduced.bed
These are the 5kb locus definitions from `chipenrich.data` run through the `mint/data/scripts/reduce_ldef.R` script. The result is that each locus should correspond to one line, except in the cases where a locus may be interrupted by another on the opposite strand. They will be used as the promoter regions.

To get original 5kb locus definitions from `chipenrich.data`:
```{r}
library(chipenrich.data)
data('locusdef.hg19.5kb')
write.table(locusdef.hg19.5kb@dframe, file='ldef_5kb_hg19.bed', sep='\t', col.names=F, row.names=F, quote=F)
```
And to get the reduced locus definitions:
```{bash}
Rscript ../scripts/reduce_ldef.R --inpath ldef_5kb_hg19.bed --outpath annot_ldefs_5kb_hg19_tmp.bed
awk -v OFS='\t' '{print $1, $2, $3, "promoter:"$4}' annot_ldefs_5kb_hg19_tmp.bed > annot_ldefs_5kb_hg19.bed
rm annot_ldefs_5kb_hg19_tmp.bed
```

## windows_5kb_hg19.bed
The hg19 genome tiled into non-overlapping 5kb windows using the following [bedtools](http://bedtools.readthedocs.org) command:
```{bash}
bedtools makewindows -g chromInfo_hg19.txt -w 5000 | awk -v OFS='\t' '{print $1, $2, $3, "windows_5kb"}' > annot_windows_5kb_hg19.bed
```
