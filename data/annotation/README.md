# Annotation Data
This file describes the data in `mint/data/annotation`.

## chromInfo_hg19.txt
These are chromosome lengths for hg19 taken directly from the UCSC Genome Browser. Datestamp 04-27-2009.

## cpg_islands_hg19_ucsc.txt
These are the CpG islands taken from the UCSC Table Browser. Group: Regulation, Track: CpG Islands. Header removed. Datestamp 10-24-2015.

```
#bin	chrom	chromStart	chromEnd	name	length	cpgNum	gcNum	perCpg	perGc	obsExp
```

## cpg_shores_hg19_ucsc.txt
These are the CpG shores, where a shore is defined as 1kb up and downstream of the CpG island start and end. The following one-liner was used to generate the shores. Datestamp 10-24-2015.

```{bash}
awk -v OFS='\t' 'NR > 1 {print $2, $3, $4}' ~/latte/mint/data/annotation/cpg_islands_hg19_ucsc.txt | bedtools flank -b 1000 -g ~/latte/mint/data/annotation/chromInfo_hg19.txt | sort -T . -k1,1 -k2,2n | bedtools merge > cpg_shores_hg19_ucsc.txt
```
