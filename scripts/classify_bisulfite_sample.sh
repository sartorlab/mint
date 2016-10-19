#!/bin/bash
set -e
set -u
set -o pipefail

# Paths from parameters
PATH_TO_BEDTOOLS=$1
PATH_TO_AWK=$2
CHROM_PATH=$3
outFile=$4
mcHighMeth=$5
mcLowMeth=$6
mcNoMethSig=$7
mcNoMethNoSig=$8
hmcHighMeth=$9
hmcLowMeth=${10}
hmcNoMethSig=${11}
hmcNoMethNoSig=${12}
ID=`basename $outFile _sample_classification.bed`

# tmp files
intTmp=classifications/sample/${ID}_tmpIntersect.txt
classTmp=classifications/sample/${ID}_tmpSampleClass.txt
joinTmp=classifications/sample/${ID}_tmpSampleJoin.txt

# Initial intersection
${PATH_TO_BEDTOOLS} multiinter \
	-header \
	-names highMeth lowMeth noMethSig noMethNoSig peak noPeakSig noPeakNoSig \
	-empty \
	-i ${mcHighMeth} ${mcLowMeth} ${mcNoMethSig} ${mcNoMethNoSig} ${hmcHighMeth} ${hmcLowMeth} ${hmcNoMethSig} ${hmcNoMethNoSig} \
	-g ${CHROM_PATH} \
> ${intTmp}

# Encode each region with the appropriate classification
# Columns 6-12 are the binary columns that should be operated on
${PATH_TO_AWK} -v OFS="\t" 'NR > 1 { \
	group1 = ($6 * 3) + ($7 * 5) + ($8 * 7) + ($9 * 11); \
	group2 = ($10 * 2) + ($11 * 4) + ($12 * 6) + ($13 * 8); \
	if (group1 * group2 > 0) print $1, $2, $3, group1, group2, group1 * group2; \
}' ${intTmp} > ${classTmp}

# Join on the classification code
# NOTE: Files must be sorted on the join field (code)
join \
	-1 2 \
	-2 6 \
	<(sort -T . -k2,2 ../../scripts/class_table_bisulfite_sample.txt) \
	<(sort -T . -k6,6 ${classTmp}) \
> ${joinTmp}

# Create the bed file
${PATH_TO_AWK} -v OFS="\t" '{
	print $4, $5, $6, $2, "1000", ".", $5, $6, $3
}' ${joinTmp} | sort -T . -k1,1 -k2,2n > ${outFile}

# Cleanup
rm -f ${intTmp} ${classTmp} ${joinTmp}
