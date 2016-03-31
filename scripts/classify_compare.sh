#!/bin/bash
set -e
set -u
set -o pipefail

# Paths from parameters
PATH_TO_BEDTOOLS=$1
PATH_TO_AWK=$2
CHROM_PATH=$3
outFile=$4
mDMup=$5
mDMdown=$6
mNoDMSignal=$7
mNoDMNoSignal=$8
hDMup=$9
hDMdown=${10}
hNoDMSignal=${11}
hNoDMNoSignal=${12}
ID=`basename $outFile _compare_classification.bed`

# tmp files
intTmp=classifications/comparison/${ID}_tmpIntersect.txt
classTmp=classifications/comparison/${ID}_tmpSampleClass.txt
joinTmp=classifications/comparison/${ID}_tmpSampleJoin.txt

# Initial intersection
${PATH_TO_BEDTOOLS} multiinter \
	-header \
	-names mDMup mDMdown mNoDMSignal mNoDMNoSignal hDMup hDMdown hNoDMSignal hNoDMNoSignal \
	-empty \
	-i ${mDMup} ${mDMdown} ${mNoDMSignal} ${mNoDMNoSignal} ${hDMup} ${hDMdown} ${hNoDMSignal} ${hNoDMNoSignal} \
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
	<(sort -T . -k2,2 ../../scripts/class_table_compare.txt) \
	<(sort -T . -k6,6 ${classTmp}) \
> ${joinTmp}

# Create the bed file
${PATH_TO_AWK} -v OFS="\t" '{
	print $4, $5, $6, $2, "1000", ".", $5, $6, $3
}' ${joinTmp} | sort -T . -k1,1 -k2,2n > ${outFile}

# Cleanup
rm -f ${intTmp} ${classTmp} ${joinTmp}
