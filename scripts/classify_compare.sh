#!/bin/bash
set -e
set -u
set -o pipefail

# Paths from parameters
CHROM_PATH=$1
outFile=$2
mDMup=$3
mDMdown=$4
mNoDMSignal=$5
mNoDMNoSignal=$6
hDMup=$7
hDMdown=$8
hNoDMSignal=$9
hNoDMNoSignal=${10}
ID=`basename $outFile _compare_classification.bed`

# tmp files
intTmp=classifications/comparison/${ID}_tmpIntersect.txt
classTmp=classifications/comparison/${ID}_tmpSampleClass.txt
joinTmp=classifications/comparison/${ID}_tmpSampleJoin.txt

# Initial intersection
bedtools multiinter \
	-header \
	-names mDMup mDMdown mNoDMSignal mNoDMNoSignal hDMup hDMdown hNoDMSignal hNoDMNoSignal \
	-empty \
	-i ${mDMup} ${mDMdown} ${mNoDMSignal} ${mNoDMNoSignal} ${hDMup} ${hDMdown} ${hNoDMSignal} ${hNoDMNoSignal} \
	-g ${CHROM_PATH} \
> ${intTmp}

# Encode each region with the appropriate classification
# Columns 6-12 are the binary columns that should be operated on
awk -v OFS="\t" 'NR > 1 { \
	group1 = ($6 * 3) + ($7 * 5) + ($8 * 7) + ($9 * 11); \
	group2 = ($10 * 2) + ($11 * 4) + ($12 * 6) + ($13 * 8); \
	if (group1 * group2 > 0) print $1, $2, $3, group1, group2, group1 * group2; \
}' ${intTmp} > ${classTmp}

# Join on the classification code
# NOTE: Files must be sorted on the join field (code)
join \
	-1 2 \
	-2 6 \
	<(sort -T . -k2,2n ../../scripts/class_table_compare.txt) \
	<(sort -T . -k6,6n ${classTmp}) \
> ${joinTmp}

# Create the bed file
awk -v OFS="\t" '{
	print $4, $5, $6, $2, "1000", ".", $5, $6, $3
}' ${joinTmp} | sort -T . -k1,1 -k2,2n > ${outFile}

# Cleanup
rm -f ${intTmp} ${classTmp} ${joinTmp}
