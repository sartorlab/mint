#!/bin/bash
set -e
set -u
set -o pipefail

# Paths from parameters
PATH_TO_BEDTOOLS=$1
PATH_TO_AWK=$2
CHROM_PATH=$3
outFile=$4
mcPeak=$5
mcNoPeakSig=$6
mcNoPeakNoSig=$7
hmcPeak=$8
hmcNoPeakSig=$9
hmcNoPeakNoSig=${10}
ID=`basename $outFile _sample_classification.bed`

# tmp files
intTmp=classifications/sample/${ID}_tmpIntersect.txt
classTmp=classifications/sample/${ID}_tmpSampleClass.txt
joinTmp=classifications/sample/${ID}_tmpSampleJoin.txt

# Initial intersection
${PATH_TO_BEDTOOLS} multiinter \
	-header \
	-names mcPeak mcNoPeakSig mcNoPeakNoSig hmcPeak hmcNoPeakSig hmcNoPeakNoSig \
	-empty \
	-i ${mcPeak} ${mcNoPeakSig} ${mcNoPeakNoSig} ${hmcPeak} ${hmcNoPeakSig} ${hmcNoPeakNoSig} \
	-g ${CHROM_PATH} \
> ${intTmp}

# Encode each region with the appropriate classification
# Columns 6-11 are the binary columns that should be operated on
${PATH_TO_AWK} -v OFS="\t" 'NR > 1 { \
	group1 = ($6 * 3) + ($7 * 5) + ($8 * 7); \
	group2 = ($9 * 2) + ($10 * 4) + ($11 * 6); \
	if (group1 * group2 > 0) print $1, $2, $3, group1, group2, group1 * group2; \
}' ${intTmp} > ${classTmp}

# Join on the classification code
# NOTE: Files must be sorted on the join field (code)
join \
	-1 2 \
	-2 6 \
	<(sort -T . -k2,2n ../../scripts/class_table_pulldown_sample.txt) \
	<(sort -T . -k6,6n ${classTmp}) \
> ${joinTmp}

# Create the bed file
${PATH_TO_AWK} -v OFS="\t" '{
	print $4, $5, $6, $2, "1000", ".", $5, $6, $3
}' ${joinTmp} | sort -T . -k1,1 -k2,2n > ${outFile}

# Cleanup
rm -f ${intTmp} ${classTmp} ${joinTmp}
