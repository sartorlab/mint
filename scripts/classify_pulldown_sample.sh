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

# tmp class files
tmp6=classifications/comparison/${ID}_tmp6.txt
tmp1014=classifications/comparison/${ID}_tmp1014.txt
tmp1218=classifications/comparison/${ID}_tmp1218.txt
tmp202830=classifications/comparison/${ID}_tmp202830.txt
tmp42=classifications/comparison/${ID}_tmp42.txt

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
	outFile = FILENAME; \
	group1 = ($6 * 3) + ($7 * 5) + ($8 * 7); \
	group2 = ($9 * 2) + ($10 * 4) + ($11 * 6); \
	if (group1 * group2 == 6) { \
		sub(/tmpIntersect/, "tmp6", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 10 || group1 * group2 == 14) { \
		sub(/tmpIntersect/, "tmp1014", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 12 || group1 * group2 == 18) { \
		sub(/tmpIntersect/, "tmp1218", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 20 || group1 * group2 == 28 || group1 * group2 == 30) { \
		sub(/tmpIntersect/, "tmp202830", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 42) { \
		sub(/tmpIntersect/, "tmp42", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} \
}' ${intTmp}

# Bedtools merge each of the files
tmpClassFiles=( ${tmp6} ${tmp1014} ${tmp1218} ${tmp202830} ${tmp42} )
for file in "${tmpClassFiles[@]}"
do
	bedtools merge -i <(sort -T . -k1,1 -k2,2n ${file}) -c 4,5,6 -o first,first,first > ${file}.merged
done

# Combine and sort in preparation for joining with the table
cat ${tmp6}.merged ${tmp1014}.merged ${tmp1218}.merged ${tmp202830}.merged ${tmp42}.merged | sort -T . -k1,1 -k2,2n > ${classTmp}

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
rm -f ${intTmp} ${classTmp} ${joinTmp} ${tmp6} ${tmp1014} ${tmp1218} ${tmp202830} ${tmp42} ${tmp6}.merged ${tmp1014}.merged ${tmp1218}.merged ${tmp202830}.merged ${tmp42}.merged
