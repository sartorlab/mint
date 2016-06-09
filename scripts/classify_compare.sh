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

# tmp class files
tmp6=classifications/comparison/${ID}_tmp6.txt
tmp10=classifications/comparison/${ID}_tmp10.txt
tmp1422=classifications/comparison/${ID}_tmp1422.txt
tmp12=classifications/comparison/${ID}_tmp12.txt
tmp20=classifications/comparison/${ID}_tmp20.txt
tmp2844=classifications/comparison/${ID}_tmp2844.txt
tmp1824=classifications/comparison/${ID}_tmp1824.txt
tmp3040=classifications/comparison/${ID}_tmp3040.txt
tmp426656=classifications/comparison/${ID}_tmp426656.txt
tmp88=classifications/comparison/${ID}_tmp88.txt

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
	outFile = FILENAME; \
	group1 = ($6 * 3) + ($7 * 5) + ($8 * 7) + ($9 * 11); \
	group2 = ($10 * 2) + ($11 * 4) + ($12 * 6) + ($13 * 8); \
	if (group1 * group2 == 6) { \
		sub(/tmpIntersect/, "tmp6", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 10) { \
		sub(/tmpIntersect/, "tmp10", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 14 || group1 * group2 == 22) { \
		sub(/tmpIntersect/, "tmp1422", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 12) { \
		sub(/tmpIntersect/, "tmp12", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 20) { \
		sub(/tmpIntersect/, "tmp20", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 28 || group1 * group2 == 44) { \
		sub(/tmpIntersect/, "tmp2844", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 18 || group1 * group2 == 24) { \
		sub(/tmpIntersect/, "tmp1824", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 30 || group1 * group2 == 40) { \
		sub(/tmpIntersect/, "tmp3040", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 42 || group1 * group2 == 66 || group1 * group2 == 56) { \
		sub(/tmpIntersect/, "tmp426656", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 88) { \
		sub(/tmpIntersect/, "tmp88", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} \
}' ${intTmp}

# Bedtools merge each of the files
tmpClassFiles=( ${tmp6} ${tmp10} ${tmp1422} ${tmp12} ${tmp20} ${tmp2844} ${tmp1824} ${tmp3040} ${tmp426656} ${tmp88} )
for file in "${tmpClassFiles[@]}"
do
	bedtools merge -i <(sort -T . -k1,1 -k2,2n ${file}) -c 4,5,6 -o first,first,first > ${file}.merged
done

# Combine and sort in preparation for joining with the table
cat ${tmp6}.merged ${tmp10}.merged ${tmp1422}.merged ${tmp12}.merged ${tmp20}.merged ${tmp2844}.merged ${tmp1824}.merged ${tmp3040}.merged ${tmp426656}.merged ${tmp88}.merged | sort -T . -k1,1 -k2,2n > ${classTmp}

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
rm -f ${intTmp} ${classTmp} ${joinTmp} ${tmp6} ${tmp10} ${tmp1422} ${tmp12} ${tmp20} ${tmp2844} ${tmp1824} ${tmp3040} ${tmp426656} ${tmp88} ${tmp6}.merged ${tmp10}.merged ${tmp1422}.merged ${tmp12}.merged ${tmp20}.merged ${tmp2844}.merged ${tmp1824}.merged ${tmp3040}.merged ${tmp426656}.merged ${tmp88}.merged
