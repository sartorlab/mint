#!/bin/bash
set -e
set -u
set -o pipefail

# Paths from parameters
PATH_TO_BEDTOOLS=$1
PATH_TO_AWK=$2
CHROM_PATH=$3
outFile=$4
highMeth=$5
lowMeth=$6
noMethSig=$7
noMethNoSig=$8
peak=$9
noPeakSig=${10}
noPeakNoSig=${11}
ID=`basename $outFile _sample_classification.bed`

# tmp files
intTmp=classifications/sample/${ID}_tmpIntersect.txt
classTmp=classifications/sample/${ID}_tmpSampleClass.txt
joinTmp=classifications/sample/${ID}_tmpSampleJoin.txt

# tmp class files
tmp6101422=classifications/sample/${ID}_tmp6101422.txt
tmp12=classifications/sample/${ID}_tmp12.txt
tmp20=classifications/sample/${ID}_tmp20.txt
tmp18=classifications/sample/${ID}_tmp18.txt
tmp30=classifications/sample/${ID}_tmp30.txt
tmp284244=classifications/sample/${ID}_tmp284244.txt
tmp66=classifications/sample/${ID}_tmp66.txt

# Initial intersection
${PATH_TO_BEDTOOLS} multiinter \
	-header \
	-names highMeth lowMeth noMethSig noMethNoSig peak noPeakSig noPeakNoSig \
	-empty \
	-i ${highMeth} ${lowMeth} ${noMethSig} ${noMethNoSig} ${peak} ${noPeakSig} ${noPeakNoSig} \
	-g ${CHROM_PATH} \
> ${intTmp}

# Encode each region with the appropriate classification
# Columns 6-12 are the binary columns that should be operated on
${PATH_TO_AWK} -v OFS="\t" 'NR > 1 { \
	outFile = FILENAME; \
	group1 = ($6 * 3) + ($7 * 5) + ($8 * 7) + ($9 * 11); \
	group2 = ($10 * 2) + ($11 * 4) + ($12 * 6); \
	if (group1 * group2 == 6 || group1 * group2 == 10 || group1 * group2 == 14 || group1 * group2 == 22) { \
		sub(/tmpIntersect/, "tmp6101422", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 12) { \
		sub(/tmpIntersect/, "tmp12", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 20) { \
		sub(/tmpIntersect/, "tmp20", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 18) { \
		sub(/tmpIntersect/, "tmp18", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 30) { \
		sub(/tmpIntersect/, "tmp30", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 28 || group1 * group2 == 42 || group1 * group2 == 44) { \
		sub(/tmpIntersect/, "tmp284244", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} else if (group1 * group2 == 66) { \
		sub(/tmpIntersect/, "tmp66", outFile); \
		print $1, $2, $3, group1, group2, group1 * group2 > outFile; \
	} \
}' ${intTmp}

# Bedtools merge each of the files
tmpClassFiles=( ${tmp6101422} ${tmp12} ${tmp20} ${tmp18} ${tmp30} ${tmp284244} ${tmp66} )
for file in "${tmpClassFiles[@]}"
do
	bedtools merge -i <(sort -T . -k1,1 -k2,2n ${file}) -c 4,5,6 -o first,first,first > ${file}.merged
done

# Combine and sort in preparation for joining with the table
cat ${tmp6101422}.merged ${tmp12}.merged ${tmp20}.merged ${tmp18}.merged ${tmp30}.merged ${tmp284244}.merged ${tmp66}.merged | sort -T . -k1,1 -k2,2n > ${classTmp}

# Join on the classification code
# NOTE: Files must be sorted on the join field (code)
join \
	-1 2 \
	-2 6 \
	<(sort -T . -k2,2 ../../scripts/class_table_hybrid_sample.txt) \
	<(sort -T . -k6,6 ${classTmp}) \
> ${joinTmp}

# Create the bed file
${PATH_TO_AWK} -v OFS="\t" '{
	print $4, $5, $6, $2, "1000", ".", $5, $6, $3
}' ${joinTmp} | sort -T . -k1,1 -k2,2n > ${outFile}

# Cleanup
rm -f ${intTmp} ${classTmp} ${joinTmp} ${tmp6101422} ${tmp12} ${tmp20} ${tmp18} ${tmp30} ${tmp284244} ${tmp66} ${tmp6101422}.merged ${tmp12}.merged ${tmp20}.merged ${tmp18}.merged ${tmp30}.merged ${tmp284244}.merged ${tmp66}.merged
