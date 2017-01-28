#!/bin/bash
set -e
set -u
set -o pipefail

peprChip1=$1
peprChip2=$2
peprCombined=$3
chip1Name=$4
chip2Name=$5

upPrefix=`basename ${peprChip1} __PePr_chip1_peaks.bed`
downPrefix=`basename ${peprChip2} __PePr_chip2_peaks.bed`

disChip1=pulldown/pepr_peaks/${upPrefix}_disjoint_chip1.bed
disChip2=pulldown/pepr_peaks/${upPrefix}_disjoint_chip2.bed

bedops --difference \
	<(cut -f 1-3 ${peprChip1}) \
	<(cut -f 1-3 ${peprChip2}) \
| sort -T . -k1,1 -k2,2n \
> $disChip1

bedops --difference \
	<(cut -f 1-3 ${peprChip2}) \
	<(cut -f 1-3 ${peprChip1}) \
| sort -T . -k1,1 -k2,2n \
> $disChip2

cat \
	<(awk -v OFS="\t" -v CHIP1=${chip1Name} '{ print $1, $2, $3, CHIP1, "1000", ".", $2, $3, "0,0,255" }' ${disChip1}) \
	<(awk -v OFS="\t" -v CHIP2=${chip2Name} '{ print $1, $2, $3, CHIP2, "1000", ".", $2, $3, "102,102,255" }' ${disChip2}) \
| sort -T . -k1,1 -k2,2n > ${peprCombined}

# Clean up
rm -f ${disChip1} ${disChip2}
