#!/bin/bash
set -e
set -u
set -o pipefail

peprUp=$1
peprDown=$2
peprCombined=$3

upPrefix=`basename ${peprUp} __PePr_up_peaks.bed`
downPrefix=`basename ${peprDown} __PePr_down_peaks.bed`

disUp=pulldown/pepr_peaks/${upPrefix}_disjoint_up.bed
disDown=pulldown/pepr_peaks/${upPrefix}_disjoint_down.bed 

bedops --difference \
	<(cut -f 1-3 ${peprUp}) \
	<(cut -f 1-3 ${peprDown}) \
| sort -T . -k1,1 -k2,2n \
> $disUp

bedops --difference \
	<(cut -f 1-3 ${peprDown}) \
	<(cut -f 1-3 ${peprUp}) \
| sort -T . -k1,1 -k2,2n \
> $disDown

cat \
  <(awk -v OFS="\t" '{ print $1, $2, $3, "hyper", "1000", ".", $2, $3, "0,0,255" }' ${disUp}) \
  <(awk -v OFS="\t" '{ print $1, $2, $3, "hypo", "1000", ".", $2, $3, "102,102,255" }' ${disDown}) \
| sort -T . -k1,1 -k2,2n > ${peprCombined}

# Clean up
rm -f ${disUp} ${disDown}
