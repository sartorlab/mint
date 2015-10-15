#!/bin/bash
set -e
set -u
set -o pipefail

# Arguments
# -project The project name giving to init_project.sh
# -mcPeaks The .narrowPeak file for mC pulldown from MACS2
# -mcZero The regions of zero input reads for mC from process_pulldown_sample.sh
# -hmcPeaks The .narrowPeak file for hmC pulldown from MACS2
# -hmcZero The regions of zero input reads for hmC from process_pulldown_sample.sh
# -humanID Some corresponding human readable ID
PROJECT=$2
mcPeaksID=$4
mcZeroID=$6
hmcPeaksID=$8
hmcZeroID=${10}
humanID=${12}

# Go to the project directory
cd ~/latte/mint/${PROJECT}

# Create appropriate file names
mcNarrowpeak=./analysis/macs_peaks/${mcPeaksID}_pulldown_macs2_peaks.narrowPeak
mcZero=./analysis/pulldown_coverages/${mcZeroID}_pulldown_zero.bdg
hmcNarrowpeak=./analysis/macs_peaks/${hmcPeaksID}_pulldown_macs2_peaks.narrowPeak
hmcZero=./analysis/pulldown_coverages/${hmcZeroID}_pulldown_zero.bdg
classBed=./analysis/classification_sample/${humanID}_sample_classification_pulldown.bed
classBB=./analysis/summary/${PROJECT}_hub/hg19/${humanID}_sample_classification_pulldown.bb

# Create tmp file names
# Need these dynamically named so we can classify samples in parallel without collisions
tmpCombinedPeaks=./analysis/classification_sample/tmp_${humanID}_combined.bed
tmpCombinedPeaksMerged=./analysis/classification_sample/tmp_${humanID}_combined_merged.bed
tmpIntersectPeaks=./analysis/classification_sample/tmp_${humanID}_intersect.bed
tmpMcPeaks=./analysis/classification_sample/tmp_${humanID}_mc.bed
tmpHmcPeaks=./analysis/classification_sample/tmp_${humanID}_hmc.bed
tmpNoPeaks=./analysis/classification_sample/tmp_${humanID}_nopeaks.bed
tmpMutuallyZero=./analysis/classification_sample/tmp_${humanID}_mutually_zero.bed
tmpNoPeaksSignal=./analysis/classification_sample/tmp_${humanID}_nopeaks_signal.bed
tmpNoPeaksNoSignal=./analysis/classification_sample/tmp_${humanID}_nopeaks_nosignal.bed
tmpMcAnnot=./analysis/classification_sample/tmp_${humanID}_mc_annot.bed
tmpHmcAnnot=./analysis/classification_sample/tmp_${humanID}_hmc_annot.bed
tmpIntersectAnnot=./analysis/classification_sample/tmp_${humanID}_intersect_annot.bed
tmpNoPSAnnot=./analysis/classification_sample/tmp_${humanID}_nopeaks_signal_annot.bed
tmpNoPNoSAnnot=./analysis/classification_sample/tmp_${humanID}_nopeaks_nosignal_annot.bed

# Combine the peaks
echo $humanID ': combining mC and hmC regions'
cat ${mcNarrowpeak} ${hmcNarrowpeak} | sort -T . -k1,1 -k2,2n > $tmpCombinedPeaks
bedtools merge -i $tmpCombinedPeaks > $tmpCombinedPeaksMerged

# Get the regions of overlap between the narrowPeak Files
# These are regions of mC + hmC according to pulldowns
echo $humanID ': intersecting mC and hmC regions'
bedtools intersect -a ${mcNarrowpeak} -b ${hmcNarrowpeak} | awk -v OFS='\t' '{print $1, $2, $3}' > $tmpIntersectPeaks

# Subtract these intersecting regions from each narrowPeak in turn
# These represent regions of only mc and only hmc
echo $humanID ': determining only mC and hmC regions'
bedtools subtract -a ${mcNarrowpeak} -b $tmpIntersectPeaks | awk -v OFS='\t' '{print $1, $2, $3}' > $tmpMcPeaks
bedtools subtract -a ${hmcNarrowpeak} -b $tmpIntersectPeaks | awk -v OFS='\t' '{print $1, $2, $3}' > $tmpHmcPeaks

# $tmpIntersectPeaks, $tmpMcPeaks, and $tmpHmcPeaks should be mutually disjoint

# Find regions where there were no peaks
echo $humanID ': determining no peak regions'
bedtools complement -i $tmpCombinedPeaksMerged -g ~/latte/Homo_sapiens/chromInfo_hg19.txt > $tmpNoPeaks

# Determine, among non-peak regions, where there is simply no peak, and where there is no signal
# If inputs are matched to pulldown, then require no signal to be intersection of both zero regions
if [ "${mcZero}" != "${hmcZero}" ]
then
  echo $humanID " Matched input files used..."
  bedtools intersect -a ${mcZero} -b ${hmcZero} | awk -v OFS='\t' '{print $1, $2, $3}' > $tmpMutuallyZero
else
  echo $humanID " Non-matched input file used..."
  cat ${mcZero} > $tmpMutuallyZero
fi

# $tmpNoPeaks - $tmpMutuallyZero gives no peaks with signal
echo $humanID ': determining no peak, but signal regions'
bedtools subtract -a $tmpNoPeaks -b $tmpMutuallyZero | awk -v OFS='\t' '{print $1, $2, $3}' > $tmpNoPeaksSignal

# $tmpNoPeaks intersect $tmpMutuallyZero gives no peaks without signal
echo $humanID ': determining no peak and no signal regions'
bedtools intersect -a $tmpNoPeaks -b $tmpMutuallyZero | awk -v OFS='\t' '{print $1, $2, $3}' > $tmpNoPeaksNoSignal

# Colors for the bigBed file
# $tmpMcPeaks - 5mC - red
# $tmpHmcPeaks - 5hmC - blue
# $tmpIntersectPeaks - 5mC and 5hmC - violet
# $tmpNoPeaksSignal - no 5mC nor 5hmC - black
# $tmpNoPeaksNoSignal - unclassifiable - gray

# Append the columns required for the track hub to each of the files listed above
echo $humanID ': annotating respective files'
awk -v OFS='\t' '{print $1, $2, $3, "mc:"$1":"$2, "1000", ".", $2, $3, "255,0,0"}' $tmpMcPeaks > $tmpMcAnnot
awk -v OFS='\t' '{print $1, $2, $3, "hmc:"$1":"$2, "1000", ".", $2, $3, "0,0,255"}' $tmpHmcPeaks > $tmpHmcAnnot
awk -v OFS='\t' '{print $1, $2, $3, "mc_hmc:"$1":"$2, "1000", ".", $2, $3, "102,0,204"}' $tmpIntersectPeaks > $tmpIntersectAnnot
awk -v OFS='\t' '{print $1, $2, $3, "no_peak:"$1":"$2, "1000", ".", $2, $3, "0,0,0"}' $tmpNoPeaksSignal > $tmpNoPSAnnot
awk -v OFS='\t' '{print $1, $2, $3, "unclassifiable:"$1":"$2, "1000", ".", $2, $3, "192,192,192"}' $tmpNoPeaksNoSignal > $tmpNoPNoSAnnot

# Combine, sort, and convert to bigBed
echo $humanID ': combining annotated pieces'
cat $tmpMcAnnot $tmpHmcAnnot $tmpIntersectAnnot $tmpNoPSAnnot $tmpNoPNoSAnnot | sort -T . -k1,1 -k2,2n > ${classBed}
echo $humanID ': creating bigBed file'
bedToBigBed ${classBed} ~/latte/Homo_sapiens/chromInfo_hg19.txt ${classBB}

# Remove tmp files created for the ${humanID}
echo $humanID ': removing temporary files'
find ./analysis/classification_sample/ -name "tmp_${humanID}*" | xargs rm
