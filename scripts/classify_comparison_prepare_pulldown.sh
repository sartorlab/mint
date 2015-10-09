#!/bin/bash
set -e
set -u
set -o pipefail

# Example usage
# sh classification_comparison_prepare_pulldown.sh -project GSE52945 -comparison IDH2vNBM -exp1cov IDH2mut_1_affinity.bdg,IDH2mut_2_affinity.bdg,IDH2mut_3_affinity.bdg -exp2cov NBM_1_affinity.bdg,NBM_2_affinity.bdg,NBM_3_affinity.bdg

# sh classification_comparison_prepare_pulldown.sh -project GSE63743 -comparison preeclamptic_v_normal_5hmc -exp1cov preeclamptic_1_5hmc_affinity.bdg,preeclamptic_2_5hmc_affinity.bdg,preeclamptic_3_5hmc_affinity.bdg,preeclamptic_4_5hmc_affinity.bdg -exp2cov normal_1_5hmc_affinity.bdg,normal_2_5hmc_affinity.bdg,normal_3_5hmc_affinity.bdg,normal_4_5hmc_affinity.bdg

# sh classification_comparison_prepare_pulldown.sh -project GSE63743 -comparison preeclamptic_v_normal_5mc -exp1cov preeclamptic_1_5mc_affinity.bdg,preeclamptic_2_5mc_affinity.bdg,preeclamptic_3_5mc_affinity.bdg,preeclamptic_4_5mc_affinity.bdg -exp2cov normal_1_5mc_affinity.bdg,normal_2_5mc_affinity.bdg,normal_3_5mc_affinity.bdg,normal_4_5mc_affinity.bdg

# Arguments
# -project The project name given to init_project.sh
# -comparison The comparison name used in process_bisulfite_comparison.sh
# -exp1cov The coverage files for group 1 (must be same same order as PePr call)
# -exp2cov The coverage files for group 2 (must be same same order as PePr call)
PROJECT=$2
COMPARISON=$4
EXP1COV=$6
EXP2COV=$8

# Go to the project directory
cd ~/latte/mint/${PROJECT}

# Split PePr output into up, down, no DM / signal, and no DM / no signal
# Also need to give bedGraphs of read coverage to determine signal

exp1=$(echo $EXP1COV | tr "," " ")
exp2=$(echo $EXP2COV | tr "," " ")

ogUpPeaks=./analysis/pepr_peaks/${COMPARISON}__PePr_up_peaks.bed
ogDownPeaks=./analysis/pepr_peaks/${COMPARISON}__PePr_down_peaks.bed

upSorted=./analysis/pepr_peaks/${COMPARISON}_PePr_up_peaks_sorted.bed
downSorted=./analysis/pepr_peaks/${COMPARISON}_PePr_down_peaks_sorted.bed

upDisjoint=./analysis/pepr_peaks/${COMPARISON}_PePr_up_peaks_disjoint.bed
downDisjoint=./analysis/pepr_peaks/${COMPARISON}_PePr_down_peaks_disjoint.bed

upPeaks=./analysis/pepr_peaks/${COMPARISON}_PePr_up_peaks.bed
downPeaks=./analysis/pepr_peaks/${COMPARISON}_PePr_down_peaks.bed

allPeaks=./analysis/pepr_peaks/${COMPARISON}_PePr_peaks.bed
allSorted=./analysis/pepr_peaks/${COMPARISON}_PePr_peaks_sorted.bed

allPeaksNoDM=./analysis/pepr_peaks/${COMPARISON}_PePr_noDM.bed

noDMSignalPeaks=./analysis/pepr_peaks/${COMPARISON}_PePr_noDM_signal.bed
noDMNoSignalPeaks=./analysis/pepr_peaks/${COMPARISON}_PePr_noDM_nosignal.bed

noDMSignalPeaksSorted=./analysis/pepr_peaks/${COMPARISON}_PePr_noDM_signal_sorted.bed
noDMNoSignalPeaksSorted=./analysis/pepr_peaks/${COMPARISON}_PePr_noDM_nosignal_sorted.bed

signal=./analysis/pepr_peaks/${COMPARISON}_signal_pulldown.bed
noSignal=./analysis/pepr_peaks/${COMPARISON}_nosignal_pulldown.bed

# Split PePr results into up, down, no DM / signal, and no DM / no signal
# These 5hmC files are the corresponding rows of the classification table

    # DM up
    awk '{ print $1 "\t" $2 "\t" $3 }' $ogUpPeaks > $upPeaks
    sort -T . -k1,1 -k2,2n $upPeaks > $upSorted
    mv $upSorted $upPeaks

    # DM down
    awk '{ print $1 "\t" $2 "\t" $3 }' $ogDownPeaks > $downPeaks
    sort -T . -k1,1 -k2,2n $downPeaks > $downSorted
    mv $downSorted $downPeaks

    # Make PePr up and down peaks disjoint
    bedops --difference $upPeaks $downPeaks > $upDisjoint
    bedops --difference $downPeaks $upPeaks > $downDisjoint

        # This should have 0 length
        bedtools intersect -a $upDisjoint -b $downDisjoint | wc -l

    mv $upDisjoint $upPeaks
    mv $downDisjoint $downPeaks

    # Combined DM are needed for the complement (noDM)
    cat $upPeaks $downPeaks > $allPeaks
    sort -T . -k1,1 -k2,2n $allPeaks > $allSorted
    mv $allSorted $allPeaks

    # No DM needs to subsequently be split into signal and no signal
    bedtools complement -i $allPeaks -g ~/latte/Homo_sapiens/chromInfo_hg19.txt > $allPeaksNoDM

    # Determine regions of signal and no signal
        # Merge all affinity signals into a single signal bed
        bedops --merge $exp1 $exp2 > $signal

        # No signal comprises the complement of the signal
        bedtools complement -i $signal -g ~/latte/Homo_sapiens/chromInfo_hg19.txt > $noSignal

    # No DM and signal
    bedtools intersect -a $signal -b $allPeaksNoDM > $noDMSignalPeaks
    sort -T . -k1,1 -k2,2n $noDMSignalPeaks > $noDMSignalPeaksSorted
    mv $noDMSignalPeaksSorted $noDMSignalPeaks

    # No DM and no signal
    bedtools intersect -a $noSignal -b $allPeaksNoDM > $noDMNoSignalPeaks
    sort -T . -k1,1 -k2,2n $noDMNoSignalPeaks > $noDMNoSignalPeaksSorted
    mv $noDMNoSignalPeaksSorted $noDMNoSignalPeaks

        rm $signal
        rm $noSignal

    # Check these four files are mutually disjoint
        # *_PePr_up_peaks.bed
        # *_PePr_down_peaks.bed
        # *_PePr_noDM_signal.bed
        # *_PePr_noDM_nosignal.bed
        echo 'All of the following values should be 0!'
        echo 'Up peaks disjoint from down peaks?'
        bedtools intersect -a $upPeaks -b $downPeaks | wc -l
        echo 'Up peaks distinct from no DM and signal regions?'
        bedtools intersect -a $upPeaks -b $noDMSignalPeaks | wc -l
        echo 'Up peaks distinct from no DM and no signal regions?'
        bedtools intersect -a $upPeaks -b $noDMNoSignalPeaks | wc -l
        echo 'Down peaks distinct from no DM and signal regions?'
        bedtools intersect -a $downPeaks -b $noDMSignalPeaks | wc -l
        echo 'Down peaks distinct from no DM and no signal regions?'
        bedtools intersect -a $downPeaks -b $noDMNoSignalPeaks | wc -l
        echo 'No DM signal regions distinct from no DM no signal regions?'
        bedtools intersect -a $noDMSignalPeaks -b $noDMNoSignalPeaks | wc -l
