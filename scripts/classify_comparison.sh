#!/bin/bash
set -e
set -u
set -o pipefail

# Steps required prior to classification

    # Example of preparation for hybrid experiment

        # sh ~/latte/Methylation/Methylation_Code/classification_comparison-wise_prepare_methylSig.sh -name IDH2vNBM

        # sh ~/latte/Methylation/Methylation_Code/classification_comparison-wise_prepare_pepr.sh -name IDH2vNBM -exp1 IDH2mut_1_affinity.bdg,IDH2mut_2_affinity.bdg,IDH2mut_3_affinity.bdg -exp2 NBM_1_affinity.bdg,NBM_2_affinity.bdg,NBM_3_affinity.bdg

    # Example of preparation for pure affinity experiment

        # sh ~/latte/Methylation/Methylation_Code/classification_comparison-wise_prepare_pepr.sh -name preeclamptic_v_normal_5hmc -exp1 preeclamptic_1_5hmc_affinity.bdg,preeclamptic_2_5hmc_affinity.bdg,preeclamptic_3_5hmc_affinity.bdg,preeclamptic_4_5hmc_affinity.bdg -exp2 normal_1_5hmc_affinity.bdg,normal_2_5hmc_affinity.bdg,normal_3_5hmc_affinity.bdg,normal_4_5hmc_affinity.bdg

        # sh ~/latte/Methylation/Methylation_Code/classification_comparison-wise_prepare_pepr.sh -name preeclamptic_v_normal_5mc -exp1 preeclamptic_1_5mc_affinity.bdg,preeclamptic_2_5mc_affinity.bdg,preeclamptic_3_5mc_affinity.bdg,preeclamptic_4_5mc_affinity.bdg -exp2 normal_1_5mc_affinity.bdg,normal_2_5mc_affinity.bdg,normal_3_5mc_affinity.bdg,normal_4_5mc_affinity.bdg

# Perform the setup to classification and then the classification_comparison-wise_regions.R will create the bed and bigBed
# the METHYLFILES and HYDROXYFILES are the two groups of four files created by the above preparation examples.

    # sh ~/latte/Methylation/Methylation_Code/classification_comparison-wise_regions.sh -name IDH2vNBM -methylFiles IDH2vNBM_methylSig_regions_DM_up.bed,IDH2vNBM_methylSig_regions_DM_down.bed,IDH2vNBM_methylSig_regions_noDM_signal.bed,IDH2vNBM_methylSig_regions_noDM_nosignal.bed -hydroxyFiles IDH2vNBM_PePr_up_peaks.bed,IDH2vNBM_PePr_down_peaks.bed,IDH2vNBM_PePr_noDM_signal.bed,IDH2vNBM_PePr_noDM_nosignal.bed

# Arguments
# -project The project name given to init_project.sh
# -comparison The comparison name used in process_bisulfite_comparison.sh
PROJECT=$2
COMPARISON=$4
METHYLFILES=$6
HYDROXYFILES=$8

# Go to the project directory
cd ~/latte/mint/${PROJECT}

# Argument expansion
methylFiles=$(echo $METHYLFILES | tr "," " ")
hydroxyFiles=$(echo $HYDROXYFILES | tr "," " ")

# Files
checkEqualTwo=./analysis/classification_comparison/${COMPARISON}_classification_test_not2.bed
checkDoubleM=./analysis/classification_comparison/${COMPARISON}_classification_test_doubleM.bed
checkDoubleH=./analysis/classification_comparison/${COMPARISON}_classification_test_doubleH.bed
multiInterFile=./analysis/classification_comparison/${COMPARISON}_classification_multiintersection.bed
firstClassFile=./analysis/classification_comparison/${COMPARISON}_classification_encoding.bed
finalClassFile=./analysis/classification_comparison/${COMPARISON}_classification_regions_final.bed
sortedFinalClassFile=./analysis/classification_comparison/${COMPARISON}_classification_regions_final_sorted.bed
classBB=./analysis/summary/ucsc_trackhub/hg19/${COMPARISON}_classification_regions.bb

    # This accomplishes the same partition and determines which file is responsible for the intersection!
    bedtools multiinter -header -names mDMup mDMdown mnDMs mnDMns hDMup hDMdown hnDMs hnDMns -empty -i  $methylFiles $hydroxyFiles -g ~/latte/Homo_sapiens/chromInfo_hg19.txt > $multiInterFile

    # Checks that all rows have one classification per 5mC + 5hmC and 5hmC each
        # Check that the number of intersections is always two. In other words, every region is always
        # categorized as DMup, DMdown, noDMsignal, or noDMnosignal with respect to 5mC (or 5mC + 5hmC) and 5hmC.
        awk -v OFS="\t" '$4 > 2 || $4 < 2 { print $1, $2, $3, $4, $5 }' $multiInterFile > $checkEqualTwo

        # Check that there is only one assignment in 5mC (or 5mC + 5hmC)
        awk -v OFS="\t" '$6 + $7 + $8 + $9 > 1 { print $1, $2, $3, $4, $5 }' $multiInterFile > $checkDoubleM

        # Check that there is only one assignment in 5hmC
        awk -v OFS="\t" '$10 + $11 + $12 + $13 > 1 { print $1, $2, $3, $4, $5 }' $multiInterFile > $checkDoubleH

    # Do the numerical encoding according to the classification table. The interpretation is the same whether
    # there is a hybrid experimental setup or a pure affinity experiment setup.
    awk -v OFS="\t" '{print $1, $2, $3, ($6 * 3) + ($7 * 5) + ($8 * 7) + ($9 * 11), ($10 * 2) + ($11 * 4) + ($12 * 6) + ($13 * 8), (($6 * 3) + ($7 * 5) + ($8 * 7) + ($9 * 11)) * (($10 * 2) + ($11 * 4) + ($12 * 6) + ($13 * 8))}' $multiInterFile > $firstClassFile

# R code to merge the encoded BED on the appropriate color scheme for the bigBed
    Rscript ~/latte/mint/scripts/classify_comparison.R --project=$PROJECT --comparison=$COMPARISON --inbed=$firstClassFile --outbed=$finalClassFile

# End R Code, now sort the result because R::merge messes it up
    sort -T . -k1,1 -k2,2n $finalClassFile > $sortedFinalClassFile
    mv $sortedFinalClassFile $finalClassFile

    # Create bigBed
    bedToBigBed $finalClassFile ~/latte/Homo_sapiens/chromInfo_hg19.txt $classBB

    # scp $classBB rcavalca@sftp.itd.umich.edu:/afs/umich.edu/user/r/c/rcavalca/Public/html/GSE52945
    # Create customtrack
    # track type=bigBed name=IDH2vNBM_class_regions description=IDH2vNBM_classification_regions db=hg19 itemRgb=on bigDataUrl=http://www-personal.umich.edu/~rcavalca/GSE52945/IDH2vNBM_classification_regions.bb
