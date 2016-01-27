#!/bin/bash
set -e
set -u
set -o pipefail

# Arguments
# -project The project name given to init_project.sh
# -comparison The comparison name used in process_bisulfite_comparison.sh
PROJECT=$2
COMPARISON=$4

# Go to the project directory
cd ~/latte/mint/${PROJECT}

# Split methylSig results into up, down, no DM / signal, and no DM / no signal
# we are using the 5th (pvalue) column instead of the 6th (qvalue) column

methylSigFile=./analysis/methylsig_calls/${COMPARISON}.txt
methylSigSmall=./analysis/methylsig_calls/${COMPARISON}.bed

dmUp=./analysis/methylsig_calls/${COMPARISON}_methylSig_DM_up.bed
dmDown=./analysis/methylsig_calls/${COMPARISON}_methylSig_DM_down.bed

dmAll=./analysis/methylsig_calls/${COMPARISON}_methylSig_DM.bed
dmAllSorted=./analysis/methylsig_calls/${COMPARISON}_methylSig_DM_sorted.bed

noDMSignal=./analysis/methylsig_calls/${COMPARISON}_methylSig_noDM_signal.bed
noDMNoSignal=./analysis/methylsig_calls/${COMPARISON}_methylSig_noDM_nosignal.bed

noDMSignalSorted=./analysis/methylsig_calls/${COMPARISON}_methylSig_noDM_signal_sorted.bed
noDMNoSignalSorted=./analysis/methylsig_calls/${COMPARISON}_methylSig_noDM_nosignal_sorted.bed

    # NOTE: For DM up and down, it may be desirable to require a minimum methylation differential
    # DM up
    awk '$5 < 0.05 && $7 > 0 { print $1 "\t" $2 "\t" $3 }' $methylSigFile > $dmUp

    # DM down
    awk '$5 < 0.05 && $7 < 0 { print $1 "\t" $2 "\t" $3 }' $methylSigFile > $dmDown

    # DM all needed for no DM and no signal
    cat $dmUp $dmDown > $dmAll
    sort -T . -k1,1 -k2,2n $dmAll > $dmAllSorted
    mv $dmAllSorted $dmAll

    # No DM but signal
    awk '$5 > 0.05 { print $1 "\t" $2 "\t" $3 }' $methylSigFile > $noDMSignal
    sed -i '1d' $noDMSignal
    sort -T . -k1,1 -k2,2n $noDMSignal > $noDMSignalSorted
    mv $noDMSignalSorted $noDMSignal

    # No DM and no signal
    # The inverse regions of $methylSigFile
    awk '{ print $1 "\t" $2 "\t" $3 }' $methylSigFile > $methylSigSmall
    sed -i '1d' $methylSigSmall

    # We make the assumption that everything not in the methylSig output is
    # not DM and has no signal.
    bedtools complement -i $methylSigSmall -g ~/latte/Homo_sapiens/chromInfo_hg19.txt.gz > $noDMNoSignal
    sort -T . -k1,1 -k2,2n $noDMNoSignal > $noDMNoSignalSorted
    mv $noDMNoSignalSorted $noDMNoSignal

        rm $methylSigSmall
        rm $dmAll

    # These four files are mutually disjoint
        # $dmUp
        # $dmDown
        # $noDMSignal
        # $noDMNoSignal

        echo 'All of the following values should be 0!'
        echo 'Up regions disjoint from down regions?'
        bedtools intersect -a $dmUp -b $dmDown | wc -l
        echo 'Up regions distinct from no DM and signal regions?'
        bedtools intersect -a $dmUp -b $noDMSignal | wc -l
        echo 'Up regions distinct from no DM and no signal regions?'
        bedtools intersect -a $dmUp -b $noDMNoSignal | wc -l
        echo 'Down regions distinct from no DM and signal regions?'
        bedtools intersect -a $dmDown -b $noDMSignal | wc -l
        echo 'Down regions distinct from no DM and no signal regions?'
        bedtools intersect -a $dmDown -b $noDMNoSignal | wc -l
        echo 'No DM signal regions distinct from no DM no signal regions?'
        bedtools intersect -a $noDMSignal -b $noDMNoSignal | wc -l
