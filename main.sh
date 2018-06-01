#!/usr/bin/env bash

#set -x

BAMFOLDER=/home/marco/temp/sanefalcon_bam  # will be /results/analysis/output/Home
BAMLINKFOLDER_BASE=/home/marco/temp/train

# preparation step
#TRAINFOLDER=$(/usr/bin/python3 ./prepare_folders.py $BAMFOLDER $BAMLINKFOLDER_BASE)

#if [[ -z $TRAINFOLDER ]]; then
#    echo "No trainfolder found. Aborting"
#    exit 1
#fi


# function definitions: 1 function per step
#function prepSamples(){
#    ./prepSamples.sh $BAMFOLDER $TRAINFOLDER
#}

function mergeSamples(){
    for subdir in `find $BAMLINKFOLDER_BASE -maxdepth 1 -mindepth 1 -type d`; do
        if [ -d ${subdir} ]; then
            ./merge.sh $subdir
        fi
    done
}

function mergeSubs(){
    ./mergeSubs.sh $BAMLINKFOLDER_BASE
}

function mergeAntiSubs(){
    for subdir in `find $BAMLINKFOLDER_BASE -maxdepth 1 -mindepth 1 -type d`; do
        if [ -d ${subdir} ]; then
            ./mergeAntiSub.sh $BAMLINKFOLDER_BASE ${subdir##*/}  # a b c d ...
        fi
    done
}

# all logic steps go there
#prepSamples && echo "passed prepSamples"
mergeSamples && echo "passed mergeSamples"
mergeSubs && echo "passed mergeSubs"
mergeAntiSubs && echo "passed mergeAntiSubs"
echo "Done"