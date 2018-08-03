#!/usr/bin/env bash

set -x

BAMFOLDER=/results/analysis/output/Home
BAMLINKFOLDER_BASE=/tmp/sanefalcontrain

# preparation step
TRAINFOLDER=$(/usr/bin/python3 ./prepare_folders.py $BAMFOLDER $BAMLINKFOLDER_BASE)

if [[ -z $TRAINFOLDER ]]; then
    echo "No trainfolder found. Aborting"
    exit 1
fi

# function definitions: 1 function per step
function prepSamples(){
    for subdir in `find $BAMLINKFOLDER_BASE -maxdepth 2 -mindepth 2 -type d`; do
        ./prepSamples.sh $subdir $subdir
    done
}

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
prepSamples && echo "passed prepSamples"
#mergeSamples && echo "passed mergeSamples"
#mergeSubs && echo "passed mergeSubs"
#mergeAntiSubs && echo "passed mergeAntiSubs"
echo "Done"