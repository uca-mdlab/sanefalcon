#!/usr/bin/env bash
BAMFOLDER=/home/marco/temp/sanefalcon_bam  # will be /results/analysis/output/Home
BAMLINKFOLDER_BASE=/home/marco/temp/train

TRAINFOLDER=$(/usr/bin/python3 ./prepare_folders.py $BAMFOLDER $BAMLINKFOLDER_BASE)

if [[ -z $TRAINFOLDER ]]; then
    echo "No trainfolder found. Aborting"
    exit 1
fi

function prepSamples(){
    ./prepSamples.sh $BAMFOLDER $TRAINFOLDER
}


# all logic steps go there
prepSamples