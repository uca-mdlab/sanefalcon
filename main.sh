#!/usr/bin/env bash
BAMFOLDER=/home/marco/temp/sanefalcon_bam  # will be /results/analysis/output/Home
BAMLINKFOLDER_BASE=/home/marco/temp/train

TRAINFOLDER=$(/usr/bin/python3 ./main.py $BAMFOLDER $BAMLINKFOLDER_BASE)



function prepSamples(){
    ./prepSamples.sh $BAMFOLDER $TRAINFOLDER
}


# all logic steps go there
prepSamples