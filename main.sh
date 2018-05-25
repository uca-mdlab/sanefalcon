#!/usr/bin/env bash
BAMFOLDER=/home/marco/temp/sanefalcon_bam  # will be /results/analysis/output/Home
BAMLINKFOLDER_BASE=/home/marco/temp/train

/usr/bin/python3 ./main.py $BAMFOLDER $BAMLINKFOLDER_BASE


exit 0

function prepSamples(){
    ./prepSamples.sh $INDIR $OUTDIR
}


# all logic steps go there
prepSamples