#!/usr/bin/env bash
BAMFOLDER=/home/marco/temp/sanefalcon_bam  # will be /results/analysis/output/Home
BAMLINKFOLDER_BASE=/home/marco/temp/train

INDIR=$1
OUTDIR=$2


function prepSamples(){
    ./prepSamples.sh $INDIR $OUTDIR
}


# all logic steps go there
prepSamples