#!/bin/bash

WORKDIR=/home/mdlab/storage/sanefalcon/training
SOURCE=$1

for d in $SOURCE/training/*/ ; do
    name=`basename $d`
    mkdir $WORKDIR/$name
done
