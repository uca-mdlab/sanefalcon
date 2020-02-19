#!/bin/bash

WORKDIR=/home/mdlab/storage/sanefalcon/training
SOURCE=$1

echo "Rebuilding training tree from ${SOURCE}"

for d in $SOURCE/training/*/ ; do
    name=`basename $d`
    mkdir $WORKDIR/$name
    pushd $WORKDIR/$name
    if [ "$name" == "profiles" ]; then
        for f in `ls $SOURCE/training/$name/p*`; do
            ln $f
        done
        for f in `ls $SOURCE/training/$name/s*`; do
            ln $f
        done
    else
        for f in `ls $SOURCE/training/$name/*.*`; do
            ln $f
        done
    fi
    popd

done
echo "Done. Training tree in ${SOURCE} restored in ${WORKDIR}"
