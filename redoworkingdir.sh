#!/bin/bash

#if [ "$#" -le 1 ]; then
#    echo "Illegal number of parameters"
#fi


NOW=`date +'%d-%m-%Y_%Hh%M'`
BKPDIR="backup${NOW}"

echo "mkdir -p $BKPDIR"
echo "Packing and deleting from app folder"
echo "tar czvf ${BKPDIR}/images.tar.gz testingNucleosome.csv* train_curvefit.png trainingNucleosome.csv*"
echo "tar czvf ${BKPDIR}/logs.tar.gz logs/*"
echo "grep 'ff(' logs/sanefalcon.log > ${BKPDIR}/results.ff"

echo "rm testingNucleosome.csv*"
echo "rm train_curvefit.png"
echo "rm trainingNucleosome.csv*"
echo "rm logs/*"


WORKDIR=/home/mdlab/storage/sanefalcon/training
SOURCE=$1

echo "Rebuilding training tree from ${SOURCE}"

for d in $SOURCE/training/*/ ; do
    name=`basename $d`
    echo "mkdir $WORKDIR/$name"
    echo "pushd $WORKDIR/$name"
    if [ "$name" == "profiles" ]; then
        for f in `ls $SOURCE/training/$name/p*`; do
            echo "ln $f"
        done
        for f in `ls $SOURCE/training/$name/s*`; do
            echo "ln $f"
        done
    else
        for f in `ls $SOURCE/training/$name/*.*`; do
            echo "ln $f"
        done
        #for f in `ls $SOURCE/training/$name/*anti*`; do
        #    echo "ln $f"
        #done
        #for f in `ls $SOURCE/training/$name/merge.*`; do
        #    echo "ln $f"
        #done
        #for f in `find $SOURCE/training/$name/ -type l`; do
        #    echo "cp $f ."
        #done
    fi
    echo "popd"
    break

done
echo "Done. Training tree in ${SOURCE} restored in ${WORKDIR}"