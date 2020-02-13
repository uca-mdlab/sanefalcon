#!/bin/bash


NOW=`date +'%d-%m-%Y_%Hh%M'`
BKPDIR="backup${NOW}"

mkdir -p $BKPDIR
echo "Packing from app folder..."
tar czf ${BKPDIR}/images.tar.gz testingNucleosome.csv* train_curvefit.png trainingNucleosome.csv*
tar czf ${BKPDIR}/logs.tar.gz logs/*
grep 'ff(' logs/sanefalcon.log > ${BKPDIR}/results.ff

echo "Deleting stuff..."
rm testingNucleosome.csv*
rm train_curvefit.png
rm trainingNucleosome.csv*
rm logs/*


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
    fi
    echo "popd"
    break

done
echo "Done. Training tree in ${SOURCE} restored in ${WORKDIR}"