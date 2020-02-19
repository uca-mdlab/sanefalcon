#!/bin/bash


NOW=`date +'%d-%m-%Y_%Hh%M'`
BKPDIR="backup${NOW}"

mkdir -p $BKPDIR
echo "Packing and deleting from app folder"
tar czvf ${BKPDIR}/images.tar.gz testingNucleosome.csv* train_curvefit.png trainingNucleosome.csv*
tar czvf ${BKPDIR}/logs.tar.gz logs/*
grep 'ff(' logs/sanefalcon.log > ${BKPDIR}/results.ff


rm testingNucleosome.csv*
rm train_curvefit.png
rm trainingNucleosome.csv*
rm logs/*
echo "Done"
