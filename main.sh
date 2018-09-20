#!/usr/bin/env bash

set -x

BAMFOLDER=/results/analysis/output/Home
#BAMLINKFOLDER_BASE=/tmp/sanefalcontrain
BAMLINKFOLDER_BASE=/results/plugins/sanefalcon/sanefalcontrain
OUTPUT_FOLDER=/home/ionadmin/tmp_david/sanefalcon/getprofile_nucleosome1

# preparation step
#TRAINFOLDER=$(/usr/bin/python3 ./prepare_folders.py $BAMFOLDER $BAMLINKFOLDER_BASE)

if [[ -z $TRAINFOLDER ]]; then
    echo "No trainfolder found. Aborting"
    #exit 1
fi

# function definitions: 1 function per step
function prepSamples(){
    for subdir in `find $BAMLINKFOLDER_BASE -maxdepth 2 -mindepth 2 -type d`; do
        ./prepSamples.sh $subdir $subdir & #add '&' to continue
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

function nuclDetectorAnti(){
    for subdir in `find $BAMLINKFOLDER_BASE -maxdepth 1 -mindepth 1 -type d`; do
        if [ -d ${subdir} ]; then
            ./nuclDetectorAnti.sh $subdir   # a b c d ...
        fi
    done
}

<<<<<<< HEAD
function getProfile(){
    for subdir in `find $BAMLINKFOLDER_BASE -maxdepth 1 -mindepth 1 -type d`; do
        if [ -d ${subdir} ]; then
            ./getProfile.sh $BAMLINKFOLDER_BASE $subdir /home/ionadmin/tmp_david/sanefalcon/getprofile_nucleosome1  # a b c d ...
        fi
    done
}
=======
#function getProfile(){
#    for subdir in `find $BAMLINKFOLDER_BASE -maxdepth 1 -mindepth 1 -type d`; do
#        if [ -d ${subdir} ]; then
#            ./nuclDetectorAnti.sh $subdir   # a b c d ...
#        fi
#    done
#}

>>>>>>> d728833fdcec39de42c1ae6e74b2778ed4005065




# all logic steps go there
<<<<<<< HEAD
#prepSamples && echo "passed prepSamples"
#wait
#mergeSamples && echo "passed mergeSamples"
#wait
#mergeSubs && echo "passed mergeSubs"
#wait
#mergeAntiSubs && echo "passed mergeAntiSubs"
#wait
#nuclDetectorAnti && echo "passed mergeAntiSubs"
getProfile && echo "passed getProfile"
=======
prepSamples && echo "passed prepSamples"
wait
mergeSamples && echo "passed mergeSamples"
wait
mergeSubs && echo "passed mergeSubs"
wait
mergeAntiSubs && echo "passed mergeAntiSubs"
wait
nuclDetectorAnti && echo "passed nuclDetectorAnti"
wait
#./getProfile.sh $BAMLINKFOLDER_BASE $subdir /home/ionadmin/tmp_david/sanefalcon/getprofile_nucleosome1
#./getProfile.sh /results/plugins/sanefalcon/sanefalcontrain /results/plugins/sanefalcon/sanefalcontrain/a /home/ionadmin/tmp_david/sanefalcon/getprofile_nucleosome2
/results/plugins/sanefalcon/sanefalcontrain/a
>>>>>>> d728833fdcec39de42c1ae6e74b2778ed4005065
echo "Done"
