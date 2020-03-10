#!/bin/bash

WORKDIR=/home/mdlab/storage/sanefalcon/training
SOURCE=$1

echo "Rebuilding training tree from ${SOURCE}"

for d in "$SOURCE"/training/*/ ; do
    name=$(basename "$d")
    mkdir $WORKDIR/"$name"
    pushd $WORKDIR/"$name" || { echo "Failure"; exit 1; }
    if [ "$name" == "profiles" ]; then
        for f in "$SOURCE"/training/"$name"/p*; do
            [[ -e "$f" ]] || break
            ln "$f" .
        done
        for f in "$SOURCE"/training/"$name"/s*; do
            [[ -e "$f" ]] || break
            ln "$f" .
        done
    else
        for f in "$SOURCE"/training/"$name"/*.*; do
            [[ -e "$f" ]] || break
            ln "$f" .
        done
    fi
    popd || { echo "Failure"; exit 1; }

done
echo "Done. Training tree in ${SOURCE} restored in ${WORKDIR}"
