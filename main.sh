#!/usr/bin/env bash
BAMFOLDER=/home/marco/temp/sanefalcon_bam  # will be /results/analysis/output/Home
BAMLINKFOLDER_BASE=/home/marco/temp/train

letters=( {a..z} )
rnd=$(( ($RANDOM % 26 ) + 1 ))

# find current subdirs in the train folder
subdir_list=()
while IFS= read -d $'\0' -r subdir ; do
    subdir_list=("${subdir_list[@]}" "$subdir")
done < <(find $BAMLINKFOLDER_BASE -mindepth 1 -maxdepth 1 -type d -print0; )

if [ ${#subdir_list[@]} -eq 0 ]; then
    last_dir=''
else
    last_dir="${subdir_list[-1]##*/}"
fi

function get_index(){
    haystack=$1
    needle=$2
    if [[ ${needle} = "" ]]; then
        return 0
    fi
    for i in "${!haystack[@]}"; do
        if [[ "${haystack[$i]}" = "${needle}" ]]; then
            break
        fi
    done

    new_index=$(( $i + 1 ))
    return $new_index

}


get_index $letters $last_dir
new_index=$?

new_folder=${letters[$new_index]}
mkdir $BAMLINKFOLDER_BASE/$new_folder

exit 0

for letter in {a..z} ; do
    if [ "$letter" == "$last_dir" ]; then
        echo 'ok'
    fi
done

exit 0

INDIR=$1
OUTDIR=$2


function prepSamples(){
    ./prepSamples.sh $INDIR $OUTDIR
}


# all logic steps go there
prepSamples