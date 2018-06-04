INDIR=`readlink -f $1`

# for testing locally only
ISQSUB=`command -v qsub`


if [[ $ISQSUB -eq 0 ]]; then
    for i in $(seq 1 3)  # for testing locally only (change to 1 22)
    do
	    echo Working on Chr: $i
        sort -n -m $INDIR/*/merge.$i > $INDIR/merge.$i
    done
else
    for i in $(seq 1 22)
    do
	    echo Working on Chr: $i
	    echo "sort -n -m $INDIR/*/merge.$i > $INDIR/merge.$i" | qsub
	done
fi
