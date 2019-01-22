export SCRIPT_NUCDEC=nuclDetector.py

export SCRIPT_PYTHON=python

export INDIR=`readlink -f $1`

for CHROM in $(seq 22 -1 1)
do
#	echo "$SCRIPT_PYTHON $SCRIPT_NUCDEC $INDIR/merge.$CHROM $INDIR/nucl_ex4.$CHROM" | qsub
	$SCRIPT_PYTHON $SCRIPT_NUCDEC $INDIR/merge.$CHROM $INDIR/nucl_ex4.$CHROM
done

