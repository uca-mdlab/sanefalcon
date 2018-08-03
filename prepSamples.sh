#!/usr/bin/env bash
INDIR=$1
OUTDIR=$2


#DIR_FETALFRAC=/home/davidpratella/fetalfrac
SCRIPT_RETRO=./retro.py

#DIR_SCRIPTS=$DIR_FETALFRAC
#export DIR_SCRIPTS=/illumina/diagnostics/bin
#export SCRIPT_PYTHON=/illumina/diagnostics/bin/Python-2.7.3/bin/python
SCRIPT_PYTHON=/usr/bin/python3

SCRIPT_SAMTOOLS=samtools

for SAMPLE in `find $INDIR -name "*.bam"`
do
#	SHORT=`echo $SAMPLE | rev | cut -d"/" -f1 | rev`
	SHORT=`echo $(basename $SAMPLE)`
	OUTFILE="$OUTDIR/${SHORT//".sort.bam"/}"
	echo "IN:$SAMPLE, OUT:$OUTFILE"
	
	for ARG_TASKID in `seq 1 2` # or "X"
		do
		    if [ ! -f $OUTFILE.$ARG_TASKID.start.fwd ]; then
			    $SCRIPT_SAMTOOLS view \
			    $SAMPLE \
			    chr$ARG_TASKID \
			    -F 20 -q 1 | \
			    $SCRIPT_PYTHON $SCRIPT_RETRO | \
				    awk '{print $4}' \
					        > $OUTFILE.$ARG_TASKID.start.fwd &

                $SCRIPT_SAMTOOLS view \
                $SAMPLE \
                chr$ARG_TASKID \
                -f 16 -F 4 -q 1 | \
                $SCRIPT_PYTHON $SCRIPT_RETRO | \
                        awk '{print ($4 + length($10) - 1)}' \
                            > $OUTFILE.$ARG_TASKID.start.rev
            fi
		done
done
