#$ -S /bin/bash
export SCRIPT_GETPRO=getProfile.py
export SCRIPT_GETPROSUB=/home/ionadmin/sanefalcon/getProfileSubmit.sh

export SCRIPT_PYTHON=python

export INDIR=`readlink -f $1`
echo $INDIR "+++++++++++++++++++"
export NUCLDIR=$2
echo $NUCLDIR "-----------------"
export OUTDIR=`readlink -f $3`

mkdir $OUTDIR

for SAMPLE in `find $INDIR -name "*.1.start.fwd"`
do
	export SIMPLE=${SAMPLE//.1.start.fwd/}
	export SIMPLER=`echo $SIMPLE | rev | cut -d"/" -f1 | rev`
	#echo $SIMPLE "THIS IS SIMPLE"
        #echo $SIMPLER "THIS IS SIMPLER"
	# qsub -t 1-22:1 -V $SCRIPT_GETPROSUB # original command, not working on this server...
	for CHROM in $(seq 1 22)
	do
	    export CHROM=$CHROM
	    $SCRIPT_GETPROSUB # à enlever

	done
done
