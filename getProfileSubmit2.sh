#CHROM=#${SGE_TASK_ID} #original command, not working on this server

#echo $CHROM

NUCLDIR=
OUTDIR=


for CHROM in $(seq 1 22)
	do
	    #export CHROM=$CHROM
	    #$SCRIPT_GETPROSUB # à enlever



$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex5.$CHROM $SIMPLE.$CHROM.start.fwd 0 $OUTDIR/$SIMPLER.$CHROM.fwd
$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex5.$CHROM $SIMPLE.$CHROM.start.rev 1 $OUTDIR/$SIMPLER.$CHROM.rev
#$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex3.$CHROM $SIMPLE.$CHROM.start.fwd 0 $OUTDIR/$SIMPLER.$CHROM.fwd
#$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex3.$CHROM $SIMPLE.$CHROM.start.rev 1 $OUTDIR/$SIMPLER.$CHROM.rev


$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex5.$CHROM $SIMPLE.$CHROM.start.fwd 1 $OUTDIR/$SIMPLER.$CHROM.ifwd
$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex5.$CHROM $SIMPLE.$CHROM.start.rev 0 $OUTDIR/$SIMPLER.$CHROM.irev
#$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex3.$CHROM $SIMPLE.$CHROM.start.fwd 1 $OUTDIR/$SIMPLER.$CHROM.ifwd
#$SCRIPT_PYTHON $SCRIPT_GETPRO $NUCLDIR/nucl_ex3.$CHROM $SIMPLE.$CHROM.start.rev 0 $OUTDIR/$SIMPLER.$CHROM.irev

done