

#for SAMPLE in ./klgen*_2/*/*.1.peaks.fwd2
#for SERIE in a b c d e f g h # i j k l;
#do
	#for SAMPLE in $1/test/*.1.rev.np
	#for SAMPLE in $1/$SERIE/*.1.rev.np
#	for SAMPLE in `find $1 -name "*.1.rev.np"`
	for SAMPLE in `find $1 -name "*.1.rev"`
#	for SAMPLE in `find $1 -name "*.1.rev"`
	do
#		SIMPLE=`echo ${SAMPLE//.1.rev.np/} | rev | cut -d"/" -f1 | rev`
		SIMPLE=`echo ${SAMPLE//.1.rev/} | rev | cut -d"/" -f1 | rev`
#		SIMPLE=`echo ${SAMPLE//rev.1.rev/} | rev | cut -d"/" -f1 | rev`

#		echo $SIMPLE,`LC_NUMERIC=C awk -F"," '{ for(i=1;i<=NF;i++){ fetal[i]+=$i } } END { for(i=1;i<=NF;i++) {printf "%s,",fetal[i]}; print "";}' \
#		${SAMPLE//.1.rev.np/}.*.irev.np ${SAMPLE//.1.rev.np/}.*.ifwd.np`
		echo $SIMPLE,`LC_NUMERIC=C awk -F"," '{ for(i=1;i<=NF;i++){ fetal[i]+=$i } } END { for(i=1;i<=NF;i++) {printf "%s,",fetal[i]}; print "";}' \
		${SAMPLE//.1.rev/}.*.irev ${SAMPLE//.1.rev/}.*.ifwd`
#		echo $SIMPLE,`LC_NUMERIC=C awk -F"," '{ for(i=1;i<=NF;i++){ fetal[i]+=$i } } END { for(i=1;i<=NF;i++) {printf "%s,",fetal[i]}; print "";}' \
#		${SAMPLE//.rev.1.rev/}.*.irev ${SAMPLE//.rev.1.rev/}.*.ifwd`
	done
#done
