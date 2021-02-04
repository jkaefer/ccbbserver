#!/bin/bash
INPUT=$1
OUTPUT=$2
LOG=$3
SAMPLE=$INPUT
WD=${LOG/log}
SCRIPT="export PATH=/home/linhai/bin/src/anaconda2/bin:/home/linhai/bin/src/bedtools2/bin:$PATH; regsnp_intron -f -s /data/regsnp_intron_web/settings.json $INPUT $OUTPUT > $LOG 2>&1"
echo $SCRIPT > $SAMPLE.job
qsub -l nodes=1:ppn=1,walltime=24:00:00 -d $PWD -o $WD -e $WD -M linhai@iupui.edu -m abe -N $SAMPLE $SAMPLE.job
# rm $SAMPLE.job

