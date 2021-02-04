#!/bin/sh

#export WEKAHOME=/opt/weka/
#export WEKAHOME=~/regsnps-splicing/source/code_package/weka-3-7-12/
#export CLASSPATH=$WEKAHOME/weka.jar

echo 'runningweka'

java  -cp ~/regsnps-splicing/source/code_package/weka-3-7-12/weka.jar weka.classifiers.trees.RandomForest -i  -l $1 -T $2 -threshold-file $4/$3.tth2.arff -threshold-label 1 -p 0 > $4/$3.prediction

#echo "$1 MCC:";
#compute MCC
#perl ./source/code_package/model/mcc.pl $4/$3.tth2.arff $4/$3.prediction
awk 'NR > 17 {print $0}' $4/$3.tth2.arff > $4/stat.$3

