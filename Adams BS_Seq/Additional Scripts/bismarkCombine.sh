#!/bin/bash

export LC_ALL=C

command -v python >/dev/null 2>&1 || { echo >&2 "I require python but it's not installed.  Aborting."; exit 1; }

echo 'usage: bismarkCombine <outputFile> <aggregate files to combine>*'

export PYTHONPATH="/mnt/50tb/repository/shared/baselib/:."

export combineReplicates="/mnt/50tb/repository/Tony/BisulfiteSeqBismark/combineReplicates.py"

outputfile=$1

shift

echo $1
cp $1 $outputfile

shift

# for each of the remaining input files
for i in "$@"
do
	echo $i
	python $combineReplicates --one $outputfile --two $i --combine "tmp"
	mv tmp $outputfile
done

