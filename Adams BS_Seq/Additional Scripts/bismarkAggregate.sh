#!/bin/bash

export LC_ALL=C

command -v python >/dev/null 2>&1 || { echo >&2 "I require python but it's not installed.  Aborting."; exit 1; }


echo 'usage: bismarkPipeline <outputFolder> <files to aggregate>*'

export PYTHONPATH="/mnt/50tb/repository/shared/baselib/:."

export aggregateBismark="/mnt/50tb/repository/Tony/BisulfiteSeqBismark/aggregateBismark.py"

outputfolder=$1

echo "Output Folder : $outputfolder"
mkdir $outputfolder

shift

if [ -a $outputfolder/CpG_context_combined.txt ];
then
	echo "No merging and sorting"
else
	# for each of the input files
	for i in "$@"
	do
		echo $i
		zcat $i | awk 'FNR > 1' >> $outputfolder/CpG_context_combined.txt
	done

	# sort for UCSC plots
	sort -k3,3 -k4,4n -S15G $outputfolder/CpG_context_combined.txt > $outputfolder/CpG_context_combined.sorted.txt

	#rm $outputfolder/CpG_context_combined.txt
fi

python $aggregateBismark --bismark $outputfolder/CpG_context_combined.sorted.txt --aggregate $outputfolder/CpG_context_aggregate.txt

rm $outputfolder/CpG_context_combined.sorted.txt
rm $outputfolder/CpG_context_combined.txt
