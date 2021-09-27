#!/bin/sh

#####################################################################################################
## This script reads all fastq files from a folder and performs quality control, genome alignment
## and read summarisation to quantify read counts 
## Script for analysis of RNA-seq data 
## Trim galore for trimming low quality reads fastq files
## Tophat/Bowtie for aligning reads to the updated reference genome (GRCh38) from ensemble 
### b - pre perfused kidney patients
### b1 - post perfused kidney patients
### DGF - Delayed graft function conditioned patients
### IGF - No Delayed graft function patients 
####################################################################################################
## Script modified : 25/03/2016
## Created: 2015 
## Copyright Suhaib Mohammed, 2015
## This script is not guaranteed to be free of bugs and/or errors.
## This script cannot be freely used and shared for any commercial purpose.

## @ARGV 
export rootdir=$1
export genomepath=$2
export gtfpath=$3
export fastqfilespath=$4
export numberofthreads=$5

echo  "\n..............Loading all modules...............\n 
..........................................................\n"

## linking path to software which are used in the analysis

    echo  "\nLoading cutadapt-1.7.1....\n...."
    PATH="/home/smm30e/Documents/Softwares/AdapterTrimming/cutadapt-1.7.1/bin:$PATH"
	. ~/.bashrc
	
    echo  "\nLoading Tophat-2.0.13.... \n...."
	PATH="/home/smm30e/Documents/Softwares/tophat-2.0.13.Linux_x86_64:$PATH"
	. ~/.bashrc

    echo  "\nLoading Bowtie2-2.2.3.... \n...."
	PATH="/home/smm30e/Documents/Softwares/bowtie2-2.2.3:$PATH"
	. ~/.bashrc
	
    echo  "\nLoading Samtools-1.1....\n...."
	PATH="/home/smm30e/Documents/Softwares/samtools-1.1:$PATH"
	. ~/.bashrc
 
	 export trimgalorepath="/home/smm30e/Documents/Softwares/AdapterTrimming/trim_galore_zip/trim_galore"
	 export featurecountspath="/home/smm30e/Documents/Softwares/subread-1.4.6-p5-Linux-x86_64/bin/featureCounts/"
 
## linking indexed human genome to the working directory 
	ln -s  $gtfpath/Homo_sapiens.GRCh38.78.gtf  Homo_sapiens.GRCh38.78.gtf
	ln -s  $genomepath/GRCh38.all2.1.bt2 GRCh38.all2.1.bt2
	ln -s  $genomepath/GRCh38.all2.2.bt2 GRCh38.all2.2.bt2
	ln -s  $genomepath/GRCh38.all2.3.bt2 GRCh38.all2.3.bt2
	ln -s  $genomepath/GRCh38.all2.4.bt2 GRCh38.all2.4.bt2
	ln -s  $genomepath/GRCh38.all2.rev.1.bt2 GRCh38.all2.rev.1.bt2
	ln -s  $genomepath/GRCh38.all2.rev.2.bt2 GRCh38.all2.rev.2.bt2
	ln -s  $genomepath/Homo_sapiens.GRCh38.dna.chromosome.all2.fa GRCh38.all2.fa


## Trim galore program used to trim paired end low quality reads and remove adapter sequences
## removing 3 bases from 3'end from all the reads as was biased to base T and A.
echo  "Creating new directory Trimmed \n"

### Trimming adapter sequences and removing low quality bases
if [ -d $1/trimmed ];
then
	echo "Trimmed folder already exists - Not trimming";
else
	mkdir $1/trimmed
	
      for i in $4/*_R1*.fastq;
       do 
       
        echo "i = $i"
		one=$i
		two=`echo $one | sed "s/_1/_2/"`
		none=`echo $one | sed "s/_1//" | sed "s/.fastq//" | sed "s/.gz//" | sed "s/.*\///"`
       
        echo "read 1 = $one"
        echo "read 2 = $two"
        
		echo $none

        outputfolder="$1/trimmed/$none"
		mkdir $outputfolder
		
       $trimgalorepath --fastqc --output_dir trimmed --phred33 -quality 20 --adapter 'AGATCGGAAGAGC' --three_prime_clip_R1 3 \
                         --three_prime_clip_R2 3 --paired $one $two                          
       done                  

	mkdir $1/trimmed/fastqc
	mkdir $1/trimmed/reports

	mv *.html $1/trimmed/fastqc/
	mv *.txt $1/trimmed/reports/
fi

########
## GTF file contains gene annotation corresponding to GRCh38
# do the alignment
if [ -d $1/alignment ];
then
	echo "Alignments folder already exists - Not aligning";
else
	mkdir $1/alignment

# Aligning to the reference genome
for i in $1/trimmed/*_1*.fq*
	do
	    echo "i = $i"
		one=$i
		two=`echo $one | sed "s/_1/_2/"`
		none=`echo $one | sed "s/_1//" | sed "s/.fq//" | sed "s/.gz//" | sed "s/.*\///"`
       
        echo "read 1 = $one"
        echo "read 2 = $two"
        
		echo $none

        outputfolder="$1/alignment/$none"
		mkdir $outputfolder
		
    echo  "\n Running Tophat...\n"
 

    tophat -r 158 --num-threads $numberofthreads --GTF Homo_sapiens.GRCh38.78.gtf GRCh38.all2 $one $two


    echo  "\n Tophat completed !!! \n"
       mv $1/tophat_out $outputfolder
     
	 done
fi

########
### Reads summarisation using featureCounts
if [ -d $1/counts ];
then
	echo "Counts folder already exists - Not counting";
else
	mkdir $1/counts
	
 echo "\n Running featureCount...!!!\n"
 
   for FILENAME in $1/alignment/**/*.bam;
       do 
           echo "Feature counts = $FILENAME"
           
               $featurecountspath -p -s 2 -M -O \
               -a Homo_sapiens.GRCh38.78.gtf  -o $1/counts/${FILENAME%_*}_counts.txt -t exon -g gene_id -T $numberofthreads $FILENAME
       done 
fi
########
