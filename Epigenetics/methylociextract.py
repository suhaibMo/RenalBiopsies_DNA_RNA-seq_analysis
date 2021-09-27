## Script modified : 16/10/2016
## Created: 2016 
## Copyright Suhaib Mohammed, 2016

from sys import argv
import getopt
import sys
import csv
import pandas as pd
import numpy as np

script, from_file, loci_file, to_file = argv

print "Script to calculate methylation status from %s and %s \n " % (from_file, loci_file)

print "The input script is called:", script
print "\n Your CpG context file is:", from_file
print " \n Methylation loci file:", loci_file
print " \n Extracted output file", to_file


# Reading Input file as dataframe
#indata = open(from_file, 'r').read()
# from_file="~/Documents/Glasgow/methyl_loci/CpG_context_74B1_NoDGF.fastq.gz_bismark_bt2.all.cov"
indata=pd.DataFrame(pd.read_csv(from_file, sep='\t', names=['chr','start','end','methcov','meth','unmeth'],header=None))

#loci_file="~/Documents/Glasgow/methyl_loci/methy_clock_loci.txt"
locidata=pd.DataFrame(pd.read_csv(loci_file, sep='\t', header=0))
locidata.rename(columns={'seqnames': 'chr'},inplace=True)

### total summary 
totalreads=indata.shape[0]
methtot=indata[:]["meth"].sum()
unmethtot=indata[:]["unmeth"].sum()
ratio_methunmeth=round((float(methtot)/float(unmethtot)),3)
ratio_index=round((float(methtot)/(float(unmethtot)+float(methtot))),3)

print "total read coverage = %d " % totalreads
print "total methylated coverage = %d" % methtot
print "total unmethylated coverage = %d" % unmethtot
print "ratio meth-unmeth = %f" % ratio_methunmeth
print "ratio index = %f" % ratio_index

##########
# extracting loci of methyl clock gene cordinates
locidata2=locidata[["start","end","strand","external_gene_name"]]


# defining empty dataframe 
dfmeth = pd.DataFrame(index=range(0,len(locidata2)), columns=('gene',"meth","unmeth","ratio","ratioindex"))
dfmeth = dfmeth.astype(float).fillna(0)
dfmeth.head()

for i in range(0,len(locidata2)):
    print "start loci = %d" % locidata2.start[i]
    print "end loci = %d " % locidata2.end[i]
    print "strand = %s" % locidata2.strand[i]
    print "gene name = %s" % locidata2.external_gene_name[i]
    
    if locidata2.strand[i]=='+' or locidata2.strand[i]=='-':
        methext=indata[(indata.start>=locidata2.start[i]) & (indata.start<=locidata2.end[i])]
        
    dfmeth.gene[i]=locidata2.external_gene_name[i]
    dfmeth.meth[i]=methext["meth"].sum()
    dfmeth.unmeth[i]=methext["unmeth"].sum()
    dfmeth.ratio[i]=round((float(dfmeth.meth[i])/float(dfmeth.unmeth[i])),3)
    dfmeth.ratioindex[i]=round((float(dfmeth.meth[i])/(float(dfmeth.meth[i])+float(dfmeth.unmeth[i]))),3)
      
##########       
#to_file=from_file.split(".")[0]
outpath=r'/media/smm30e/Backup\ Data/WGBS/Analysis/HovarthMethloci'

out_file = open(to_file, 'w')
print "writing output to %s" % out_file
out_file.write(outpath+str(dfmeth))

print "Alright, all done."
out_file.close()    
