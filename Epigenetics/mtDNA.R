
## Script modified : 16/10/2016
## Created: 2016 
## Copyright Suhaib Mohammed, 2016

library(Biostrings)

mt.dna<-readDNAStringSet("/media/smm30e/Seagate Backup Plus Drive/WGBS_seq/Genome2/Homo_sapiens.GRCh38.dna.chromosome.MT.fa")
seq<-paste(mt.dna)
m1 <- as.data.frame(unlist(vmatchPattern("CG", mt.dna,fixed=TRUE)))
seq = paste(seq, collapse="")
write.table(m1,file="/media/smm30e/Backup Data/WGBS/Analysis/MTmethylation/MtDNA_CpGs.txt",sep="\t",col.names=NA)

Views(seq,start=c(1213), end=c(1218))
Views(seq, start=c(1313), end=c(1315))

# "CpG start at 1216, 1 based cordinate ~ 1215 zero based"
subseq(mt.dna, start=c(1213), end=c(1218))

# "CpG start at 1314, 1 based cordinate ~ 1313 zero based"
subseq(mt.dna, start=c(1313), end=c(1315))

## frequency of dinucleotide
dinucleotideFrequency(mt.dna)

fasta2dataframe=function(fastaFile){
  s = readDNAStringSet(fastaFile)
  RefSeqID = names(s)
  RefSeqID = sub(" .*", "", RefSeqID) 
  #erase all characters after the first space: regular expression matches a space followed by any sequence of characters and sub replaces that with a string having zero  characters 
  
  for (i in 1:length(s)){
    seq[i]=toString(s[i])
  }
  
  RefSeqID_seq=data.frame(RefSeqID,seq)
  return(RefSeqID_seq)
}

myFastaFile.fasta<-"/media/smm30e/Seagate Backup Plus Drive/WGBS_seq/Genome2/Homo_sapiens.GRCh38.dna.chromosome.MT.fa"

mydf = fasta2dataframe(myFastaFile.fasta)