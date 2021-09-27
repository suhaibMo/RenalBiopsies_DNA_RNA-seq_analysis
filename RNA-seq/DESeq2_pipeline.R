#####################################################################################################
#### This script reads all counts.txt files from a folder and extracts gene ids and respective counts 
#### from each sample and performes differential gene expression analysis against disease and control
### b - pre perfused kidney patients
### b1 - post perfused kidney patients
### DGF - Delayed graft function conditioned patients
### IGF - No Delayed graft function patients 
####################################################################################################
##  Script modified : 04/04/2016
##  Author : Suhaib Mohammed
## Copyright Suhaib Mohammed, 2015
## This script is not guaranteed to be free of bugs and/or errors.
## This script cannot be freely used and shared for any commercial purpose.
#####################################################################################################

require(foreach)
path<-"/media/smm30e/New Volume1/RNA_seq/Counts_overlap2/"
all_files <-list.files(path,pattern="\\counts.txt$", full.names=T)
dset<-foreach(i=1:length(all_files), .combine =cbind) %do% read.table(all_files[i], header=TRUE,sep="\t")
dset<-dset[, !names(dset) %in% c("Chr","Start","End","Length","Strand","X")]
dset<-dset[,-grep("Geneid.", colnames(dset))]

### DGF- Delayed graft function patients samples id
DGF.samples<-c("97b","106b","128b","136b","138b","140b","150b","156b","169b","181b","36b","70b",
               "97b1","106b1","128b1","136b1","138b1","140b1","150b1","156b1","169b1","181b1","36b1","70b1")

dset.counts<-dset
head(dset.counts)

## trimming extra annotions from the colums for simplicity
colnames(dset.counts)<-gsub(pattern = ".sm.|.bam|bam", replacement = "",
                            x = colnames(dset.counts))
## trimming extra annotions from the colums for simplicity
colnames(dset.counts)<-gsub(pattern = "\\..", replacement = "",
                            x = colnames(dset.counts))
colnames(dset.counts)<-tolower(colnames(dset.counts))

########################################
## appending colum names with condition type: DGF or IGF 
for ( i in 2:length(colnames(dset.counts)))
{
  
  # if the merged dataset doesn't exist, create it
  if (colnames(dset.counts)[i] %in%  DGF.samples) 
    colnames(dset.counts)[i] <- paste("DGF",colnames(dset.counts)[i],sep="_")
  else 
    colnames(dset.counts)[i] <- paste("IGF",colnames(dset.counts)[i],sep="_")
}

#### maing first colum as rownames for simplicity
dset.counts<-data.frame(dset.counts[,-1], row.names=dset.counts[,1])
head(dset.counts)
dim(dset.counts)

## separting all samples based on conditions
# b- pre-perfused kidney transplant patients
# b1 - post-perfused kidney transplant patients

b.all<- dset.counts[,grep("b$",colnames(dset.counts))]
b1.all<-dset.counts[,grep("1$",colnames(dset.counts))]

## ordering the table of matrix to keep b and b1 separately
dset.order<-dset.counts[,c(colnames(b.all),colnames(b1.all))]
head(dset.order)
dim(dset.order)

###########################
# Stratifying samples
DGF.len<-length(grep("^DGF",colnames(dset.order)))
IGF.len<-length(grep("^IGF",colnames(dset.order)))

b.len<-length(grep("b$",colnames(dset.order)))
b1.len<-length(grep("b1$",colnames(dset.order)))

## set working directory
setwd("/home/smm30e/Documents/Scripts/smallRNA/MiRNA/SubRead/DGFvsNoDGF")

###########################
## function to remove any spaces
nospaces<-function(sig)
{
  sig.nospaces<-matrix(nrow=nrow(sig))
  
  for (i in 1:nrow(sig))
  {
    str.len<-length(strsplit(noquote(gsub("\\s","",rownames(sig))), "")[[i]])
    sig.nospaces[i]<-paste(strsplit(noquote(gsub("\\s","",rownames(sig))), "")[[i]][2:str.len],sep="",collapse="")
  }
  
  return(sig.nospaces)
  
}
rownames(dset.order)<-nospaces(dset.order)
write.table(dset.order, sep="\t", col.names=NA, file ="Read_counts.tsv")

##################################################################
### DESeq2 method for differential expression ####################

require(DESeq2)

## creating phenotypic meta data for group index for kidney condition and type of pateints 
#  condition indicates state of kidney transplant (B-perperfused) or (B1-post perfused)  
#  type indicates No delayed graft function (No DGF) or delayed graft function (DGF)
cond.idx<-  rep(c("Pre", "Post"), c(b.len,b1.len))
type.idx<-sapply(strsplit(as.character(colnames(dset.order)), split="_"), "[", 1)
sample.idx<-sapply(strsplit(as.character(colnames(dset.order)), split="_"), "[", 2)


## Donor gender identification
male.samples<-c("34b","36b","39b","70b","128b","129b","138b","150b","151b","153b","156b","181b",
                "34b1","36b1","39b1","70b1","128b1","129b1","138b1","150b1","151b1","153b1","156b1","181b1") 

gender<-matrix()
for (i in 1:length(sample.idx))
{
  if (sample.idx[i] %in%  male.samples) 
    gender[i] <-"Male"
  else 
    gender[i] <-"Female"
} 

####  GSK patients database 
path.gsk <- "~/Glasgow/GSK/DifferentialGenes/gsk database.csv"
data.gsk <- read.csv(path.gsk, header = TRUE,sep=",")
samples<-gsub('b.*', '', sample.idx)

patients.db<-(data.gsk[match(samples,data.gsk$labStudy.ID),])

head(patients.db)   
##
# recepient gender
recep.gender<-patients.db$Recip.sex.0.m.1.f
recep.gender<-replace(recep.gender,which(patients.db$Recip.sex.0.m.1.f=="0"),"Male")
recep.gender<-replace(recep.gender,which(recep.gender=="1"),"Female")

coldat=data.frame(colnames=factor(colnames(dset.order)),condition=factor(cond.idx),
                  type=factor(type.idx),sample=factor(sample.idx), D.gender=factor(gender),R.gender=factor(recep.gender), row.names=1)
coldat

#################################################################
# Constructing data object from matrix of counts
ddsMat <- DESeqDataSetFromMatrix(dset.order, colData=coldat, design =~type)

dds <- estimateSizeFactors(ddsMat)

# Transforming data
rld <- rlogTransformation(ddsMat)
head(assay(rld))


# Plot
png(filename="rlogTransformation.png",width=1000,height=1000,res=150)
par( mfrow = c( 1,2 ) )
dds <- estimateSizeFactors(ddsMat)
plot( log2( 1 + counts(dds, normalized=TRUE)[ , c(1,25)]),
      col=rgb(0,0,0,.2), pch=16, cex=0.3,main=" log2Transformed" )
plot( assay(rld)[ , c(1,25)],
      col=rgb(0,0,0,.2), pch=16, cex=0.3, main=" rlogTransformed" )
dev.off()

##  dist calculates distances between data rows
##  and our samples constitute the columns.
sampledist <- dist(t(assay(rld)))

## simple Hierarchical clustering
tiff(filename="Sample_Clust.tiff",width=1000,height=1000,res=100,compression="lzw")
par( mfrow = c( 1,1 ) )
plot(hclust(sampledist),xlab="samples ID", main= "Sample cluster dendogram")
dev.off()

#####################################################
## Quality control of samples
#####################################################
# heatmap for sample dsitance
library("RColorBrewer")
library("gplots")

sampledistmat <- as.matrix(sampledist)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(555)
hc <- hclust(sampledist)

###
coldat$type.col<-  replace(coldat$type.col,which(coldat$type=="DGF"),"red")
coldat$type.col<-  replace(coldat$type.col,which(coldat$type=="IGF"),"blue")

coldat$D.gender.col<-replace(coldat$D.gender.col,which(coldat$D.gender=="Female"),"indianred1")
coldat$D.gender.col<-replace(coldat$D.gender.col,which(coldat$D.gender=="Male"),"deepskyblue")

coldat$cond.col<-replace(coldat$cond.col,which(coldat$condition=="Pre"),"seagreen")
coldat$cond.col<-replace(coldat$cond.col,which(coldat$condition=="Post"),"yellow")


hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

###
require(pheatmap)

df.cols <- as.data.frame(colData(rld)[,c("condition","D.gender")])
df.rows <- as.data.frame(colData(rld)[,c("R.gender","type")])

# Specify colors
ann_colors = list(
  type = c(DGF="red", IGF="blue"),
  condition = c(Pre="orange", Post="darkgreen"),
  D.gender = c(Male="deepskyblue", Female="indianred1"),
  R.gender = c(Male="darkolivegreen", Female="darkmagenta")
)

hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(100)

tiff(filename="pheatmap.sampleheatmap.samples_colorcodes2.tiff",width=1200,height=1000,res=120)
pheatmap(sampledistmat,
         clustering_distance_rows=sampledist,
         clustering_distance_cols=sampledist,
         col=rev(hmcol), annotation_col = df.cols, annotation_row = df.rows,
         annotation_colors = ann_colors)
dev.off()

tiff(filename="heatmap.samples.tiff",width=1000,height=1000,res=110,compression="lzw")
heatmap.2(sampledistmat, col = rev(hmcol), key = TRUE,
          #  Rowv=as.dendrogram(hc),symm=TRUE,
          trace = "none",margins=c(5.5,10))
dev.off()

### PCA

data <- plotPCA(rld, intgroup=c("condition", "type", "D.gender","R.gender"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
## pca plot gender

tiff(filename="PCA_Samples_D.gender2.tiff",width=700,height=700,res=100)
sp.D.gen<-ggplot(data, aes(PC1, PC2, color=D.gender)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

sp.D.gen + scale_color_manual(values=c("indianred1","deepskyblue"))
dev.off()


tiff(filename="PCA_Samples_R.gender2.tiff",width=700,height=700,res=100)
sp.R.gen<-ggplot(data, aes(PC1, PC2, color=R.gender)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

sp.R.gen + scale_color_manual(values=c("darkmagenta","darkolivegreen"))
dev.off()

## pca plot type
tiff(filename="PCA_Samples_type2.tiff",width=700,height=700,res=100)
sp.type<-ggplot(data, aes(PC1, PC2, color=type, shape=type)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

sp.type + scale_color_manual(values=c("red","blue"))
dev.off()

## pca plot condition
tiff(filename="PCA_Samples_cond2.tiff",width=700,height=700,res=100)
sp.type<-ggplot(data, aes(PC1, PC2, color=condition)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

sp.type + scale_color_manual(values=c("darkgreen","orange"))
dev.off()


# Principal components analysis using DESeq2
tiff(filename="PCA_Samples_condition_D.gender.tiff",width=700,height=700,res=100)
plotPCA(rld, intgroup = c("condition","D.gender"))
dev.off()

# Principal components analysis using DESeq2
tiff(filename="PCA_Samples_type_D.gender.tiff",width=700,height=700,res=100)
plotPCA(rld, intgroup = c("type","D.gender"))
dev.off()


#####################################################
## Diffential expression alaysis
#####################################################
### DGF vs IGF

# making sure that IGF is the first level in the type factor, 
# so that the default log2 fold changes are calculated as DGF over IGF
dds$type<-relevel(dds$type, "IGF")

dds<- DESeq(dds)
## comparing differentuial genes for DGF over No DGF
res<- results(dds,contrast=c("type","DGF","IGF"))

mcols(res, use.names=TRUE)
summary(res)

res.padj <- res[order(res$padj), ]

## FDR set to 10%
sig <- subset(res.padj, padj < 0.1)

sig

## MA plot
tiff(filename="MA_plot_DGFIGF.tiff",width=1000,height=1000,res=180)
DESeq2::plotMA(res, main="DGF vs IGF, FDR=10%", alpha=0.1, ylim=c(-2,2))
dev.off()

#### annotating gene signatures from differentially expressed genes
gene_ids<-read.table(file ="~/Glasgow/GSK/DifferentialGenes/Collated_counts/EnsemblIDs_GeneNames.tsv", sep="\t", header=TRUE)
gene.sig<-(gene_ids[(match(rownames(sig),gene_ids$ensembl_gene_id)),])


###########################
gene.sig_DGFIGF_v3<-gene.sig
save(gene.sig_DGFIGF_v3, file="sig_DGFvsIGFv3.rda")

tbl.gene.sig_DGF_IGF_v3<-cbind(sig,gene.sig[,c(2:ncol(gene.sig))])
write.table(tbl.gene.sig_DGF_IGF_v3, sep="\t", col.names=NA, file ="Diff_Exp_DGF_IGF_v3.xls")

rank.gene.DGF_v3<-cbind(as.data.frame(sig),gene.sig$external_gene_name,gene.sig$transcript_biotype)

png(filename="RankedSig.DGF_v3.png",width=2000,height=1000,res=150)
par(mar=c(10,5,4,2))
p<-barplot(rank.gene.DGF_v3$log2FoldChange,col="brown3",ylab="log fold change",xaxt="n", 
           main= "Ranked signatures",cex.axis=1.5)
axis(1, at=p, labels=gene.sig$external_gene_name,las=2,cex=0.5)
text(p,rank.gene.DGF_v3$log2FoldChange, labels = seq(1,length(rank.gene.DGF_v3$log2FoldChange),1), pos=1, cex=1)
dev.off()

save(rank.gene.DGF_v3,file="rank.gene.DGF_v3.rda")
write.table(rank.gene.DGF_v3,file="rank.gene.DGF.xls",sep="\t",col.names=NA,quote=FALSE)

###########################

# Diagnostic plots
topGene <- rownames(res)[which.min(res$padj)]
top.gename<-(gene_ids[(match(topGene,gene_ids$ensembl_gene_id)),3])
png(filename="DiagnosticPlot.png",width=1000,height=1000,res=200)
plotCounts(dds, gene=topGene, main=top.gename, intgroup=c("type"),col="red",pch=16)
dev.off()


#######################
# Bioage genes heatmap
library("genefilter")
library("pheatmap")

load("/home/smm30e/Documents/Scripts/smallRNA/MiRNA/SubRead/DGFspecific/DGF.unique.rda")
rld.DGFspec<-(assay(rld))[DGF.unique,]

topVarGenes <- order(rowMeans(rld.DGFspec),decreasing=TRUE)

mat <- rld.DGFspec[topVarGenes, ]
mat <- mat - rowMeans(mat)

df <- as.data.frame(colData(rld)[,c("type","condition","D.gender","R.gender")])

## separting all samples based on type
# b- pre-perfused kidney transplant patients
# b1 - post-perfused kidney transplant patients

mat.IGF<- mat[,grep("^IGF",colnames(mat))]
mat.DGF<-mat[,grep("^DGF",colnames(mat))]

## ordering the table of matrix to keep IGF and DGF separately
mat.order<-mat[,c(colnames(mat.IGF),colnames(mat.DGF))]

## ordering the table of matrix to keep b and b1 separately
dset.order<-dset.order[,c(colnames(b.all),colnames(b1.all))]

## heatmap for differentially expressed genes
tiff(filename="heatmap.DGFspecific.tiff",width=800,height=1000,res=80,compression = "lzw")
pheatmap(mat.order, annotation_col=df,scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(13),
         annotation_colors = ann_colors)
dev.off()

