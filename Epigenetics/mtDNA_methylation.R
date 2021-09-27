## Script modified : 16/10/2016
## Created: 2016 
## Copyright Suhaib Mohammed, 2016

setwd("/media/smm30e/Backup Data/WGBS/Analysis/MTmethylation/Analysis/")

require(foreach)


path<-"/media/smm30e/Backup Data/WGBS/Analysis/MTmethylation/methloci"


all_files <-list.files(path,pattern="^CpG_context", full.names=T)
dset<-foreach(i=1:length(all_files), .combine =cbind) %do% read.csv(all_files[i], header=TRUE)
colnames(dset)<-c("CpGloc","meth","unmeth","perc","ratioindex")
dset$sample<-gsub("rpt","",dset$sample)
head(dset)

dset<-foreach(i=1:length(all_files), .combine =cbind) %do% read.table(all_files[i], header=TRUE,sep="\t")
rownames(dset)<-paste("MT",dset[,2],sep="")
dset<-dset[,-grep("X|^CpGloc", colnames(dset))]

## extracting sample names associated with the CpG context files
col.names<-gsub("CpG_context_","",grep("^CpG",unlist(strsplit(all_files,"/",fixed=TRUE)),value=TRUE))
col.names<-gsub("rpt","",col.names)

### renaming columns
window=seq(1,length(colnames(dset)),by=4)

for (i in 1:length(col.names))
{
  print(paste("i=",i,",window=",window[i]))
  colnames(dset)[window[i]:(window[i]+3)] <- paste(rep(col.names[i],4),sep="_",colnames(dset)[window[i]:(window[i]+3)])
}

####  GSK patients database 
path.gsk <- "~/Dropbox/Glasgow/GSK/DifferentialGenes/gsk database.csv"
data.gsk <- read.csv(path.gsk, header = TRUE,sep=",")
sample.idx<-sapply(strsplit(as.character(colnames(dset)), split="_"), "[", 1)
samples<-unique(gsub('B.*', '', sample.idx))

patients.db<-(data.gsk[match(samples,data.gsk$labStudy.ID),])
head(patients.db) 


####
#MT1215
## extracting percentage of methyalion from pre samples

plotmtDNAmeth<-function(dset,CpG="MT1321",age=patients.db$D.age){
  
MT.perc<-as.numeric(dset[CpG,intersect(grep("perc",colnames(dset)),grep("B_",colnames(dset)))])

tiff(filename=paste("CpG=",CpG,".tiff",sep=""),width=1000,height=1000,res=200)

plot(MT.perc,age,col="black",ylab="Chronological age (years)", xlab="mtDNAMethylation %",pch=16,main=paste("CpG=",CpG)) 
lm.nuc<-lm(age~MT.perc)
abline(lm.nuc,col="black",lwd=3)

r2.nuc<-round(summary(lm.nuc)$adj.r.squared,digit=4)
nuc.pval<- summary(lm.nuc)$coefficients[2,4]

nuc.rp = vector('expression',2)
nuc.rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                       list(MYVALUE = format(r2.nuc,dig=3)))[2]
nuc.rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                       list(MYOTHERVALUE = format(nuc.pval, digits = 2)))[2]
legend('topright', legend = nuc.rp, bty = 'n',cex=1.25)

dev.off()

}

## hyper methylation
plotmtDNAmeth(dset,CpG="MT1321",age=patients.db$D.age)
plotmtDNAmeth(dset,CpG="MT1176",age=patients.db$D.age)
plotmtDNAmeth(dset,CpG="MT1215",age=patients.db$D.age)
plotmtDNAmeth(dset,CpG="MT1301",age=patients.db$D.age)
plotmtDNAmeth(dset,CpG="MT1215",age=patients.db$D.age)
plotmtDNAmeth(dset,CpG="MT1313",age=patients.db$D.age)


## hypo
plotmtDNAmeth(dset,CpG="MT1225",age=patients.db$D.age)
plotmtDNAmeth(dset,CpG="MT1261",age=patients.db$D.age)


points(MT1313.perc,patients.db$D.age,col="red") 
lm.1313<-lm(patients.db$D.age~MT1313.perc)
abline(lm.1313,col="red",lwd=2)
#########

### percentage chronoage and dnamethylation pre samples

CpG_ids<-as.matrix(read.table(file ="/media/smm30e/Backup Data/WGBS/Analysis/MTmethylation/mtDNA_54CpG.csv", sep="\t", header=FALSE))
pre.perc<-dset[,intersect(grep("perc",colnames(dset)),grep("B_",colnames(dset)))]
dim(pre.perc)
pre.perc.54CpGs<-pre.perc[CpG_ids[,1],]
dim(pre.perc.54CpGs)

pre.perc<-dset[,intersect(grep("perc",colnames(dset)),grep("B_",colnames(dset)))]
nozeropre.perc<-pre.perc[which(apply(pre.perc,1,mean)>0),]
samples.mean.perc<-apply(nozeropre.perc,2,mean)

plotmtDNAmeth2(X=samples.mean.perc,Y=patients.db$D.age,xlab="mtDNAMethylation(%)",ylab="Chronological age (years)",color="blue",main="mtDNAmeth vs Age")

## 
load("/media/smm30e/Backup Data/WGBS/Analysis/MTmethylation/Analysis/Horvath_DNAmAge.rda")

## Horvath
mAge.cdkn2a

## mtDNAmeth and HorvathAge
plotmtDNAmeth2(X=samples.mean.perc,Y=mAge.cdkn2a[,2],xlab="mtDNAMethylation(%)",ylab="Horvath's DNAmAge",color="blue",main="mtDNAmeth vs Horvath's DNAmAge")

plotmtDNAmeth2(X=samples.mean.perc,Y=mAge.cdkn2a[,1],xlab="mtDNAMethylation(%)",ylab="CDKN2A expression",color="blue",main="mtDNAmeth vs CDKN2A")


dim(nozeropre.perc)
samples.mean.perc<-apply(pre.perc.54CpGs,2,mean)
length(samples.mean.perc)


tiff(filename="54CpGscoverage.tiff",width=2000,height=500,res=100)
boxplot(t(pre.perc.54CpGs),col="brown3",las=2,cex=1,ylab="Methylation (%)",xlab="CpGs",ylim=c(0,100))
dev.off()

plotmtDNAmeth2<-function(X,Y,xlab="mtDNAMethylation(%)",ylab="Chronological age (years)",color="blue",main=""){
  
  tiff(filename=paste(main,".tiff",sep=""),width=1000,height=1000,res=200)
  
  plot(X,Y,col="black",ylab=ylab, xlab=xlab,pch=16,main=paste(main)) 
  lm.nuc<-lm(Y~X)
  abline(lm.nuc,col=color,lwd=3)
  
  r2.nuc<-round(summary(lm.nuc)$adj.r.squared,digit=4)
  nuc.pval<- summary(lm.nuc)$coefficients[2,4]
  
  nuc.rp = vector('expression',2)
  nuc.rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                         list(MYVALUE = format(r2.nuc,dig=3)))[2]
  nuc.rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                         list(MYOTHERVALUE = format(nuc.pval, digits = 2)))[2]
  legend('topleft', legend = nuc.rp, bty = 'n',cex=1.25)
  
  dev.off()
  
}

  
### starifying DGF and IGF
pre.perc.IGF<-dset[,intersect(intersect(grep("perc",colnames(dset)),grep("B_",colnames(dset))),grep("_NoDGF",colnames(dset)))]
pre.perc.IGF.mean<-apply(pre.perc.IGF,2,mean)
pre.perc.DGF<-dset[,intersect(intersect(grep("perc",colnames(dset)),grep("B_",colnames(dset))),grep("_DGF",colnames(dset)))]
pre.perc.DGF.mean<-apply(pre.perc.DGF,2,mean)

##
patients.age<-function(db,IGF){
  sample.idx<-sapply(strsplit(as.character(colnames(IGF)), split="_"), "[", 1)
  samples<-unique(gsub('B.*', '', sample.idx))
  patients.db<-(data.gsk[match(samples,data.gsk$labStudy.ID),])
  return(patients.db) 
}

## DGF IGF age
chron.age.IGF<-patients.age(db=data.gsk,IGF=pre.perc.IGF)
chron.age.DGF<-patients.age(db=data.gsk,IGF=pre.perc.DGF)

##
mAge.cdkn2a.IGF<-mAge.cdkn2a[grep("NoDGF",rownames(mAge.cdkn2a)),]
mAge.cdkn2a.DGF<-mAge.cdkn2a[grep("^DGF_",rownames(mAge.cdkn2a)),]


plotmtDNAmeth3(X=pre.perc.IGF.mean,Y=chron.age.IGF$D.age,A=pre.perc.DGF.mean,B=chron.age.DGF$D.age,
               xlab="mtDNAMethylation(%)", ylab="Chronological age (years)", color=c("blue","red"), main="mtDNAmeth vs Age (IGF-DGF)")
  
plotmtDNAmeth3(X=pre.perc.IGF.mean,Y=mAge.cdkn2a.IGF[,2],A=pre.perc.DGF.mean,B=mAge.cdkn2a.DGF[,2],
               xlab="mtDNAMethylation(%)", ylab="Horvath's DNAmAge", color=c("blue","red"), main="mtDNAmeth vs Horvath's DNAmAge (IGF-DGF)")

plotmtDNAmeth3(X=pre.perc.IGF.mean,Y=mAge.cdkn2a.IGF[,1],A=pre.perc.DGF.mean,B=mAge.cdkn2a.DGF[,1],
               xlab="mtDNAMethylation(%)", ylab="CDKN2A expression", color=c("blue","red"), main="mtDNAmeth vs CDKN2A (IGF-DGF)")

plotmtDNAmeth3<-function(X,Y,A,B,xlab="mtDNAMethylation(%)",ylab="Chronological age (years)",color,main=""){
  
  tiff(filename=paste(main,".tiff",sep=""),width=1000,height=1000,res=140)
  ## IGF
  plot(X,Y,col="blue",ylab=ylab, xlab=xlab,xlim=c(min(X,A),max(X,A)),ylim=c(min(Y,B),max(Y,B)), pch=16,main=paste(main)) 
  lm.nuc1<-lm(Y~X)
  abline(lm.nuc1,col=color[1],lwd=3)
  
  r2.nuc1<-round(summary(lm.nuc1)$adj.r.squared,digit=4)
  nuc.pval1<- summary(lm.nuc1)$coefficients[2,4]
  
  ## DGF
  points(A,B,col="red",pch=8) 
  lm.nuc2<-lm(B~A)
  abline(lm.nuc2,col=color[2],lwd=3)
  
  r2.nuc2<-round(summary(lm.nuc2)$adj.r.squared,digit=4)
  nuc.pval2<- summary(lm.nuc2)$coefficients[2,4]
  
  legend("topright",c(paste("IGF"," p=", round(nuc.pval1,digit=2),sep=""),
                      paste("DGF"," p=",round(nuc.pval2,digit=2),sep="")),
     pch=c(16,8), col=color,lwd=1.3)
  
  dev.off()
}


lm.1215<-lm(patients.db$D.age~samples.mean.perc)
summary(lm(patients.db$D.age~samples.mean.perc))
abline(lm.1215,col="blue",lwd=2)


## fraction of methyaltion calculator
meth<-dset[,intersect(grep("_meth",colnames(dset)),grep("B_",colnames(dset)))]
unmeth<-dset[,intersect(grep("_unmeth",colnames(dset)),grep("B_",colnames(dset)))]

samples.mean.meth<-apply(meth,2,mean)
samples.mean.unmeth<-apply(unmeth,2,mean)


ratio<-samples.mean.meth/samples.mean.unmeth
plot(ratio,patients.db$D.age,col="blue") 
lm.1215<-lm(patients.db$D.age~ratio)
abline(lm.1215,col="blue",lwd=2)


#####################################
## function to convert df to list gene wise
dftolist<-function(dset,row){
  t.dset<-t(dset[row,])
  sample.idx<-paste(sapply(strsplit(as.character(colnames(dset)), split="_"), "[", 1),sapply(strsplit(as.character(colnames(dset)), split="_"), "[", 2),sep="_")
  
  df<- as.data.frame(matrix(nrow=length(unique(sample.idx)),ncol=5,dimnames=list(c(paste(unique(sample.idx))),c("meth","unmeth","ratio","ratioindex","gender"))))
  df$meth<-t.dset[grep("_meth",rownames(t.dset)),]
  df$unmeth<-t.dset[grep("_unmeth",rownames(t.dset)),]
  df$ratio<-t.dset[grep("_ratio$|_ratio.[0-9]", rownames(t.dset)),]
  df$ratioindex<-t.dset[grep("_ratioindex",rownames(t.dset)),]

  
  ## extracting sample ids
  samples<-(gsub('B.*', '', rownames(df))) 
  df$gender<-(data.gsk[match(samples,data.gsk$labStudy.ID),]$D.sex.0.M)
  
  ## converting dataframe to list
  df.list<-list(df)
  names(df.list)<-colnames(t.dset)
  return(as.list(df.list))
}
##########################################

gene.list<-list()

for (i in 1:nrow(dset))
{
  gene.list[[i]]<-dftolist(dset,i)
}

###########

#####################################
kruskal.sig.test<-function(list)
{
  Y<-as.matrix(unlist(list))
  Data=data.frame(val=Y[,1],freq=gsub("[0-9]","",rownames(Y)))
  res <- kruskal.test(val~freq, data=Data)
  print("#######################")
  print(res)
  print("#######################")
  return(res$p.val)
}

#####################################
## fucntion to add sig stars
require(pander)
sig.stars<-function(data)
{ 
  data<-paste(data,gsub("(\\_)|\\s","",add.significance.stars(data)),sep="")
  return(data)
}
######################################  
#### function to plot boxplots and test for significance
plot.boxplot<-function(list,color,main,ylab,ylim)
{
  a<-boxplot(list,col=color,las=1,ylab=paste(ylab),
             main=paste(main,sep=""),cex.axis=1,cex.lab=1,ylim=ylim)
 #text(1:length(a$n), a$stats[1,]-500, paste("n =",a$n))
  legend('topright',legend=paste("p-val =",sig.stars(round(kruskal.sig.test(list),digit=3))) ,bty = "n",cex=1)
  n.pval<- data.frame(n=length(unlist(list)),p.val=sig.stars(round(kruskal.sig.test(list),digit=4)))
  return(as.character(n.pval$p.val))
}

########
## all
## DGF NoDGF
########
## all
## DGF NoDGF
strat.samples<-function(dset,type,condition,column,scale){
  if (condition=="all"){
    strat<-dset[grep(paste("_",type,"$",sep=""),rownames(dset)),column]/scale
    print(paste("stratified samples="))
    print(dset[grep(paste("_",type,"$",sep=""),rownames(dset)),])
  }
  else if (condition=="B" || condition=="B1"){
    strat<-dset[intersect(grep(paste("_",type,"$",sep=""),rownames(dset)),grep(paste(condition,"_",sep=""),rownames(dset))),column]/scale
    print(paste("stratified samples="))
    print(dset[intersect(grep(paste("_",type,"$",sep=""),rownames(dset)),grep(paste(condition,"_",sep=""),rownames(dset))),])
  }
  return(strat)
}

##
## function to remove names of list and convert to data frame
dfnull<-function(list){
  names(list)<-NULL
  df.list<-as.data.frame(list)
  return(df.list)
}

#######################


Promoter=FALSE
Intragenic=TRUE
scale=1
gene.mat<-as.data.frame(matrix(nrow=nrow(dset),ncol=3,dimnames=list(c(rownames(dset)),c("all","pre","post"))))

for (i in 1:length(gene.list)) {
  print(names(gene.list[[i]]))
  rownames(gene.mat)[i]<-names(gene.list[[i]])
  
  dfnull(gene.list[[i]])
  ############
  ## all
  # DGF NoDGF
  DGF.all.ratio<-strat.samples(dfnull(gene.list[[i]]),type="DGF",condition="all",column="meth",scale=scale)
  NoDGF.all.ratio<-strat.samples(dfnull(gene.list[[i]]),type="NoDGF",conditio="all",column="meth",scale=scale)
  DGF.NoDGF.all.list<-list(DGF=DGF.all.ratio,NoDGF=NoDGF.all.ratio)
  
  ## startify Pre
  # DGF NoDGF
  DGF.B.ratio<-strat.samples(dfnull(gene.list[[i]]),type="DGF",condition="B",column="meth",scale=scale)
  NoDGF.B.ratio<-strat.samples(dfnull(gene.list[[i]]),type="NoDGF",condition="B",column="meth",scale=scale)
  DGF.NoDGF.B.list<-list(DGF.Pre=DGF.B.ratio,NoDGF.Pre=NoDGF.B.ratio)
  
  ## Stratify Post
  # DGF NoDGF
  DGF.B1.ratio<-strat.samples(dfnull(gene.list[[i]]),type="DGF",condition="B1",column="meth",scale=scale)
  NoDGF.B1.ratio<-strat.samples(dfnull(gene.list[[i]]),type="NoDGF",condition="B1",column="meth",scale=scale)
  DGF.NoDGF.B1.list<-list(DGF.Post=DGF.B1.ratio,NoDGF.Post=NoDGF.B1.ratio)
  
  
  ## y axis scale edges adjustment 
  if ((max(unlist(DGF.NoDGF.all.list))-min(unlist(DGF.NoDGF.all.list)))>10^5){
    edge<-20^4
  }
  else if ((max(unlist(DGF.NoDGF.all.list))-min(unlist(DGF.NoDGF.all.list)))>10^4){
    edge<-20^3
  }
  else if ((max(unlist(DGF.NoDGF.all.list))-min(unlist(DGF.NoDGF.all.list)))>10^3){
    edge<-20^2
  }
  else {edge<-20}
  
  ### for promoters
  if (Promoter) {
    # pdf("Phosphate_food_consumption.pdf",width=8,height=10)
    png(paste("PromoterMethyl_",names(gene.list[[i]]),".png",sep=""),width=1500,height=700,res=180)
    par(mfrow=c(1,3))
    par(mar=c(5, 5, 4.1, 1))
    gene.mat$all[i]<-plot.boxplot(DGF.NoDGF.all.list,color=c("red","grey"),main=paste(names(gene.list[[i]]),"- all"),
                                  ylim=c(min(unlist(DGF.NoDGF.all.list))-edge,max(unlist(DGF.NoDGF.all.list))+edge),ylab=paste("Promoter methylation (mCpGs)"))
    gene.mat$pre[i]<-plot.boxplot(DGF.NoDGF.B.list,color=c("red","grey"),main=paste(names(gene.list[[i]]),"- Pre"),
                                  ylim=c(min(unlist(DGF.NoDGF.all.list))-edge,max(unlist(DGF.NoDGF.all.list))+edge),ylab=paste("Promoter methylation (mCpGs)"))
    gene.mat$post[i]<-plot.boxplot(DGF.NoDGF.B1.list,color=c("red","grey"),main=paste(names(gene.list[[i]]),"- Post"),
                                   ylim=c(min(unlist(DGF.NoDGF.all.list))-edge,max(unlist(DGF.NoDGF.all.list))+edge),ylab=paste("Promoter methylation (mCpGs)"))
    
    dev.off()
  }
  
  ### for intragenic
  if (Intragenic) {
    png(paste("IntraMethyl_",names(gene.list[[i]]),".png",sep=""),width=1500,height=700,res=180)
    par(mfrow=c(1,3))
    par(mar=c(5, 5, 4.1, 1))
          
    gene.mat$all[i]<- plot.boxplot(DGF.NoDGF.all.list,color=c("red","grey"),main=paste(names(gene.list[[i]]),"- all"),
                                   ylim=c(min(unlist(DGF.NoDGF.all.list))-edge,max(unlist(DGF.NoDGF.all.list))+edge),ylab=paste("Intragenic methylation (mCpGs)"))
    gene.mat$pre[i]<-plot.boxplot(DGF.NoDGF.B.list,color=c("red","grey"),main=paste(names(gene.list[[i]]),"- Pre"),
                                  ylim=c(min(unlist(DGF.NoDGF.all.list))-edge,max(unlist(DGF.NoDGF.all.list))+edge),ylab=paste("Intragenic methylation (mCpGs)"))
    gene.mat$post[i]<- plot.boxplot(DGF.NoDGF.B1.list,color=c("red","grey"),main=paste(names(gene.list[[i]]),"- Post"),
                                    ylim=c(min(unlist(DGF.NoDGF.all.list))-edge,max(unlist(DGF.NoDGF.all.list))+edge),ylab=paste("Intragenic methylation (mCpGs)"))
    
    dev.off()
  }
  
}


