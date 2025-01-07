args=commandArgs(T)
#data augmentation for spatial RNA sequencing data
library(Seurat)
library("parallel")
library("foreach")
library(reshape2)
library(tidyverse)

dir=args[1]
samid=args[2]
#augment_expmean=as.numeric(args[3])
samdir=paste0(dir,"/",samid,"/");#/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/151673"
#samid="151673"
augment_expmean=1

rdsfile=paste0(samdir,"/",samid,".rds");
seurat=readRDS(rdsfile)
#ground_truth=seurat@meta.data$cell_type #cellxgene_breast

aug.file=paste0(samdir,"/augment.counts.csv")
#if (file.exists(aug.file)){
#  stop("aug.csv file.exists")
#}

counts <- seurat@assays$RNA@counts

sta_counts <- function(sub){
  sub=sub[rowSums(sub)>5,]
  sparsity=mean(apply(sub,1,function(x) length(x[x==0])/length(x)))
  sta=c(samid,ncol(sub),nrow(sub),sparsity,mean(colMeans(sub)),mean(apply(sub,1,function(x) sd(x)^2)))
  sta=data.frame(sta);colnames(sta)=samid;
  row.names(sta)=c("samid","spot.num","gene.num","sparsity","exp.mean","exp.var")
  return (sta)
}
sam.sta=sta_counts(counts); exp.mean=as.numeric(sam.sta["exp.mean",])

#calculate the augment ratio 
aug.ratio=1
if (exp.mean<augment_expmean){aug.ratio=round(augment_expmean/exp.mean,0)+2}

#dist between each single cell

cells=colnames(counts)
genes=rownames(counts)

cell.exp.augment <- function(cellid){
  cell.neighbours=names(cell.dist[cellid,][order(cell.dist[cellid,])[1:aug.ratio]])
  aug.cell.exp <- apply(counts[,cell.neighbours],1,function(x) sum(x))
  return(aug.cell.exp)
}

if (aug.ratio>1){
  cell.dist <- as.matrix(dist(t(as.matrix(counts))))
  aug.cell.exp.all<-mclapply(cells, cell.exp.augment, mc.cores = 20)
  #counts_augment =counts;
  counts_augment = array( unlist( aug.cell.exp.all ) , dim = c( nrow(counts) , ncol(counts) ))
  #for (i in 1:length(cells)){
  #  counts_augment[,i]=aug.cell.exp.all[[i]]
  #}
  sam.sta.aug=sta_counts(counts_augment);
  colnames(sam.sta.aug)=paste0("aug.size",aug.ratio)
  sam.sta=cbind(sam.sta,sam.sta.aug)
  write.csv(counts_augment,file=paste0(samdir,"/augment.counts.csv"))
}else{
  sam.sta.aug=sam.sta
  colnames(sam.sta.aug)=paste0("aug.size",aug.ratio)
  sam.sta=cbind(sam.sta,sam.sta.aug)
}
sam.sta=t(sam.sta)
write.table(sam.sta,file=paste0(samdir,"/sam.exp.sta"),sep="\t",col.names = T,row.names = T,quote=F)

