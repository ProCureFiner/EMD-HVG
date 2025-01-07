args=commandArgs(T)
#
gene_num= 20000; #as.numeric(args[1])

library(scran)
library(scuttle)
library("edgeR")
library("DESeq2")

library(ggplot2)
library(Seurat)
library(SpatialExperiment)
library(spatialLIBD)
library(scater)
library(ggspavis)
library(SCENIC)
library(SCopeLoomR)
library(foreach)
library(tidyverse)
library(cluster)
library(cowplot)
library(ggpubr)
library("reshape2")
library("transport")

library("splatter")

set.seed(1)
sce <- mockSCE()
out_path="/media/disk1/SYSU/project/splatter_simulation/"
setwd(out_path)

# Estimate parameters from mock data 利用已有的数据估计模拟的参数；这里测试151673空间转录组的数据用来模拟数据的效果
# Estimate parameters from mock data 利用已有的数据估计模拟的参数；这里测试151673空间转录组的数据用来模拟数据的效果
param_celln=2000
param_clusters=2
param_genen=10000
param_degr=0.2
param_facloc=0.2
param_facscale=0.2
param_libloc=8
param_libscale=0.2
cells=c(2000,4000,6000,8000);groups=c(2,4,6,8);locs=c(7,7.6,8,8.3,8.7,9);scales=c(0.1,0.5,1);degrs=c(0.1,0.2,0.3);

cell_clusterdir=paste0(out_path,"/",gene_num,"gene")
if (! file.exists(cell_clusterdir)) {
  dir.create(cell_clusterdir)
}

for(cellbin in 1:length(cells)){
#cellbin=1
param_clusters=groups[cellbin]
param_celln=cells[cellbin]

# 检查是否成功创建了目录

samnum=length(locs)*length(scales)*length(degrs)
seurat_sta = data.frame(matrix(NA,samnum,12))
colnames(seurat_sta)=c("sample","libloc","libscale","degr","celltype_num","spot_num","spot_sum","spot_mean","gene_num","gene_sum","gene_mean","gene_var")
idsort=0;
for (i in 1:length(locs)){
  param_libloc=locs[i];sum=0;print(i);
  for (j in 1:length(scales)){
    param_libscale=scales[j];
    for (k in 1:length(degrs)){
      idsort=idsort+1;
      print(idsort)
      param_degr=degrs[k]
      sample=paste0("splatter","_",gene_num,"_",param_celln,"_",param_clusters,"_",param_libloc,"_",param_libscale,"_",param_degr)
      samdir=paste0(cell_clusterdir,"/",sample,"/")
      if (! dir.exists(samdir)) {
        dir.create(samdir)
      }
      sim.groups <- splatSimulate(group.prob =rep(round(1/param_clusters,3),param_clusters), method = "groups", de.prob = param_degr, verbose = FALSE,nGenes = param_genen,batchCells = param_celln,de.facLoc = param_facloc, de.facScale = param_facloc,lib.loc=param_libloc,lib.scale=param_libscale)
      sim.groups <- logNormCounts(sim.groups)
      sim.groups <- runPCA(sim.groups)
      plotPCA(sim.groups, colour_by = "Group")
      
      counts <- sim.groups@assays@data$counts
      #spotsum=colSums(counts)
      #counts <- counts[order(rowSums(counts),decreasing = T),]
      #genesta <- data.frame(matrix(NA,nrow(counts),4))
      #colnames(genesta)=c("geneid","genesum","genemean","genevar")
      #genes=rownames(counts);genesta[,1]=genes
      #genesum=rowSums(counts);genesta[,2]=genesum
      #genemean=genesum/ncol(counts);genesta[,3]=genemean
      #genevar=apply(counts,1,function(x) sd(x)^2);genesta[,4]=genevar
      #genesta <- genesta[genesta$genemean>=quantile(genesta$genemean,0.001) & genesta$genemean<=quantile(genesta$genemean,0.999),]
      seurat <- CreateSeuratObject(counts = counts, project = sample,assay = "RNA")
      meta.data <- data.frame(sim.groups@colData)
      seurat@meta.data=meta.data
      ground_truth=seurat@meta.data$Group
      seurat@meta.data=cbind(seurat@meta.data,ground_truth)
      saveRDS(seurat,file=paste0(samdir,"/",sample,".rds"))
    }
  }
}
}


