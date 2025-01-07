args=commandArgs(T)
#data augmentation for spatial RNA sequencing data
library(Seurat)
library("scran")
library(ggplot2)
library(Seurat)
library(tidyverse)
library(cluster)
library("SingleCellExperiment")
library("BASiCS")
library("splatter")
library("monocle")

dir=args[1]
samid=args[2]
#augment_expmean=as.numeric(args[3])
samdir=paste0(dir,"/",samid,"/");#/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/151673"
#samid="151673"
#augment_expmean=1
aug.counts.file=paste0(samdir,"/augment.counts.csv")
wasfile=paste0(samdir,"/gene.Wasserstein.distance.sta")
sam.sta.file=aug.counts.file=paste0(samdir,"/sam.exp.sta")
rdsfile=paste0(samdir,"/",samid,".rds");
seurat=readRDS(rdsfile)
ground_truth <-seurat@meta.data$ground_truth
seurat <- NormalizeData(seurat, normalization.method = 'LogNormalize',scale.factor = 10000,verbose = FALSE)
seurat <- NormalizeData(seurat, verbose = FALSE, assay = "Spatial")

seurat <-  FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 20000)
counts <- seurat@assays$Spatial@counts
sce <- SingleCellExperiment(assays=list(counts=counts))
wasdata= read.table(wasfile,sep="\t",header = T)


#准备BASICS的数据
Chain <- BASiCS_MCMC(sce, N = 10, Thin = 2, Burn = 2, Regression = FALSE,WithSpikes = FALSE)
# Highly and lowly variable genes detection (within a single group of cells)
BASiCS_HVG <- BASiCS_DetectHVG(Chain, VarThreshold = 0.60, EFDR = 0.10, Plot = TRUE)
BASiCS.sort.gene  = BASiCS_HVG@Table[order(BASiCS_HVG@Table$Sigma,decreasing = T),]$GeneName

#准备monocle的数据
gene_ann <- data.frame(
  gene_short_name = row.names(counts),
  row.names = row.names(counts)
)
sample_ann <- data.frame(
  sample_short_name = colnames(counts),
  col.names = colnames(counts)
)
pd <- new("AnnotatedDataFrame",
          data=sample_ann)
rownames(pd)=colnames(counts)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
sparse_counts <- Matrix(counts, sparse = TRUE)
monocle_cds<- newCellDataSet(sparse_counts,phenoData = pd,featureData =fd,expressionFamily = negbinomial.size(),lowerDetectionLimit=1)

monocle_cds <- estimateSizeFactors(monocle_cds) #估计每个样本的size factor
monocle_cds <- estimateDispersions(monocle_cds)

#Removing 270 outliers
disp_table <- dispersionTable(monocle_cds)
disp.genes <- subset(disp_table,mean_expression >=0.01 & dispersion_empirical>=1* dispersion_fit)$gene_id
monocle_cds <- setOrderingFilter(monocle_cds,disp.genes)
#unsup_clustering_genes <- subset(disp_table, mean_expression = 0.1) #基本没过滤 最开始的做法
unsup_clustering_genes <- subset(disp_table, mean_expression >=0.01 & dispersion_empirical>=1* dispersion_fit)
#monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id[1:100]) #function marks genes that will be used for clustering in subsequent calls to clusterCells. The
#这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将他嵌入cds对象
#这些基因被存储在monocle_cds@featureData@data[['use_for_ordering']]；可以通过table(monocle_cds@featureData@data[['use_for_ordering']])查看
#plot_ordering_genes(monocle_cds)
disp_table$dispersion_diff=disp_table$dispersion_empirical - disp_table$dispersion_fit
#筛选基因后，按照表达量大小对基因排序
monocle_hvg <- unsup_clustering_genes[order(unsup_clustering_genes$mean_expression, decreasing=TRUE),][,1] #原来的方法相当于直接用基因的表达量从高到低排序
##准备scran的数据
scran_sce <- sce
scran_sce <- computeSumFactors(scran_sce)
scran_sce  <- logNormCounts(scran_sce )
allf <- modelGeneVar(scran_sce)


washvg=read.table(wasfile,sep="\t",header=T)
wasvalue=washvg$KLD_value
if ( !file.exists(wasfile)){stop("was file not exists")}

Idents(seurat) <- ground_truth
cluster_tag <- "cell_type"
clusters <- unique(ground_truth);

clusternum <- table(ground_truth)[order(table(ground_truth),decreasing = T)]
clusters=names(clusternum[clusternum>20])

DEG_result <- data.frame(matrix(0,nrow(seurat@assays$Spatial@data),length(clusters)))
row.names(DEG_result)<- row.names(seurat@assays$Spatial@data)
colnames(DEG_result) <- paste0("C_",clusters)
DEG_result_value <- DEG_result

sample_id=1
cluster.gene.list=list()
for (i in 1:length(clusters)){
  cid=clusters[i]
  DEG=FindMarkers(seurat,ident.1 = cid,slot="data",min.cells.feature =3,logfc.threshold = 0.25,min.pct = 0.25,test.use = "wilcox",min.cells.group = 3)
  upgene <- row.names(DEG[DEG$avg_log2FC>0.25 & DEG$p_val<1e-10,])
  downgene <- row.names(DEG[DEG$avg_log2FC< (-0.25) & DEG$p_val<1e-10,])
  DEG_result[upgene,i]=1
  DEG_result[downgene,i]=-1
  DEG_result_value[upgene,i]<-DEG[DEG$avg_log2FC>0.25 & DEG$p_val<1e-10,]$avg_log2FC
  DEG_result_value[downgene,i]<-DEG[DEG$avg_log2FC< (-0.25) & DEG$p_val<1e-10,]$avg_log2FC
  cluster.gene.list[[i]]=c(upgene,downgene)
}


degs= cluster.gene.list[[1]]
for(i in 2:length(clusters)){
  degs=c(degs,  cluster.gene.list[[i]])
}
degs=unique(degs)  #所有差异基因的lists



bins=c(100,250,500,1000,1500,2000)
staHV=data.frame(matrix(NA,length(bins),12))
row.names(staHV)=bins
methods=c("WAS","seurat","scran","Monocle","BASICS")
colnames(staHV)=c("sample","topgene",paste0(methods,"_overlap"),paste0(methods,"_coverage"))


for (binid in 1:length(bins)){
  #选定top基因的数目
  HVNUM=bins[binid]
  #获取Was方法的top基因
  mincutoff <- quantile(wasvalue,1-HVNUM/length(wasvalue))
  washvgsub <- washvg[wasvalue>mincutoff,]; 
  gene_was=washvgsub$geneid
  #gene_was="NA"
  #获取scran方法的top基因
  gene_scran <- getTopHVGs(allf, n=HVNUM)
  #获取monocle方法的top基因
  gene_monocle <-  monocle_hvg[1:HVNUM]
  #获取seurat方法的top基因
  gene_seurat <- VariableFeatures(seurat)[1:HVNUM]
  gene_BASICS <- BASiCS.sort.gene[1:HVNUM]
  sta=data.frame(matrix(NA,length(clusters),12))
  row.names(sta)=clusters
  colnames(sta)=c("topgene","clusters",paste0(methods,"_overlap"),paste0(methods,"_coverage"))
  #统计类别的差异表达基因被各种HV基因覆盖的基因数和比例
  for (i in 1:length(clusters)){
    cid=clusters[i]
    deg=cluster.gene.list[[i]];degnum=length(deg);degnum=degnum+degnum;
    all=c(gene_was,deg);tnum=length(all);overnum=length(table(c(gene_was,deg))[table(c(gene_was,deg))>1]);overatio=overnum/tnum;cover=overnum/length(deg)
    sta[i,1:2]=c(HVNUM,cid);sta[i,c(3,8)]=c(overatio,cover)
    all=c(gene_seurat,deg);tnum=length(all);overnum=length(table(c(gene_seurat,deg))[table(c(gene_seurat,deg))>1]);overatio=overnum/tnum;cover=overnum/length(deg)
    sta[i,c(4,9)]=c(overatio,cover)
    all=c(gene_scran,deg);tnum=length(all);overnum=length(table(c(gene_scran,deg))[table(c(gene_scran,deg))>1]);overatio=overnum/tnum;cover=overnum/length(deg)
    sta[i,c(5,10)]=c(overatio,cover)
    all=c(gene_monocle,deg);tnum=length(all);overnum=length(table(c(gene_monocle,deg))[table(c(gene_monocle,deg))>1]);overatio=overnum/tnum;cover=overnum/length(deg)
    sta[i,c(6,11)]=c(overatio,cover)
    all=c(gene_BASICS,deg);tnum=length(all);overnum=length(table(c(gene_BASICS,deg))[table(c(gene_BASICS,deg))>1]);overatio=overnum/tnum;cover=overnum/length(deg)
    sta[i,c(7,12)]=c(overatio,cover)
    if (i==1){
      degs=deg;wass=gene_was;seurats=gene_seurat;monocles=gene_monocle;scrans=gene_scran;basics=gene_BASICS
    }else{
      degs=c(degs,deg);wass=c(wass,gene_was);seurats=c(seurats,gene_seurat);
      monocles=c(monocles,gene_monocle);scrans=c(scrans,gene_scran);basics=c(basics,gene_BASICS)
    }
  }
  degs=unique(degs);wass=unique(wass);seurats=unique(seurats);monocles=unique(monocles);scrans=unique(scrans);basics=unique(basics)
  degnum <- length(degs);wasnum=length(wass);seuratnum=length(seurats);monoclenum=length(monocles);scrannum=length(scrans);basicsnum=length(basics)
  wasnum2=length(unique(c(degs,wass)));seuratnum2=length(unique(c(degs,seurats)));monoclenum2=length(unique(c(degs,monocles)));scrannum2=length(unique(c(degs,scrans)));basicsnum2=length(unique(c(degs,basics)));
  wasnum3=length(table(c(degs,wass))[table(c(degs,wass))>1]);seuratnum3=length(table(c(degs,seurats))[table(c(degs,seurats))>1]);
  monoclenum3=length(table(c(degs,monocles))[table(c(degs,monocles))>1]);scrannum3=length(table(c(degs,scrans))[table(c(degs,scrans))>1]);basicsnum3=length(table(c(degs,basics))[table(c(degs,basics))>1]);
  wasover=wasnum3/wasnum2;seuratover=seuratnum3/seuratnum2;scranover=scrannum3/scrannum2;monoover=monoclenum3/monoclenum2;basicsover=basicsnum3/basicsnum2;
  wascover=wasnum3/degnum;seuratcover=seuratnum3/degnum;scrancover=scrannum3/degnum;monocover=monoclenum3/degnum;basicscover=basicsnum3/degnum;
  staHV[binid,] <- c(samid,HVNUM,wasover,seuratover,scranover,monoover,basicsover,wascover,seuratcover,monocover,scrancover,basicscover)
  sta=cbind(samid,sta)
  if (sample_id==1 & binid==1 ){
    staall=sta
  }else{
    staall=rbind(staall,sta)
  }
}
write.table(staall,file=paste0(samdir,"/HVgene.DEG.recall.tmp.sta"),sep="\t",row.names = F,col.names = T)
write.table(staHV,file=paste0(samdir,"/HVgene.DEG.recall.sta"),sep="\t",row.names = F,col.names = T)


############################DEG基因和非DEG基因在不同指标的分布差异￥￥￥￥￥￥￥￥￥￥￥

mainsize=1.5
plot_was <- function(maintext){
  indata=wasdata
  valuedeg=indata[indata$geneid %in% degs,]$KLD_value
  valueother=indata[!(indata$geneid %in% degs),]$KLD_value
  plot(density(indata$KLD_value),xlim=c(0,quantile(indata$KLD_value,0.998)),col="black",main="WASHV",cex.main=mainsize)
  lines(density(valuedeg),col="red")
  lines(density(valueother),col="green")
  legend("topright",legend = c("DEG","Not DEG","All"),col=c("red","green","black"),lwd=c(1,1,1),cex=1)
}

plot_seurat <- function(maintext){
  indata=seurat@assays$Spatial@meta.features
  valuedeg=indata[row.names(seurat@assays$Spatial@meta.features) %in% degs,]$vst.variance.standardized
  valueother=indata[!(row.names(seurat@assays$Spatial@meta.features) %in% degs),]$vst.variance.standardized
  plot(density(indata$vst.variance.standardized),xlim=c(0,quantile(indata$vst.variance.standardized,0.998)),col="black",main="seurat",cex.main=mainsize)
  lines(density(valuedeg),col="red")
  lines(density(valueother),col="green")
  legend("topright",legend = c("DEG","Not DEG","All"),col=c("red","green","black"),lwd=c(1,1,1),cex=1)
}

plot_scran <- function(maintext){
  allf_ori=allf
  allf <- allf[allf$bio>=quantile(allf$bio,0.001) & allf$bio <= quantile(allf$bio,0.999), ]
  valuedeg=allf[row.names(allf) %in% degs,]$bio
  valueother=allf[!(row.names(allf) %in% degs),]$bio
  plot(density(allf$bio),xlim=c(quantile(allf$bio,0.001),quantile(allf$bio,0.999)),col="black",main="scran",cex.main=mainsize)
  lines(density(valuedeg),col="red")
  lines(density(valueother),col="green")
  legend("topright",legend = c("DEG","Not DEG","All"),col=c("red","green","black"),lwd=c(1,1,1),cex=1)
}

plot_monocle <- function(maintext){
  #disp_table$dispersion_diff=disp_table$dispersion_empirical - disp_table$dispersion_fit
  #筛选基因后，按照表达量大小对基因排序
  #monocle_hvg <- unsup_clustering_genes[order(unsup_clustering_genes$mean_expression, decreasing=TRUE),][,1] #原来的方法相当于直接用基因的表达量从高到低排序
  ##准备scran的数据
  valuedeg=unsup_clustering_genes[unsup_clustering_genes$gene_id %in% degs,]$mean_expression
  valueother=unsup_clustering_genes[!(unsup_clustering_genes$gene_id %in% degs),]$mean_expression
  plot(density(unsup_clustering_genes$mean_expression),xlim=c(quantile(unsup_clustering_genes$mean_expression,0.01),
                                                              quantile(unsup_clustering_genes$mean_expression,0.999)),col="black",main="monocle",cex.main=mainsize)
  lines(density(valuedeg),col="red")
  lines(density(valueother),col="green")
  legend("topright",legend = c("DEG","Not DEG","All"),col=c("red","green","black"),lwd=c(1,1,1),cex=1)
}

plot_BASICS <- function(maintext){
  #disp_table$dispersion_diff=disp_table$dispersion_empirical - disp_table$dispersion_fit
  #筛选基因后，按照表达量大小对基因排序
  #monocle_hvg <- unsup_clustering_genes[order(unsup_clustering_genes$mean_expression, decreasing=TRUE),][,1] #原来的方法相当于直接用基因的表达量从高到低排序
  ##准备scran的数据
  # BASiCS.sort.gene  = BASiCS_HVG@Table[order(BASiCS_HVG@Table$Sigma,decreasing = T),]$GeneName
  allf=BASiCS_HVG@Table
  allf <- allf[allf$Sigma>=quantile(allf$Sigma,0.001) & allf$Sigma <= quantile(allf$Sigma,0.999), ]
  valuedeg=allf[allf$GeneName %in% degs,]$Sigma
  valueother=allf[!(allf$GeneName %in% degs),]$Sigma
  plot(density(allf$Sigma),xlim=c(quantile(allf$Sigma,0.001),quantile(allf$Sigma,0.999)),col="black",main="BASICS",cex.main=mainsize)
  lines(density(valuedeg),col="red")
  lines(density(valueother),col="green")
  legend("topright",legend = c("DEG","Not DEG","All"),col=c("red","green","black"),lwd=c(1,1,1),cex=1)
}
pdf.file=paste0(samdir,"/HVgene.DEG.density.compare",samid,".pdf")
pdf(pdf.file,width = 15,height=4)
par(mfrow=c(1,5))
plot_was("WAS");
plot_seurat("Seurat")
plot_scran("Scran")
plot_monocle("Moncle")
plot_BASICS("BASICS")
dev.off()
write.table(staall,file=paste0(samdir,"/HVgene.DEG.recall.tmp.sta"),sep="\t",row.names = F,col.names = T)
write.table(staHV,file=paste0(samdir,"/HVgene.DEG.recall.sta"),sep="\t",row.names = F,col.names = T)






