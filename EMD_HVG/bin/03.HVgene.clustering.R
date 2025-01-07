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

calDBI <- function(x=data,labels=labesls)
  ##DBI：任意两类别的类内样本到类中心平均距离之和除以两类中心点之间的距离，取最大值。DBI越小意味着类内距离越小，同时类间距离越大。
  ##data必须行为样本，列为特征 戴维森堡丁指数（DBI）
{
  labels2=labels
  labels=labels2[labels2 %in% names(table(labels2))[table(labels)>2]]
  clusters_n <- length(unique(labels))
  cluster_k <- list()
  for (i in c(1:clusters_n)) {
    cluster_k[[i]] <- x[which(labels==i),]
  }
  
  centroids <- list()
  for (i in c(1:clusters_n)) {
    centroids[[i]] <- apply(cluster_k[[i]],2,mean)
  }
  
  s <- list()
  for (i in c(1:clusters_n)) {
    a <- c()
    if (nrow(cluster_k[[i]])>5){
      for (j in c(1:nrow(cluster_k[[i]]))) {
        b <- dist(rbind(cluster_k[[i]][j,],centroids[[i]]),method = "euclidean")
        a <- c(a,b)
      }
      s[[i]] <- mean(a)
    }else{
      s[[i]] <- 0
    }
  }
  Ri <- list()
  for (i in c(1:clusters_n)){
    if (nrow(cluster_k[[i]])>5){
      r <- c()
      for (j in c(1:clusters_n)){
        if (j!=i){
          h <- (s[[i]]+s[[j]])/dist(rbind(centroids[[i]],centroids[[j]]),method = "euclidean")
          r <- c(r,h)
        }
      }
      Ri[[i]] <-  max(r,na.rm = T)
    }else{
      Ri[[i]] <- "NA"
    }
  }
  ris <- unlist(Ri);ris=ris[! ris=="NA"]
  dbi <-mean(as.numeric(ris[! ris=="NA"]))
  return(dbi)
}


calCH <- function(X,labels){
  ##X必须行为样本，列为特征
  labels2=labels
  labels=labels2[labels2 %in% names(table(labels2))[table(labels)>2]]
  labels_n <- length(unique(labels))
  samples_n <- nrow(X)
  X_mean <- apply(X,2,mean)
  ex_disp <- c()
  in_disp <- c()
  for (i in c(1:labels_n)) {
    cluster_k <- X[which(labels==i),]
    mean_k <- apply(cluster_k,2,mean)
    a1 <- nrow(cluster_k)*sum((mean_k-X_mean)^2)
    ex_disp <- c(ex_disp,a1)
    a2 <- sum((t(t(cluster_k)-mean_k))^2)
    in_disp <- c(in_disp,a2)
  }
  k1<- sum(ex_disp,na.rm=T)
  k2<- sum(in_disp,na.rm=T)
  if(k2==0)
  {
    return(1)
  }
  else
  {
    return((k1*(samples_n-labels_n))/(k2*(labels_n-1)))
  }
}




#augment_expmean=as.numeric(args[3])
samdir=paste0(dir,"/",samid,"/");#/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/151673"
#samid="151673"
#augment_expmean=1
aug.counts.file=paste0(samdir,"/augment.counts.csv")
WASfile=paste0(samdir,"/gene.Wasserstein.distance.sta")
sam.sta.file=aug.counts.file=paste0(samdir,"/sam.exp.sta")
rdsfile=paste0(samdir,"/",samid,".rds");
seurat=readRDS(rdsfile)
ground_truth <-seurat@meta.data$ground_truth
seurat <- NormalizeData(seurat,normalization.method = 'LogNormalize',scale.factor = 10000,verbose = FALSE)
#sct 与log-归一化相比，结果如何?
#  为了探究规范化方法的不同，我们研究了sctransform和log规范化结果如何与UMIs的数量相关。
# also run standard log normalization for comparison

if (!(is.null (seurat@assays$Spatial))){
  counts <- seurat@assays$Spatial@counts
  seurat <- NormalizeData(seurat, verbose = FALSE, assay = "Spatial")
}else{
  counts <- seurat@assays$RNA@counts
  seurat <- NormalizeData(seurat, verbose = FALSE, assay = "RNA")
}

seurat <-  FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 20000)
sce <- SingleCellExperiment(assays=list(counts=counts))
WASdata = read.table(WASfile,sep="\t",header = T)

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

WASdata = read.table(WASfile,sep="\t",header = T)

compareHVgeneair <- function (type,HVNUM){
  
  if (type=="WASHV"){
    #KLD散度排序方法的高变异基因
    KLD_gene <- WASdata[order(WASdata$KLD_value,decreasing = T)[1:HVNUM],1]
    #genesub=rownames(seurat)[genesum>100]
    KLD_gene  <- KLD_gene [KLD_gene  %in% rownames(seurat) ]
    # KLD_gene  <- KLD_gene [KLD_gene  %in% rownames(seurat) & KLD_gene  %in% genesub]
    seurat <- ScaleData(seurat,
                        verbose = FALSE,
                        #features = rownames(seurat)) %>% 
                        features=KLD_gene) %>%
      Seurat::RunPCA() %>%
      Seurat::RunUMAP(reduction = 'pca', dims = 1:30, verbose = FALSE)
  }
  if (type=="TOP_EXP"){
    #文件已经按照表达量高低排好序
    highexp_gene <- WASdata[1:HVNUM,1]
    #genesub=rownames(seurat)[genesum>100]
    highexp_gene  <- highexp_gene [highexp_gene  %in% rownames(seurat) ]
    # KLD_gene  <- KLD_gene [KLD_gene  %in% rownames(seurat) & KLD_gene  %in% genesub]
    seurat <- ScaleData(seurat,
                        verbose = FALSE,
                        #features = rownames(seurat)) %>% 
                        features=highexp_gene) %>%
      Seurat::RunPCA() %>%
      Seurat::RunUMAP(reduction = 'pca', dims = 1:30, verbose = FALSE)
  }
  if (type=="seurat"){
    seurat <- ScaleData(seurat,
                        verbose = FALSE,
                        #features = rownames(seurat)) %>% 
                        features=VariableFeatures(seurat)[1:HVNUM]) %>% #0.2522867和所有基因结果一样;只是取前100个基因，准确性提高到0.2797;取前1000个 0.2133098
      #features=KLD_gene) %>% #前1000个KLD散度的基因 0.2797476 #2000 0.2744991 #100 0.2377724
      Seurat::RunPCA() %>%
      Seurat::RunUMAP(reduction = 'pca', dims = 1:30, verbose = FALSE)
  }
  
  
  if (type=="BASICS"){  
    seurat <- ScaleData(seurat,
                        verbose = FALSE,
                        #features = rownames(seurat)) %>% 
                        features=BASiCS.sort.gene[1:HVNUM]) %>% #0.2522867和所有基因结果一样;只是取前100个基因，准确性提高到0.2797;取前1000个 0.2133098
      #features=KLD_gene) %>% #前1000个KLD散度的基因 0.2797476 #2000 0.2744991 #100 0.2377724
      Seurat::RunPCA() %>%
      Seurat::RunUMAP(reduction = 'pca', dims = 1:30, verbose = FALSE)
  }
  
  
  if (type=="Scran"){ 
    top.hvgs2 <- getTopHVGs(allf, n=HVNUM)
    seurat <- ScaleData(seurat,
                        verbose = FALSE,
                        #features = rownames(seurat)) %>% 
                        features=top.hvgs2) %>% #0.2522867和所有基因结果一样;只是取前100个基因，准确性提高到0.2797;取前1000个 0.2133098
      #features=KLD_gene) %>% #前1000个KLD散度的基因 0.2797476 #2000 0.2744991 #100 0.2377724
      Seurat::RunPCA() %>%
      Seurat::RunUMAP(reduction = 'pca', dims = 1:30, verbose = FALSE)
  }
  
  if (type=="Monocle"){ 
    
    monocle_hvg2 <- monocle_hvg[1:HVNUM]
    seurat <- ScaleData(seurat,
                        verbose = FALSE,
                        #features = rownames(seurat)) %>% 
                        features=monocle_hvg2) %>% #0.2522867和所有基因结果一样;只是取前100个基因，准确性提高到0.2797;取前1000个 0.2133098
      #features=KLD_gene) %>% #前1000个KLD散度的基因 0.2797476 #2000 0.2744991 #100 0.2377724
      Seurat::RunPCA() %>%
      Seurat::RunUMAP(reduction = 'pca', dims = 1:30, verbose = FALSE)
  }  
  
  if (type=="KMT"){ 
    KLD_gene <- WASdata[order(WASdata$KLD_value,decreasing = T)[1:HVNUM],1]
    #highexp_gene <- WASdata[1:HVNUM,1]
    #KMgene=unique(c(KLD_gene,highexp_gene))
    monocle_hvg2 <- monocle_hvg[1:HVNUM]
    KMgene=unique(c(KLD_gene,monocle_hvg2))
    KMgene  <- KMgene [KMgene  %in% rownames(seurat) ]
    monocle_hvg2 <- monocle_hvg[1:HVNUM]
    seurat <- ScaleData(seurat,
                        verbose = FALSE,
                        #features = rownames(seurat)) %>% 
                        features=KMgene) %>% #0.2522867和所有基因结果一样;只是取前100个基因，准确性提高到0.2797;取前1000个 0.2133098
      #features=KLD_gene) %>% #前1000个KLD散度的基因 0.2797476 #2000 0.2744991 #100 0.2377724
      Seurat::RunPCA() %>%
      Seurat::RunUMAP(reduction = 'pca', dims = 1:30, verbose = FALSE)
  }  
  
  seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)
  seq_res <- seq(0.5, 1.5, 0.1)
  seurat <- FindClusters(seurat,
                         resolution = seq_res,
                         verbose = F)
  
  
  cell_dists <- dist(seurat@reductions$pca@cell.embeddings,method = "euclidean")
  cluster_info <- seurat@meta.data[,grepl(paste0(DefaultAssay(seurat), "_snn_res"),colnames(seurat@meta.data))]
  cluster_info <- dplyr::mutate_all(cluster_info,as.character)
  cluster_info <- dplyr::mutate_all(cluster_info,as.numeric)
  
  #根据给定的墨迹聚类计算轮廓信息。 
  for(i in 1:ncol(cluster_info)){
    if (length(table(cluster_info[,i])==1)){
      cluster_info[1:2,i]=c(0,1)
    }
  }
  #根据给定的墨迹聚类计算轮廓信息。 
  silhouette_res <- apply(cluster_info, 2, function(x){
    si <- silhouette(x, cell_dists)
    mean(si[, 'sil_width'])
  })
  optm_res <- names(which.max(silhouette_res)) #选择可以另silhouette_res值最大化的分辨率
  seurat[["opt_clust"]] <- seurat[[optm_res]]
  seurat_smooth.ari <- mclust::adjustedRandIndex(ground_truth,seurat[[optm_res]][,1])
  ariresult <- rep(0, times=(length(silhouette_res)+1))
  
  x<-seurat@reductions$pca@cell.embeddings
  labels=as.numeric(seurat[[optm_res]][,1])
  SC <- mean(silhouette(labels, cell_dists)[,'sil_width']) # silhouette
  DBI <-  calDBI(x,labels) # 戴维森堡丁指数
  ch<- calCH(x,labels) #Calinski-Harabaz
  clusnum=length(unique(labels))
  cresult=list()
  cresult[[1]] <- c(samid,HVNUM,type,NA,optm_res,clusnum,SC,DBI,ch,seurat_smooth.ari);
  cresult[[2]] <- seurat[[optm_res]]
  return (cresult)
}


WAS_hv_top500 <- compareHVgeneair("WASHV",500)
WAS_hv_top1000 <- compareHVgeneair("WASHV",1000)
WAS_hv_top2000 <- compareHVgeneair("WASHV",2000)
WAS_hv_top5000 <- compareHVgeneair("WASHV",5000)
seurat_hv_top500 <- compareHVgeneair("seurat",500)
seurat_hv_top1000 <- compareHVgeneair("seurat",1000)
seurat_hv_top2000 <- compareHVgeneair("seurat",2000)
seurat_hv_top5000 <- compareHVgeneair("seurat",5000)
Scran_hv_top500 <- compareHVgeneair("Scran",500)
Scran_hv_top1000 <- compareHVgeneair("Scran",1000)
Scran_hv_top2000 <- compareHVgeneair("Scran",2000)
Scran_hv_top5000 <- compareHVgeneair("Scran",5000)
Monocle_hv_top500 <- compareHVgeneair("Monocle",500)
Monocle_hv_top1000 <- compareHVgeneair("Monocle",1000)
Monocle_hv_top2000 <- compareHVgeneair("Monocle",2000)
Monocle_hv_top5000 <- compareHVgeneair("Monocle",5000)
BASICS_hv_top500 <- compareHVgeneair("BASICS",500)
BASICS_hv_top1000 <- compareHVgeneair("BASICS",1000)
BASICS_hv_top2000 <- compareHVgeneair("BASICS",2000)
BASICS_hv_top5000 <- compareHVgeneair("BASICS",5000)

was_cr_result=rbind(WAS_hv_top500[[1]],WAS_hv_top1000[[1]],WAS_hv_top2000[[1]],WAS_hv_top5000[[1]],
                    seurat_hv_top500[[1]],seurat_hv_top1000[[1]],seurat_hv_top2000[[1]],seurat_hv_top5000[[1]],
                    Scran_hv_top500[[1]],Scran_hv_top1000[[1]],Scran_hv_top2000[[1]],Scran_hv_top5000[[1]],
                    Monocle_hv_top500[[1]],Monocle_hv_top1000[[1]],Monocle_hv_top2000[[1]],Monocle_hv_top5000[[1]],
                    BASICS_hv_top500[[1]],BASICS_hv_top1000[[1]],BASICS_hv_top2000[[1]],BASICS_hv_top5000[[1]]
                    )
cluster.all = cbind(ground_truth,WAS_hv_top500[[2]],WAS_hv_top1000[[2]],WAS_hv_top2000[[2]],WAS_hv_top5000[[2]],
                    seurat_hv_top500[[2]],seurat_hv_top1000[[2]],seurat_hv_top2000[[2]],seurat_hv_top5000[[2]],
                    Scran_hv_top500[[2]],Scran_hv_top1000[[2]],Scran_hv_top2000[[2]],Scran_hv_top5000[[2]],
                    Monocle_hv_top500[[2]],Monocle_hv_top1000[[2]],Monocle_hv_top2000[[2]],Monocle_hv_top5000[[2]],
                    BASICS_hv_top500[[2]],BASICS_hv_top1000[[2]],BASICS_hv_top2000[[2]],BASICS_hv_top5000[[2]]
                    )
hvbins=c(500,1000,2000,5000)
list.names=c(paste0("WASHV",hvbins),paste0("seurat",hvbins),paste0("Scran",hvbins),paste0("Monocle",hvbins),paste0("BASICS",hvbins))
was_cr_result=data.frame(was_cr_result)
colnames(was_cr_result)=c("samid","HV.gene.num","method","optimal","optm.res","cluster.num","SC","DBI","CH","ARI")
row.names(was_cr_result)=list.names
methods=c("WASHV","seurat","Scran","Monocle","BASICS")

for (i in 1:length(methods)){
  method=methods[i]
  tmp=was_cr_result[was_cr_result$method==method,]
  loc=(i-1)*length(hvbins)+order(tmp$SC,decreasing = T)[1]
  was_cr_result$optimal[loc]="optimal"
}

colnames(cluster.all)=c("ground_truth",list.names)
rownames(cluster.all)=colnames(seurat)

out.sta.file=paste0(samdir,"/HVgene.clustering.sta")
clustering.file=paste0(samdir,"/HVgene.clustering.result")
write.table(was_cr_result,file=out.sta.file,sep="\t",row.names=T,col.names=T,quote=F)
write.table(cluster.all,file=clustering.file,sep="\t",row.names=T,col.names=T,quote=F)



