args=commandArgs(T)
#data augmentation for spatial RNA sequencing data
library(Seurat)
library("parallel")
library("foreach")
library(reshape2)
library(tidyverse)
library("transport")

dir=args[1]
samid=args[2]
#augment_expmean=as.numeric(args[3])
samdir=paste0(dir,"/",samid,"/");#/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/151673"
#samid="151673"
#augment_expmean=1
WASfile=paste0(samdir,"/gene.Wasserstein.distance.sta")
if (file.exists(WASfile)){
  stop("aug.csv file.exists")
}
rdsfile=paste0(samdir,"/",samid,".rds");
seurat=readRDS(rdsfile)
if (!(is.null (seurat@assays$Spatial))){
  counts <- seurat@assays$Spatial@counts
}else{
  counts <- seurat@assays$RNA@counts
}
aug.counts.file=paste0(samdir,"/augment.counts.csv")
if (file.exists(aug.counts.file)){
  counts.aug=read.csv(aug.counts.file,row.names = );row.names(counts.aug)=counts.aug[,1];counts.aug=counts.aug[,2:ncol(counts.aug)]
  counts=as.matrix(counts.aug)
  if (is.null(seurat@assays$RNA)){
    rownames(counts)=rownames(seurat@assays$Spatial@counts)
    colnames(counts)=colnames(seurat@assays$Spatial@counts)
  }else{
    rownames(counts)=rownames(seurat@assays$RNA@counts)
    colnames(counts)=colnames(seurat@assays$RNA@counts)
  }
  write.csv(counts,aug.counts.file,row.names = T)
}
counts <- round(counts,0)
spotsum=colSums(counts)
#counts rank gene by total counts of each gene.
counts <- counts[order(rowSums(counts),decreasing = T),]
genesta <- data.frame(matrix(NA,nrow(counts),4))
colnames(genesta)=c("geneid","genesum","genemean","genevar")
genes=rownames(counts);genesta[,1]=genes
genesum=rowSums(counts);genesta[,2]=genesum
genemean=genesum/ncol(counts);genesta[,3]=genemean
genevar=apply(counts,1,function(x) sd(x)^2);genesta[,4]=genevar
####随着总样本sample size的变化，求单个基因内离散的counts的概率分布，用smoothing的方法求
######固定一次方系数为1，截距为0，求二次方系数
#系统自动根据样本数据量总counts建立合适的bins,
csum <- spotsum
#binsize=500 #设置sample size (count) bins的窗口大小，默认为500
binsize=ceiling((quantile(csum,0.75)-quantile(csum,0.25))/10/100)*200#500 #设置sample size (count) bins的窗口大小，默认为500
#maxevent=500 #设置泊松分布待考察的events最大值，默认为20
maxevent = round(mean(counts[genesta$geneid[10],])/50,0)*50 #使用表达值最大的基因的平均表达值圆整100设置maxevent。
if (maxevent==0){
  maxevent=50
}
print (paste0("maxevent:", maxevent))
mincount <- quantile(csum,0.01);maxcount <- quantile(csum,0.99)
minbin <- floor(mincount/binsize)*binsize; maxbin=ceiling(maxcount/binsize)*binsize
quantiles <- seq(minbin,maxbin, binsize)
as=quantiles[1:length(quantiles)-1]
bs=quantiles[2:length(quantiles)]
binsta <- data.frame(matrix(NA,length(as),4))

rm_hv_genes <- genesta[genesta$genevar>quantile(genesta$genevar,0.999),1]

for(i in 1:length(as)){
  a=as[i];b=bs[i]
  num=length(csum[csum>=as[i] & csum<=bs[i]])
  binsta[i,]=c(i,a,b,num)
}
colnames(binsta)=c("id","count_min","count_max","sam_num")
#去除落在bin中样本量低于20个的bin
binsta <- binsta[binsta$sam_num >=20,]

genesta_count_proportion <- function(gene){  
  #统计每个基因在不同sample size bin中检测到event=n的成概率值。
  uniqcounts=seq(0,maxevent,1)
  gsta <- data.frame(matrix(NA,length(binsta),length(uniqcounts)))
  colnames(gsta)=uniqcounts
  ys=counts[gene,]
  for(i in 1:nrow(binsta)){
    ysub=ys[csum>=binsta$count_min[i] & csum<=binsta$count_max[i]];
    for (j in 1:length(uniqcounts)){
      uc=uniqcounts[j];ysubsub=ysub[ysub==uc]
      if (length(ysub)>0 & length(ysub[ysub==uc])>0){
        gsta[i,j]=length(ysubsub)/length(ysub)
      }else{
        gsta[i,j]=0
      }
    }
  }
  gsta=data.frame(cbind(round((binsta$count_min +binsta$count_max)/2,0),gsta)); colnames(gsta)=c("samsize",uniqcounts)
  gsta.melt <- melt(gsta,id="samsize",variable.name = "event",value.name = "proportion")
  #ggplot(gsta.melt, aes(x=factor(samsize), y=proportion, colour=count,group=count)) + 
  #  geom_point(size=2) +
  #  geom_line(size=1) 
  return(gsta.melt)
  # plot(spotsum,counts[gene,])
  #  SpatialFeaturePlot(seurat,features = gene)
}

examp <-genesta_count_proportion(genes[1])
#利用每个基因表达量作为参考线，计算当前这个基因和其他基因的KLD散度
#设置用来对比的基因个数
kld_gene_num <- 16 #为了对称需要设置为偶数
library("philentropy")

genesta <- cbind(genesta,"NA")
colnames(genesta)[5] = "KLD_value"

#为了弥补高表达的基因缺少ref的问题，采用抽样合并的方法得到从5~30的基因的Refgene的分布
foraddgenes <- genesta[genesta$genemean>=0.8 & genesta$genemean<=1.2,]$geneid
genesta_count_proportion_simulate <- function(expmean){  
  samgenes <- sample(foraddgenes,round(expmean,0),replace = T)
  fakecounts <- colSums(counts[samgenes,])
  #统计每个基因在不同sample size bin中检测到event=n的成概率值。
  uniqcounts=seq(0,maxevent,1)
  gsta <- data.frame(matrix(NA,length(binsta),length(uniqcounts)))
  colnames(gsta)=uniqcounts
  ys=fakecounts
  for(i in 1:nrow(binsta)){
    ysub=ys[csum>=binsta$count_min[i] & csum<=binsta$count_max[i]];
    for (j in 1:length(uniqcounts)){
      uc=uniqcounts[j];ysubsub=ysub[ysub==uc]
      if (length(ysub)>0 & length(ysub[ysub==uc])>0){
        gsta[i,j]=length(ysubsub)/length(ysub)
      }else{
        gsta[i,j]=0
      }
    }
  }
  gsta=data.frame(cbind(round((binsta$count_min +binsta$count_max)/2,0),gsta)); colnames(gsta)=c("samsize",uniqcounts)
  gsta.melt <- melt(gsta,id="samsize",variable.name = "event",value.name = "proportion")
  #ggplot(gsta.melt, aes(x=factor(samsize), y=proportion, colour=count,group=count)) + 
  #  geom_point(size=2) +
  #  geom_line(size=1) 
  return(gsta.melt)
  # plot(spotsum,counts[gene,])
  #  SpatialFeaturePlot(seurat,features = gene)
}

result_ref <- mclapply(genes, genesta_count_proportion, mc.cores = 10)

melt.data.all <- data.frame(matrix(NA,nrow(examp),nrow(genesta)))
colnames(melt.data.all)=genesta$geneid
for (i in 1:nrow(genesta)){
  #print(i)
  gene=genesta$geneid[i]
  melt=result_ref[[i]]#genesta_count_proportion(counts,gene)
  #print (sum(melt[,3]))
  if (nrow(examp)==nrow(melt)){
    melt.data.all[,i]=melt[,3]
  }
}


was_function <- function(i){
  if (i<=kld_gene_num/2){
    aid=1;bid=kld_gene_num+1
  }else{
    aid=i-round(kld_gene_num/2); bid=i+round(kld_gene_num/2)
  }
  if (bid>=nrow(genesta)){
    aid=nrow(genesta)-kld_gene_num; bid=nrow(genesta)
  }
  
  
  genesets=genes[aid:bid]
  genea=genes[i]
  
  #去除高变异基因
  if (length(genesets [genesets %in% rm_hv_genes ])>0){
    rms = genesets [genesets %in% rm_hv_genes ] 
    for (mid in 1:length(rms)){
      if (rms[mid] != genea){
        genesets=genesets[- which (genesets == rms[mid])]
      }
    } 
  }
  
  stagenes <- data.frame(cbind (examp[,1:2],melt.data.all[,colnames(melt.data.all) %in% genesets]))
  
  if (length(genesets)<=5){
    temp=0.01;
  }else{
    
    ##20220904 当表达量高于99.5%的基因时候，用simulate的数据来当做ref
    #if ( genesta$genemean[i]>quantile(genesta$genemean,0.995)){
    #  for(j in 1:length(genesets)){
    #    gene=genesets[j]
    #    if (genea != gene){
    #      melt=genesta_count_proportion_simulate(genesta$genemean[i])
    #    }else{
    #      melt=genesta_count_proportion(gene)
    #    }
    #    if (nrow(examp)==nrow(melt)){
    #      stagenes[,j+2]=melt[,3]
    #    }
    #  }
    #}
    
    colnames(stagenes)=c("samsize","event",genesets)
    #write.table(stagenes,file="gene_count_poisson_proportion.test.sta",sep="\t",col.names = T,row.names = F)
    #转换成每个基因的概率加和为1
    psum <- colSums(data.frame(stagenes[,3:ncol(stagenes)]))
    stagenes.p=stagenes
    stagenes.p[,3:ncol(stagenes)]=t(t(stagenes[,3:ncol(stagenes)])/psum)
    #转换成长矩阵，用于画图
    stagenes.p.melt <- melt(stagenes.p,id.vars=c("samsize","event"),variable.name = "gene",value.name = "pvalue")
    id = paste0(stagenes.p.melt$samsize,"_",stagenes.p.melt$event)
    stagenes.p.melt <- cbind(id,stagenes.p.melt)
    stagenes.p.melt$id <- factor(stagenes.p.melt$id,levels=unique(stagenes.p.melt$id))
    
    x=t(stagenes.p[,3:ncol(stagenes)])
    #  x.JSD <- JSD(x);row.names(x.JSD)=row.names(x);colnames(x.JSD)=row.names(x);
    #write.table(x.JSD,file="gene.JSD.sta",sep="\t",row.names = T,col.names = T)
    x.KL <- KL(x);row.names(x.KL)=row.names(x);colnames(x.KL)=row.names(x);
    
    #观察KL均值排序（从大到小）的基因在空间中的表达分布差异
    #SpatialFeaturePlot(seurat,features = rownames(x.KL)[order(colSums(x.KL),decreasing = T)][7:12])
    #colSums(x.KL)[order(colSums(x.KL),decreasing = T)]
    
    #利用两两KL散度的聚类结果，分成两类
    #clus <- kmeans(as.matrix(x.KL),centers=2)[1]$cluster
    #获得聚类结果中样本数较多的那一组作为参考集合。测试了发现聚类效果不好，修改成为直接将KL散度均值最高的两个基因+自己从ref移除
    refgenes <- names(sort(colSums(x.KL), decreasing = T)[5:length(colSums(x.KL))])
    if (length(refgenes)<5){
      refgenes <- names(sort(colSums(x.KL), decreasing = T)[2:length(colSums(x.KL))])
    }
    p <- apply(stagenes.p[,refgenes[refgenes != genea]],1,function(x) mean(x))
    q <- stagenes.p[,genea];
    temp    <- try(wasserstein(pp(matrix(p,nrow(binsta),length(p)/nrow(binsta))),pp(matrix(q,nrow(binsta),length(p)/nrow(binsta))),p=1),silent=T)
    #计算当前基因和其他基因的KL散度均值作为得分
    if('try-error' %in% class(temp)){          # 判断当前循环的try语句中的表达式是否运行正确
      temp=wasserstein(pp(matrix(p,nrow(binsta),length(p)/nrow(binsta))),pp(matrix(q,nrow(binsta),length(p)/nrow(binsta))),p=2)                              # 此处可以对运行错误的情况进行处理应对
    }
  }
  return(temp)
}

if ( ! file.exists(WASfile)){
  ids=seq(1,nrow(genesta),1)
  result_gene<-mclapply(ids, was_function, mc.cores = 10)
  for (i in 1:nrow(genesta)){
    genesta[i,5]=result_gene[[i]]
  }
  write.table(genesta,file=  WASfile,sep="\t",row.names = T,col.names = T)
}
