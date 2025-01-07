args=commandArgs(T)
#data augmentation for spatial RNA sequencing data
library(Seurat)

dir=args[1]
samid=args[2]
#augment_expmean=as.numeric(args[3])
samdir=paste0(dir,"/",samid,"/");#/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/151673"
#samid="151673"
augment_expmean=1

rdsfile=paste0(samdir,"/",samid,".rds");
seurat=readRDS(rdsfile)
coordinates <- seurat@images$image@coordinates
x <- seurat@images$image@coordinates$row
y <- seurat@images$image@coordinates$col

#ground_truth=seurat@meta.data$ground_truth
saminfo <- data.frame(cbind(x,y,rownames(seurat@meta.data)))
colnames(saminfo)=c("x","y","spot") #,"ground_truth")
rownames(saminfo)=rownames(seurat@meta.data)
saminfo$x <- as.numeric(saminfo$x); saminfo$y <- as.numeric(saminfo$y)

counts <- seurat@assays$Spatial@counts

sta_counts <- function(sub){
  sub=sub[rowSums(sub)>5,]
  sparsity=mean(apply(sub,1,function(x) length(x[x==0])/length(x)))
  sta=c(samid,ncol(sub),nrow(sub),sparsity,mean(colMeans(sub)),mean(apply(sub,1,function(x) sd(x)^2)))
  sta=data.frame(sta);colnames(sta)=samid;
  row.names(sta)=c("samid","spot.num","gene.num","sparsity","exp.mean","exp.var")
  return (sta)
}
sam.sta=sta_counts(counts); exp.mean=as.numeric(sam.sta["exp.mean",])
#if the average gene expression lower than sam.sta=data.frame(sam.sta), then used the kernel methods to merge the expression neighbor spots
creat_IDW_W_matrix <- function(k,dp){
  #winsize=8;k=floor((winsize-1)/2) #这个设计往左右各扩展3个窗口
  W=data.frame(matrix(0,nrow(coordinates),nrow(coordinates)))
  row.names(W)=saminfo$spot
  colnames(W)=saminfo$spot
  for(i in 1:nrow(coordinates)){
    #  print(i)
    x <- coordinates$row[i]
    y <- coordinates$col[i]
    spot=row.names(coordinates[i,])
    xs<-x-k;xe<-x+k
    ys<-y-k;ye<-y+k
    block <-  coordinates[coordinates$row>=xs & coordinates$row<=xe & coordinates$col>=ys & coordinates$col<=ye ,] #因为旁边有空点，实际有效的邻居一般>只有4个或者更少。
    #block <- block[rownames(block) !=spot,] #k=3 最多24个邻居；k=2最多12个邻居
    #求得这些邻居的距离
    dis <- sqrt((block$row-x)^2 + (block$col-y)^2);dis[dis==0]=1
    #这里，自己和自己dis为0，这样最终权重会无穷大。我们需要处理一下，设定自己和自己的dis等于1
    w <- 1/dis^dp;  w <- w#/sum(w) #归一化成加和为1 这里不归一化成1，否则counts会变成非整数
    W[i, rownames(block)]=w ##往当前点往左右各扩展k个窗口，k个窗口内的邻居权重都是一样，行加和为1
  }
  return(W)
}

kbins=c(1,2,3)
augbins=c(5,12,24) #注意，这里的augbins是和空间转录组测序时候的孔排列有关的。10X genomics DLFPC的测序孔，每往周围扩一圈分别包括自己可以获得5个邻居；扩两圈12个，以此类推
if (exp.mean<augment_expmean){
  augr=augment_expmean/exp.mean
  k=kbins[which((round(augment_expmean/exp.mean,0)-augbins)<0)[1]]
  kersize=k*2+1
  p=0
  IDW <-  creat_IDW_W_matrix(k,p)
  counts_augment <- counts %*% as.matrix(t(IDW))
  sam.sta.aug=sta_counts(counts_augment);
  colnames(sam.sta.aug)=paste0("kersize",kersize)
  sam.sta=cbind(sam.sta,sam.sta.aug)
  write.csv(counts_augment,file=paste0(samdir,"/augment.counts.csv"))
}else{
  sam.sta=cbind(sam.sta,sam.sta)
  colnames(sam.sta)[2]="kersize1"
}
sam.sta=t(sam.sta)
write.table(sam.sta,file=paste0(samdir,"/sam.exp.sta"),sep="\t",col.names = T,row.names = T,quote=F)

