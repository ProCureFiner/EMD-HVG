args=commandArgs(T)
#data augmentation for spatial RNA sequencing data
library(Seurat)
library("parallel")
library("foreach")
library(tidyverse)
library(cowplot)
library(ggpubr)
library("reshape2")

dir=args[1]
out_path=args[2]
project=args[4]
samfile=args[3]

#out_path="/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/integrating_12samples/"
setwd(out_path)
dir="/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/"
out_path="/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/integrating_12samples/"  
samfile="/media/disk1/SYSU/project/bin/WASHV_pipe_line/mksh_DLFPC/sample"
project="DLFPC"
setwd(out_path)
sams=read.table(samfile,sep="\t",header=F)
#dir="/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/"
#samfile="/media/disk1/SYSU/project/bin/WASHV_pipe_line/mksh_DLFPC/sample"
sams=sams[,1]

for (id in 1:length(sams)){
  samid=sams[id]
  samdir=paste0(dir,"/",samid,"/");#/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/151673"
  clustering.sta.file=paste0(samdir,"/HVgene.clustering.sta");
  data=read.table(clustering.sta.file,sep="\t",header=T)
  if (id==1){
    cdata=data;
  }else{
    cdata=rbind(cdata,data);
  }
}
write.table(cdata,file="HVgene.clustering_dataset.eachsample.eachHVgene.sta",sep="\t",col.names = T,row.names = T,quote=F)

gene.bins= c(500,1000,2000,5000)
sta.data <- function(indata,tag,colid){
  #indata=data
  gene.bins= c(500,1000,2000,5000)
  methods=unique(data$method)
  sta.mean=data.frame(matrix(NA,length(gene.bins),length(methods)))
  for (i in 1:length(gene.bins)){
    genenum=gene.bins[i]
    for(j in 1:length(methods)){
      method=methods[j]
      sta.mean[i,j]=round(mean(data[data$HV.gene.num==genenum & data$method==method,colid],na.rm=T),3)
    }
  }
  row.names(sta.mean)=gene.bins
  colnames(sta.mean)=paste0(methods,".",tag)
  return(sta.mean)
}

ins=c("ARI","SC","DBI","CH")
indexs=ins
incols=c(10,7,8,9)
sta.optimal.data <- function(data){
  methods=unique(data$method)
  sta.mean=data.frame(matrix(NA,length(methods),length(indexs)))
  for(i in 1:4){
    colid=incols[i];
    for(j in 1:length(methods)){
      method=methods[j]
      meanv=round(mean(data[data$method==method ,colid],na.rm=T),3)
      sdv=round(sd(data[data$method==method ,colid],na.rm=T),3)
      #sta.mean[j,i]=round(mean(data[data$method==method ,colid],na.rm=T),3)
      sta.mean[j,i]=paste0(meanv,"±",sdv)
    }
  }
  row.names(sta.mean)=methods
  colnames(sta.mean)=paste0(ins)
  return(sta.mean)
}

sta.optimal.data.all <- function(data){
gene.bins= c(500,1000,2000,5000)
for (i in 1:length(gene.bins)){
  genenum=gene.bins[i]
  datasub=cdata[cdata$HV.gene.num==genenum,]
  tmp=sta.optimal.data(datasub)
  tmp=cbind(methods,genenum,tmp)
  if (i==1){
    sta.all=tmp
  }else{
    sta.all=rbind(sta.all,tmp)
  }
}
return(sta.all)
}


boxplotsta <- function(data){
  data$HV.gene.num=factor(data$HV.gene.num,levels = c(500,1000,2000,5000))
  data$samid=factor(data$samid,levels = unique(data$samid))
  data$method=factor(data$method,levels = unique(data$method))
  #top gene数的影响
  p1 <- ggplot(data, aes(x=method,y=ARI,color=method)) +
    geom_boxplot() + #geom_jitter(width = 0.01,alpha = 0.2) +
    facet_wrap(~HV.gene.num,ncol=4) +
    theme(#legend.position = "none",
      legend.text = element_text(face = "bold"),legend.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45,hjust = 1,size = 10,face = "bold"),axis.text.y = element_text(size = 10,face = "bold")  ,axis.title=element_text(size=10,face="bold"))+
    ylim(quantile(data$ARI,0.01),quantile(data$ARI,0.99))+
    ggtitle(paste0(project," datasets (ARI)"))
  p2 <- ggplot(data, aes(x=method,y=SC,color=method)) +
    geom_boxplot(show.legend = FALSE) + #geom_jitter(width = 0.01,alpha = 0.2) +
    facet_wrap(~HV.gene.num,ncol=4) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 11,face = "bold"),axis.text.y = element_text(size = 11,face = "bold")
          ,axis.title=element_text(size=11,face="bold"))+
    ylim(quantile(data$SC,0.01),quantile(data$SC,0.99)) +
    ggtitle(paste0(project," datasets (SC)"))
  p3 <- ggplot(data,  aes(x=method,y=DBI,color=method)) +
    geom_boxplot(show.legend = FALSE) + #geom_jitter(width = 0.01,alpha = 0.2) +
    facet_wrap(~HV.gene.num,ncol=4) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 11,face = "bold"),axis.text.y = element_text(size = 11,face = "bold")  ,axis.title=element_text(size=11,face="bold"))+
    ylim(quantile(data$DBI,0.01),quantile(data$DBI,0.99))+
    ggtitle(paste0(project," datasets (DBI)"))
  p4 <- ggplot(data,  aes(x=method,y=CH,color=method)) +
    geom_boxplot(show.legend = FALSE) + #geom_jitter(width = 0.01,alpha = 0.2) +
    facet_wrap(~HV.gene.num,ncol=4) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 11,face = "bold"),axis.text.y = element_text(size = 11,face = "bold")  ,axis.title=element_text(size=11,face="bold"))+
    ylim(quantile(data$CH,0.01),quantile(data$CH,0.99))+
    ggtitle(paste0(project," datasets (CH)"))
  cresult=list()
  #sta.ARI=sta.data(data,"ARI",10)
  #sta.SC=sta.data(data,"SC",7)
  #sta.DBI=sta.data(data,"DBI",8)
  #sta.CH=sta.data(data,"CH",9)
  clu.sta <-sta.optimal.data.all(cdata)#cbind(project,gene.bins,sta.ARI,sta.SC,sta.DBI,sta.CH)
  cresult[[1]]=p1; cresult[[2]]=p2; cresult[[3]]=p3; cresult[[4]]=p4; cresult[[5]]=clu.sta;
  return(cresult)
}


cresult.all= boxplotsta(cdata)
cdata$method[cdata$method=="WASHV"]="EMD-HVG"
cresult.all.plot= cresult.all[[1]] + cresult.all[[2]] + cresult.all[[3]] +cresult.all[[4]] 
cresult.all.sta= cresult.all[[5]]
#print(cresult.all.plot)
#cresult.all.sta



boxplotsta.optimal <- function(data){
  data$samid=factor(data$samid,levels = unique(data$samid))
  data$method=factor(data$method,levels = unique(data$method))
  #top gene数的影响
  p1 <- ggplot(data, aes(x=method,y=ARI,color=method)) +
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1))+
    ylim(quantile(data$ARI,0.01),quantile(data$ARI,0.99))+
    ggtitle(paste0(project," datasets (ARI)"))
    
  p2 <- ggplot(data, aes(x=method,y=SC,color=method)) +
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1))+
    ylim(quantile(data$SC,0.01),quantile(data$SC,0.99))+
    ggtitle(paste0(project," datasets (SC)"))
  p3 <- ggplot(data,  aes(x=method,y=DBI,color=method)) +
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1))+
    ylim(quantile(data$DBI,0.01),quantile(data$DBI,0.99))+
    ggtitle(paste0(project," datasets (DBI)"))
  p4 <- ggplot(data,  aes(x=method,y=CH,color=method)) +
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1))+
    ylim(quantile(data$CH,0.01),quantile(data$CH,0.99))+
    ggtitle(paste0(project," datasets (CH)"))
  cresult=list()
  clu.sta <-sta.optimal.data(cdata.optimal)
  cresult[[1]]=p1; cresult[[2]]=p2; cresult[[3]]=p3; cresult[[4]]=p4; cresult[[5]]=clu.sta;
  return(cresult)
}


cdata.optimal = cdata[!is.na(cdata$optimal),]
cresult.all.optimal= boxplotsta.optimal(cdata.optimal)
cresult.all.plot.optimal= cresult.all.optimal[[1]] + cresult.all.optimal[[2]] + cresult.all.optimal[[3]] +cresult.all.optimal[[4]] 
cresult.all.sta.optimal= cresult.all.optimal[[5]]
#print(cresult.all.plot.optimal)
cresult.all.sta.optimal

pdf("HVgene.clustering_dataset_combind.pdf",width=12,height=7)
print(cresult.all.plot)
print(cresult.all.plot.optimal)
dev.off()
write.table(cresult.all.sta,file="HVgene.clustering_dataset_combind.eachHVgene.sta",sep="\t",col.names = T,row.names = T,quote=F)
write.table(cresult.all.sta.optimal,file="HVgene.clustering_dataset_combind.optimal.sta",sep="\t",col.names = T,row.names = T,quote=F)


