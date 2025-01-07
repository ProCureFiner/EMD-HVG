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
library("pROC")

dir="/media/disk1/SYSU/project/splatter_simulation/20000gene/"
bid="splatter_combind" #合并后的结果
bdir=paste0(dir,"/",bid,"/");#
if (! dir.exists(bdir)){
  dir.create(bdir)
}

samfile="/media/disk1/SYSU/project/bin/WASHV_pipe_line/mksh_splatter_simulation/sample"
sams=read.table(samfile,header=F)
#augment_expmean=as.numeric(args[3])
sams=as.character(t(sams))
setwd(bdir)

rmsams=c("splatter_20000_2000_2_9_0.1_0.2","splatter_20000_2000_2_9_0.1_0.3","splatter_20000_6000_6_7_0.1_0.1","splatter_20000_6000_6_7_0.1_0.3")
sams=sams[sams %in% rmsams == FALSE]

score.list=list()
precision.list=list()
cluster.list=list()
for (sid in 1:length(sams)){
  print(sid)
  samid=sams[sid]
  samdir=paste0(dir,"/",samid,"/");#/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/151673"

  hvg.score.file=paste0(dir,"/",samid,"/gene.HVG.score.sta");
  precision.file=paste0(dir,"/",samid,"/",samid,"HVG.precision.sta");
  cluster.file=paste0(dir,"/",samid,"/","HVgene.clustering.sta");
  gene.score=read.table(hvg.score.file,sep="\t",header=T)
  pre.data=read.table(precision.file,sep="\t",header=T)
  clu.data=read.table(cluster.file,sep="\t",header=T)
  score.list[[sid]]=gene.score
  precision.list[[sid]]=pre.data
  cluster.list[[sid]]=clu.data
}

#计算top2000基因的overlap，同一种方法之间计算
#reproducibility calculation
methods=c("EMD-HVG","Seurat","Scran","Monocle","BASICS")
cols=c("EMD_HVG","Seurat","Scran","Monocle","BASICS")
Jaccard_similarity <- data.frame(matrix(NA,length(methods),(length(genebins)+1)));
genebins=c(100,250,500,750,1000,1500,2000)
Jaccard_similarity [,1]=methods;

#绘制同一套数据中多个样本各项评估指标的boxplot图,这个是关于hvg基因灵敏性、特异性等指标的图示。
#将多个样本的准确性指标合并。

for(ida in 1:length(sams)){
  if (ida==1){
    precision.combind=precision.list[[ida]];
  }else{
    precision.combind=rbind(precision.combind,precision.list[[ida]]);
  }
}
#统计合并的结果，输出均值
precision.mean=precision.list[[1]];
precision.mean[,1]="All"
for (i in 1:nrow(precision.mean)){
  if (is.na(precision.mean$Top.gene[i])){
    sub=precision.combind[precision.combind$measure==precision.mean$measure[i] ,]
    precision.mean[i,4:ncol(sub)]=apply(sub[,4:ncol(sub)],2,function(x) mean(x,na.rm=T))
  }else{
    sub=precision.combind[precision.combind$measure==precision.mean$measure[i] & precision.combind$Top.gene==precision.mean$Top.gene[i],]
    precision.mean[i,4:ncol(sub)]=apply(sub[,4:ncol(sub)],2,function(x) mean(x,na.rm=T))
  }
}
precision.all=rbind(precision.combind,precision.mean)
write.table(precision.all,file="combind.HVG.precision.sta",row.names = F,col.names = T, quote=F,sep="\t")

library("reshape2")
#画图
#data.m = melt(precision.combind[precision.combind$measure !="AUC",],id.vars=c("sample","measure","Top.gene"),variable.name = "Method", value.name = "Value")
topnum=2000;
gbins=unique(precision.combind$Top.gene)[1:6]

for (i in 1:length(gbins)){
  topnum=gbins[i]  

data.m = melt(precision.combind[(precision.combind$measure =="AUC") | precision.combind$Top.gene==topnum,],id.vars=c("sample","measure","Top.gene"),variable.name = "Method", value.name = "Value")
data.m$sample <- factor(data.m$sample,levels=unique(data.m$sample))
data.m$measure=factor(data.m$measure,levels=unique(data.m$measure))
data.m$Method <- factor(data.m$Method,levels = unique(data.m$Method))
data.m$Value=as.numeric(data.m$Value)
data.m$Top.gene=factor(data.m$Top.gene,levels=as.numeric(unique(data.m$Top.gene))[order(as.numeric(unique(data.m$Top.gene)))])
#plot_A <-
 # p1=ggplot(data.m[data.m$measure=="Sensitivity",], aes(x=Method,y=Value,fill=Method)) +
#    geom_boxplot() +xlab("Top gene numbers")+
#    theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 11,face = "bold"),axis.text.y = element_text(size = 11,face = "bold")  ,axis.title=element_text(size=11,face="bold"))+
#    facet_wrap(~ Top.gene,ncol=3,scales = "free")+ylab("Sensitivity")
  p2=ggplot(data.m, aes(x=Method,y=Value,fill=Method)) +
    geom_boxplot() +xlab("Methods")+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 11,face = "bold"),axis.text.y = element_text(size = 11,face = "bold")  ,axis.title=element_text(size=11,face="bold"))+
    facet_wrap(~ measure,ncol=5,scales = "free")+ylab("Value")
  ggsave(paste0("HVG.precision.top",topnum,".combind.png"),p2,width = 14,height=3)
}



############################聚类结果合并

#绘制同一套数据中多个样本各项评估指标的boxplot图,这个是关于聚类结果指标的图示。
#将多个样本的准确性指标合并。

for(ida in 1:length(sams)){
  if (ida==1){
    cluster.combind=cluster.list[[ida]];
  }else{
    cluster.combind=rbind(cluster.combind,cluster.list[[ida]]);
  }
}
#统计合并的结果，输出均值
cluster.mean=cluster.list[[1]];
cluster.mean[,1]="All"
for (i in 1:nrow(cluster.mean)){
    sub=cluster.combind[cluster.combind$method==cluster.mean$method[i] & cluster.combind$HV.gene.num==cluster.mean$HV.gene.num[i],]
    cluster.mean[i,6:ncol(sub)]=apply(sub[,6:ncol(sub)],2,function(x) mean(x,na.rm=T))
}
cluster.all=rbind(cluster.combind,cluster.mean)
write.table(cluster.all,file="combind.cluster.precision.sta",row.names = F,col.names = T, quote=F,sep="\t")

library("reshape2")
#画图
#data.m = melt(precision.combind[precision.combind$measure !="AUC",],id.vars=c("sample","measure","Top.gene"),variable.name = "Method", value.name = "Value")
gbins=unique( unique(cluster.all$HV.gene.num))
cluster.combind=cluster.combind[,c(1,2,3,7:10)]
for (i in 1:length(gbins)){
  topnum=gbins[i]  
  data.m = melt(cluster.combind[cluster.combind$HV.gene.num==topnum,],id.vars=c("samid","HV.gene.num","method"),variable.name = "measure", value.name = "Value")
  colnames(data.m)=c("sample","Top.gene","Method","measure","Value")
  data.m$sample <- factor(data.m$sample,levels=unique(data.m$sample))
  data.m$measure=factor(data.m$measure,levels=c("ARI","SC","DBI","CH"))
  data.m$Method <- factor(data.m$Method,levels = unique(data.m$Method))
  data.m$Value=as.numeric(data.m$Value)
  data.m$Top.gene=factor(data.m$Top.gene,levels=as.numeric(unique(data.m$Top.gene))[order(as.numeric(unique(data.m$Top.gene)))])
  #plot_A <-
  # p1=ggplot(data.m[data.m$measure=="Sensitivity",], aes(x=Method,y=Value,fill=Method)) +
  #    geom_boxplot() +xlab("Top gene numbers")+
  #    theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 11,face = "bold"),axis.text.y = element_text(size = 11,face = "bold")  ,axis.title=element_text(size=11,face="bold"))+
  #    facet_wrap(~ Top.gene,ncol=3,scales = "free")+ylab("Sensitivity")
  p2=ggplot(data.m, aes(x=Method,y=Value,fill=Method)) +
    geom_boxplot() +xlab("Methods")+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 11,face = "bold"),axis.text.y = element_text(size = 11,face = "bold")  ,axis.title=element_text(size=11,face="bold"))+
    facet_wrap(~ measure,ncol=5,scales = "free")+ylab("Value")
  ggsave(paste0("cluster.precision.top",topnum,".combind.png"),p2,width = 12,height=3)
}





##################################################################################################
#20240921 modify
##################################################################################################

cells=c(2000,4000,6000,8000);groups=c(2,4,6,8);locs=c(7,7.6,8,8.3,8.7,9);scales=c(0.1,0.5,1);degrs=c(0.1,0.2,0.3);

idsort=0

for(cellbin in 1:length(cells)){
  param_clusters=groups[cellbin]
  param_celln=cells[cellbin]
  for (i in 1:length(locs)){
    param_libloc=locs[i];sum=0;print(i);
    for (j in 1:length(scales)){
      param_libscale=scales[j];
      for (k in 1:length(degrs)){
        idsort=idsort+1;
        param_degr=degrs[k]
        samid=paste0("splatter_20000_",param_celln,"_",param_clusters,"_",param_libloc,"_",param_libscale,"_",param_degr)
        print (samid)
        degfile=paste0(dir,"/",samid,"/",samid,"HVG.precision.sta");
        arifile=paste0(dir,"/",samid,"/HVgene.clustering.sta");
        
        if (!(file.exists(degfile) & file.exists(arifile) )){next}
        degdata = read.table(degfile,sep="\t",header=T)
        aridata = read.table(arifile,sep="\t",header=T)
        
        if (idsort==1){
          colnameari=row.names(aridata)
          degall=cbind(samid,param_celln,param_clusters,param_libloc,param_libscale,param_degr,degdata);
          ariall=cbind(samid,param_celln,param_clusters,param_libloc,param_libscale,param_degr,aridata);
        }else{
          tmpdeg=cbind(samid,param_celln,param_clusters,param_libloc,param_libscale,param_degr,degdata);
          tmpari=cbind(samid,param_celln,param_clusters,param_libloc,param_libscale,param_degr,aridata);
          degall=rbind(degall,tmpdeg)
          ariall=rbind(ariall,tmpari)
        }
      }
    }
  }
}

outfile.ari <-  paste0("combind.cluster.precision.v2.sta")
outfile.deg <-  paste0("combind.HVG.precision.v2.sta")
write.table(ariall,file=outfile.ari,sep="\t",row.names = F,col.names = T)
write.table(degall,file=outfile.deg,sep="\t",row.names = F,col.names = T)



#进一步统计结果
for(k in 1:4){
  if (k==1){sub=sta.by.cell};if (k==2){sub=sta.by.libloc};if (k==3){sub=sta.by.libscale};if (k==4){sub=sta.by.degr};
  cn=colnames(sub)[1]
  methods=c("EMD-HVG","Seurat","Scran","Monocle","BASICS")
  measures=c("Sensitivity","Specificity","Accuracy","Jaccard","AUC")
  parabins=unique(sub[,1])
  for (m in 1:length(measures)){
    measure=measures[m];
    sta.result=data.frame(matrix(NA,length(parabins),length(methods)+3))
    for (n in 1:length(parabins)){
      para=parabins[n];sta.result[n,1:3]=c(cn,measure,para)
      for (l in 1:length(methods)){
        method=methods[l];
        sta.result[n,l+3]=mean(sub[sub[,1]==para & sub$methods==method,m+3],na.rm=T)
      }
    }
    if (k==1 & m==1){
      sta.all=sta.result
    }else{
      sta.all=rbind(sta.all,sta.result)
    }
  }
}
colnames(sta.all)=c("para","measure","para.value",methods)

#############统计
#############

cells=c(2000,4000,6000,8000);groups=c(2,4,6,8);locs=c(7,7.6,8,8.3,8.7,9);scales=c(0.1,0.5,1);degrs=c(0.1,0.2,0.3);

cdata=ariall#degall[is.na(degall$measure)==FALSE,]
cdata$method[cdata$method=="WASHV"]="EMD-HVG"
cdata=cdata[,c(1:6,8:9,13:16)]
#colnames(cdata)=c(colnames(degall)[1:7],"HV.gene.num","WASHV","seurat","scran","Monocle","BASICS")
cdata.ori=cdata
cdata=melt(cdata.ori,value.name = "Value",variable.name = "Measure", id.vars=c(colnames(cdata[1:8])))
colnames(cdata)[7:9]=c("Top.gene","Method","measure")
cdata=cdata[is.na(cdata$samid)==FALSE,]
#cdata$Top.gene[is.na(cdata$Top.gene)]=2000 #AUC指标和topgene数量无关，原来的程序无法统计，赋值2000用于统计画图。
#incols=10
#indexs=c("DEG.recall")
measures=c("SC","DBI","CH","ARI")
sta.optimal.data.v2 <- function(data){
  data=data[is.na(data$samid)==FALSE,]
  methods=unique(data$Method)
  sta.mean=data.frame(matrix(NA,length(methods),length(measures)))
  #i=1
  #colid=incols[i];
  for(m in 1:length(methods)){
    method=methods[m]
    for (n in 1:length(measures)){
      measure=measures[n]
      values=data[data$Method==method  & data$measure==measure,]$Value
      values=values[is.na(values)==FALSE]
      if (length(values)>1){
        meanv=round(mean(values,na.rm=T),3)
        sdv=round(sd(values,na.rm=T),3)
        if (meanv=="-Inf" | meanv=="Inf" ){meanv=NA}
        sta.mean[m,n]=meanv
      }
    }
  }
  
  row.names(sta.mean)=methods
  colnames(sta.mean)=measures
  return(sta.mean)
}

sta.optimal.data.all.v2 <- function(data){
  gene.bins= c(500,1000,2000,5000)
  for (i in 1:length(gene.bins)){
    genenum=gene.bins[i]
    datasub=data[data$Top.gene==genenum,]
    tmp=sta.optimal.data.v2(datasub)
    tmp=cbind(methods,genenum,tmp)
    if (i==1){
      sta.all=tmp
    }else{
      sta.all=rbind(sta.all,tmp)
    }
  }
  return(sta.all)
}


for(cellbin in 1:length(cells)){
  param_clusters=groups[cellbin]
  param_celln=cells[cellbin]
  tmp.sta = sta.optimal.data.all.v2(cdata[cdata$param_celln==param_celln,])
  tmp.sta=cbind(param_celln,tmp.sta)
  if (cellbin==1){
    sta.by.cell=tmp.sta
  }else{
    sta.by.cell=rbind(sta.by.cell,tmp.sta)
  }
}
write.table(sta.by.cell,file=paste0("combind.cluster.precision.v2.stabycellnum.sta"),sep="\t",col.names = T,row.names = F)

for(i in 1:length(locs)){
  param_libloc=locs[i]
  tmp.sta = sta.optimal.data.all.v2(cdata[cdata$param_libloc==param_libloc,])
  tmp.sta=cbind(param_libloc,tmp.sta)
  if (i==1){
    sta.by.libloc=tmp.sta
  }else{
    sta.by.libloc=rbind(sta.by.libloc,tmp.sta)
  }
}
write.table(sta.by.libloc,file=paste0("combind.cluster.precision.v2.stabylibloc.sta"),sep="\t",col.names = T,row.names = F)

for(i in 1:length(scales)){
  param_libscale=scales[i]
  tmp.sta = sta.optimal.data.all.v2(cdata[cdata$param_libscale==param_libscale,])
  tmp.sta=cbind(param_libscale,tmp.sta)
  if (i==1){
    sta.by.libscale=tmp.sta
  }else{
    sta.by.libscale=rbind(sta.by.libscale,tmp.sta)
  }
}
write.table(sta.by.libscale,file=paste0("combind.cluster.precision.v2.stabylibscale.sta"),sep="\t",col.names = T,row.names = F)

for(i in 1:length(degrs)){
  param_degr=degrs[i]
  tmp.sta = sta.optimal.data.all.v2(cdata[cdata$param_degr==param_degr,])
  tmp.sta=cbind(param_degr,tmp.sta)
  if (i==1){
    sta.by.degr=tmp.sta
  }else{
    sta.by.degr=rbind(sta.by.degr,tmp.sta)
  }
}
write.table(sta.by.degr,file=paste0("combind.cluster.precision.v2.stabydegr.sta"),sep="\t",col.names = T,row.names = F)

#进一步统计结果 聚类结果 top 2000个基因进行聚类的结果
for(k in 1:4){
  if (k==1){sub=sta.by.cell};if (k==2){sub=sta.by.libloc};if (k==3){sub=sta.by.libscale};if (k==4){sub=sta.by.degr};
  cn=colnames(sub)[1]
  methods=c("EMD-HVG","Seurat","Scran","Monocle","BASICS")
  measures=c("SC","DBI","CH","ARI")
  parabins=unique(sub[,1])
  for (m in 1:length(measures)){
    measure=measures[m];
    sta.result=data.frame(matrix(NA,length(parabins),length(methods)+3))
    for (n in 1:length(parabins)){
      para=parabins[n];sta.result[n,1:3]=c(cn,measure,para)
      for (l in 1:length(methods)){
        method=methods[l];
        sta.result[n,l+3]=mean(sub[sub[,1]==para & sub$methods==method & sub$genenum==2000,m+3],na.rm=T)
      }
    }
    if (k==1 & m==1){
      sta.all.top2000=sta.result
    }else{
      sta.all.top2000=rbind(sta.all.top2000,sta.result)
    }
  }
}
colnames(sta.all.top2000)=c("para","measure","para.value",methods)

#write.table(sta.all,file="combind.HVG.precision.v3.sta",sep="\t",col.names = T,row.names = F,quote=F)
write.table(sta.all.top2000,file="combind.cluster.precision.top2000.v3.sta",sep="\t",col.names = T,row.names = F,quote=F)





plot_each_para <- function(pdin,paran){
  #pdin=sta.by.cell;paran="param_celln"
  pdin$methods=factor(pdin$methods,levels=unique(pdin$methods))
  #画图展示
  data=pdin
  pdata=melt(data,id.vars=c(paran,"methods","genenum"),variable.name = "Measure", value.name = "Value")
  pdata=pdata[is.na(pdata$Value)==FALSE,]
  pdata[,1]=factor(pdata[,1],levels=unique(pdata[,1])[order(unique(pdata[,1]))])
  pdata$genenum=factor(pdata$genenum,levels=c(500,1000,2000,5000))
  pdata$methods=factor(pdata$methods,levels=unique(pdata$methods))
  measures=c("SC","DBI","CH","ARI")
  for (mid in 1:length(measures)){
    measure=measures[mid]
    p1 <- ggplot(pdata[pdata$Measure==measure,], aes(x=methods,y=Value,color=methods)) +
      geom_boxplot() + #geom_jitter(width = 0.01,alpha = 0.2) +
      facet_wrap(as.formula(paste("~", paran)),ncol=4) +
      theme(axis.text.x = element_text(angle = 45,hjust = 1))+
      ggtitle(paste0(project," (",measure,")"))
  
    p12 <- ggplot(pdata[pdata$Measure==measure,], aes(x=methods,y=Value,color=methods)) +
      geom_boxplot() + #geom_jitter(width = 0.01,alpha = 0.2) +
      facet_wrap(~genenum,ncol=6) +
      theme(axis.text.x = element_text(angle = 45,hjust = 1))+
      ggtitle(paste0(project," (",measure,")"))
    ggsave(paste0("cluster.precision.",measure,".",paran,".png"),p1+p12,width=12,height=4)
  }
}
plot_each_para (sta.by.cell,"param_celln")
plot_each_para (sta.by.libloc,"param_libloc")
plot_each_para (sta.by.libscale,"param_libscale")
plot_each_para (sta.by.degr,"param_degr")



pdatain=ariall
genebins=c(500,1000,2000,5000)
for (genebinid in 1:length(genebins)){
  genenum=genebins[genebinid]
  pfile=paste0("combind.cluster.precision.top",genenum,".v3.sta.plot.pdf")
  pdf(pfile,width=12,height=8)
  par(mfrow=c(2,2))
  measures=c("ARI","SC","DBI","CH")
  for (i in 1:length(measures)){
    measure=measures[i]
    pdsub=pdatain[ pdatain$HV.gene.num==genenum,2:ncol(pdatain)]
    pdata=pdsub#melt(pdsub,id.vars=colnames(pdsub)[1:9],variable.name = "method", value.name = "Value")
    pdata$method[pdata$method=="WASHV"]="EMD-HVG"
    #colnames(pdata)=c("param","Measure","para","method","Value")
    #pdata$param=factor(pdata$param,levels=unique(pdata$param))
    paras=colnames(pdata)[c(1,3,4,5)]
    parasn=c("cell.number","lib.loc","gene.var","deg.ratio")
    pdata$param_degr=factor(pdata$param_degr,levels=unique(pdata$param_degr))
    pdata$param_celln=factor(pdata$param_celln,levels=unique(pdata$param_celln))
    pdata$param_libloc=factor(pdata$param_libloc,levels=unique(pdata$param_libloc))
    pdata$param_libscale=factor(pdata$param_libscale,levels=unique(pdata$param_libscale))
    pdata$method=factor(pdata$method,levels=unique(pdata$method))
    plotlist=list()
    for (paraid in 1:length(paras)){
      m="method";v=measure;pn=paras[paraid]
      p= ggplot(pdata, aes_string(x=pn,y=v,fill=m)) +
        geom_boxplot() + #geom_jitter(width = 0.01,alpha = 0.2) +
        labs(x=NULL)+ylab(measure)+
        theme(#legend.position = "none",
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10, face = "bold"),
          plot.title = element_text(face = "bold"),legend.text = element_text(face = "bold"),legend.title = element_text(face = "bold"))+
        ggtitle(parasn[paraid])
      plotlist[[paraid]]=p
    }
    grid.arrange(grobs = plotlist, ncol = 2)
  }
  dev.off()
}



pdatain=ariall
genebins=c(500,1000,2000,5000)
for (genebinid in 1:length(genebins)){
  genenum=genebins[genebinid]
  pfile=paste0("combind.cluster.precision.top",genenum,".v3.sta.lineplot.pdf")
  pdf(pfile,width=12,height=8)
  par(mfrow=c(2,2))
  measures=c("ARI","SC","DBI","CH")
  for (i in 1:length(measures)){
    measure=measures[i]
    pdsub=pdatain[ pdatain$HV.gene.num==genenum,2:ncol(pdatain)]
    pdata=pdsub#melt(pdsub,id.vars=colnames(pdsub)[1:9],variable.name = "method", value.name = "Value")
    pdata$method[pdata$method=="WASHV"]="EMD-HVG"
    
    #colnames(pdata)=c("param","Measure","para","method","Value")
    #pdata$param=factor(pdata$param,levels=unique(pdata$param))
    paras=colnames(pdata)[c(1,3,4,5)]
    parasn=c("cell.number","lib.loc","gene.var","deg.ratio")
    pdata$param_degr=factor(pdata$param_degr,levels=unique(pdata$param_degr))
    pdata$param_celln=factor(pdata$param_celln,levels=unique(pdata$param_celln))
    pdata$param_libloc=factor(pdata$param_libloc,levels=unique(pdata$param_libloc))
    pdata$param_libscale=factor(pdata$param_libscale,levels=unique(pdata$param_libscale))
    pdata$method=factor(pdata$method,levels=unique(pdata$method))
    pdata[,colnames(pdata)==measure]=as.numeric(pdata[,colnames(pdata)==measure])
    plotlist=list()
    for (paraid in 1:length(paras)){
      m="method";v=measure;pn=paras[paraid]
      p= ggplot(pdata, aes_string(x=pn,y=v,group=m,color=m)) +
        geom_line() + #geom_jitter(width = 0.01,alpha = 0.2) +
        labs(x=NULL)+ylab(measure)+
        theme(#legend.position = "none",
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10, face = "bold"),
          plot.title = element_text(face = "bold"),legend.text = element_text(face = "bold"),legend.title = element_text(face = "bold"))+
        ggtitle(parasn[paraid])
      plotlist[[paraid]]=p
    }
    grid.arrange(grobs = plotlist, ncol = 2)
  }
  dev.off()
}
