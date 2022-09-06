

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")

library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

setwd("D:\\biowolf\\panCancer\\25.GSEA")      
gene="ABCG2"                              
gmtFile="c5.all.v7.1.symbols.gmt"          


gmt=read.gmt(gmtFile)


files=dir()
files=grep("^symbol.",files,value=T)

for(i in files){
	
	rt=read.table(i,sep="\t",header=T,check.names=F)
	CancerType=gsub("^symbol\\.|\\.txt$","",i)
	
	
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	#É¾³ýÕý³£ÑùÆ·
	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	group=gsub("2","1",group)
	data=data[,group==0]
	
	
	dataL=data[,(data[gene,]<=median(data[gene,]))]
	dataH=data[,(data[gene,]>median(data[gene,]))]
	meanL=rowMeans(dataL)
	meanH=rowMeans(dataH)
	meanL[meanL<0.00001]=0.00001
	meanH[meanH<0.00001]=0.00001
	logFC=log2(meanH/meanL)
	logFC=sort(logFC,decreasing=T)

   
    kk=GSEA(logFC,TERM2GENE=gmt, nPerm=100,pvalueCutoff = 1)
	kkTab=as.data.frame(kk)
	kkTab=kkTab[kkTab$pvalue<0.05,]
	write.table(kkTab,file=paste0("Term.",CancerType,".txt"),sep="\t",quote=F,row.names = F)
	
	
	termNum=5
	if(nrow(kkTab)>=termNum){
		gseaplot=gseaplot2(kk, row.names(kkTab)[1:termNum],base_size =8)
		pdf(file=paste0("Term.",CancerType,".pdf"),width=8,height=6)
		print(gseaplot)
		dev.off()
	}
}


