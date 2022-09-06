######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("corrplot")


#引用包
library(limma)
library(corrplot)

expFile="panGeneExp.txt"                            #表达输入文件
scoreFile="StemnessScores_RNAexp_20170127.2.tsv"    #RNAss文件
scoreType="RNAss"                                   #干细胞指数
setwd("D:\\bioinforamtion\\PANCANCER\\8.RNAss")     #设置工作目录

#读取表达文件
exp=read.table(expFile, header=T,sep="\t",row.names=1,check.names=F)
exp=exp[(exp[,"Type"]=="Tumor"),]

#读取RNAss文件
STEM=read.table(scoreFile, header=T,sep="\t",row.names=1,check.names=F)
STEM=t(STEM)

#相关性检验
outTab=data.frame()
pTab=data.frame()
#按肿瘤类型循环
for(i in levels(factor(exp[,"CancerType"]))){
    exp1=exp[(exp[,"CancerType"]==i),]
    exp1=as.matrix(exp1[,1:(ncol(exp1)-2)])
    row.names(exp1)=gsub(".$","",row.names(exp1))
    exp1=avereps(exp1)

	#样品取交集
	sameSample=intersect(row.names(STEM),row.names(exp1))
	STEM1=STEM[sameSample,]
	exp1=exp1[sameSample,]
	
    x=as.numeric(STEM1[,scoreType])
    pVector=data.frame(i)
    outVector=data.frame(i)
	#按基因循环
	genes=colnames(exp1)
	for(j in genes){
		y=as.numeric(exp1[,j])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pValue=corT$p.value
		pVector=cbind(pVector,pValue)
		outVector=cbind(outVector,cor)
	}
	pTab=rbind(pTab,pVector)
	outTab=rbind(outTab,outVector)
}
#输出相关性的表格
colNames=c("CancerType",colnames(exp)[1:(ncol(exp)-2)])
colnames(outTab)=colNames
write.table(outTab,file="RNAssCor.cor.txt",sep="\t",row.names=F,quote=F)
#输出相关性检验p值的表格
colnames(pTab)=colNames
write.table(pTab,file="RNAssCor.pval.txt",sep="\t",row.names=F,quote=F)

#RNAss相关性图形
pdf("RNAssCor.pdf",height=7,width=12)
row.names(outTab)=outTab[,1]
outTab=outTab[,-1]
corrplot(corr=as.matrix(t(outTab)),
         col=colorRampPalette(c("blue", "white", "red"))(50),
         cl.ratio = 0.6,cl.cex = 0.5,cl.offset = 1)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
