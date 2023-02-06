#!/bin/bash
#################################################
#  File Name:diff_peak_V2.r
#  Author: Minglu.Xie
#  Mail: minglu.xie.1023@student.uu.se,mingluxie@gmail.com
#  Created Time: Fri 03 Feb 2023 10:32:17 PM CET
#################################################

#command to run differential peak analysis
#Usage:
#ls *.rmdup.bam|awk -F "." '{if(NR<=2)print $1"\tknockout"}{if(NR>2)print $1"\tcontrol"}' > groupinfo
#/disk1/minglu/software/miniconda3/envs/r420/bin/Rscript /disk1/minglu/pipeline/Minglu_Atac/diff_peak_V2.r merged.multicov.bed groupinfo

library(limma)
library(edgeR)
library(DESeq2)
args=commandArgs(T)
data = read.table(args[1],header = F)
name = read.table(args[2],header = F)
name2=c("chr","start","end",as.character(name$V1))
colnames(data)=name2

data2 = data[,4:ncol(data)]
data3=as.data.frame(apply(data2,2,as.integer))
#rownames(data3)==NULL!!!!!careful!
mynames=paste(rep("peak_",length(data3[,1])),row.names(data2),sep="")
row.names(data3)=mynames
#colnames(data3)<- name$V1
row.names(data2)=mynames
coldata = data.frame(row.names = colnames(data3),factor(as.character(name$V2)))
colnames(coldata)=c("type")

dds <- DESeqDataSetFromMatrix(data3, coldata, design= ~ type)
dds2 <- DESeq(dds)
res= results(dds2)
###the relationship of the fold change

##tips about the comparision
resultsNames(dds2)

plus1<- strsplit(strsplit(resultsNames(dds2)[2],split = "type_")[[1]][2],split = "_vs_")[[1]][1]
minus1<-strsplit(strsplit(resultsNames(dds2)[2],split = "type_")[[1]][2],split = "_vs_")[[1]][2]
res = res[order(res$padj),]
diff_res <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
diff = as.data.frame(diff_res)
row.names(data)=mynames

peak=data[,1:3]
all=rownames(diff)
all_peak=peak[match(all,rownames(data)),]
diff_all_peak=cbind(all_peak,diff)
diff_all_peak$peak_name<- rownames(diff_all_peak)
write.table(diff_all_peak,"diff_peak_0.05.xls",sep="\t",quote=FALSE,row.names = F)

#here, we include mean value from the count matrix for each group.
#MARGIN: a vector giving the subscripts which the function will be applied over. E.g., for a matrix 1 indicates rows, 2 indicates columns, c(1, 2) indicates rows and columns. Where X has named dimnames, it can be a character vector selecting dimension names.
mymean1 = apply(data[,4:5],1,mean)
mymean2 = apply(data[,6:7],1,mean)
peak$mean1=mymean1
peak$mean2=mymean2

colnames(peak)=c("chr","start","end",as.character(name[1,2]),as.character(name[3,2]))
all_res = as.data.frame(res)
all_names = rownames(all_res)
all_peak=peak[match(all_names,rownames(data)),]
diff_all_peak2 = cbind(all_peak,all_res)
diff_all_peak2$peak_name<- rownames(diff_all_peak2)

diff_all_peak2[which(diff_all_peak2$log2FoldChange>1 & diff_all_peak2$padj<0.05),'sig']<- plus1
diff_all_peak2[which(diff_all_peak2$log2FoldChange< -1 & diff_all_peak2$padj<0.05),'sig']<- minus1
diff_all_peak2[which(abs(diff_all_peak2$log2FoldChange)<=1 | diff_all_peak2$padj>=0.05),'sig']<- 'no diff'

####remove NA value of padj
nrow(diff_all_peak2[which(is.na(diff_all_peak2$padj)),])
diff_all_peak3<- diff_all_peak2[!is.na(diff_all_peak2$padj),]
diff_up_plus1<-diff_all_peak3[which(diff_all_peak3$sig == plus1),]
diff_up_minus1<-diff_all_peak3[which(diff_all_peak3$sig == minus1),]

table(diff_all_peak3$sig)

##output will be name **up and **down
write.table(diff_all_peak3,"all_peak_diff_mean.xls",quote=FALSE,sep = "\t",row.names = F)
write.table(x = diff_up_plus1,file = paste0(plus1,'_up.xls'),quote = F,sep = "\t",row.names = F)
write.table(x = diff_up_minus1,file = paste0(plus1,'_down.xls'),quote = F,sep = "\t",row.names = F)
##output the name
write.table(data.frame(name=c(plus1,minus1)),file = "nameforanno.txt",quote = F,sep = "\t",row.names = F,col.names = F)
#save.image('my.Rdata')