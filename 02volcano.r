#!/bin/bash
#################################################
#  File Name:volcano_03Feb.r
#  Author: Minglu.Xie
#  Mail: minglu.xie.1023@student.uu.se,mingluxie@gmail.com
#  Created Time: Fri 03 Feb 2023 10:29:18 PM CET
#################################################

#Usage:
#/disk1/minglu/software/miniconda3/envs/r420/bin/Rscript volcano_specific.r all_peak_diff_mean.xls

library(ggplot2)
library(ggpubr)
library(ggthemes)
args = commandArgs(T)
deg.data <- read.table(args[1],header=T,sep="\t")
deg.data$logPadj <- -log10(deg.data$padj)
#deg.data$Group = "not-significant"
#deg.data$Group[which((deg.data$padj < 0.05) & (deg.data$log2FoldChange > 1))] = paste0(args[2])
#deg.data$Group[which((deg.data$padj < 0.05) & (deg.data$log2FoldChange < -1))] = paste0(args[3])
#deg.data$Group=factor(deg.data$Group,levels = c(args[2],"not-significant",args[3]))
##show the counts for each group under "sig" row.
mytable=as.data.frame(table(deg.data$sig))

#information for double check 
mytable
table(deg.data$sig)
# A data.frame: 3 × 2
# Var1	Freq
# <fct>	<int>
# mous_KO	324
# mous_NULL	481
# no diff	38371

#delete "no diff", only keep control and treated group information
mytable<- mytable[!mytable$Var1=="no diff",]

#match the color with the group information
deg.data$sig<- factor(deg.data$sig,levels =c(as.character(mytable[1,1]),as.character(mytable[2,1]),"no diff") )

p=ggscatter(deg.data, x = "log2FoldChange", y = "logPadj", color = "sig", 
	palette = c("#2f5688","#CC0000","#BBBBBB"),size=1,
	#order for the color is same as order of letters.
	#label = deg.data$Label, font.label = 8, repel = T, 
	xlab = "log2FoldChange", ylab = "-log10(Adjust P-value)",)+
theme_base()+geom_hline(yintercept = 1.3, linetype="dashed")+
geom_vline(xintercept = c(-1,1), linetype="dashed")+annotate("text", x = 5, y = max(deg.data$logPadj)*0.8, label = mytable[2,2])+ 
#annotate("text", x = 0, y = max(deg.data$logP)*0.8, label = mytable[1,2])+
annotate("text", x = -5, y = max(deg.data$logP)*0.8, label = mytable[1,2])

if(is.na(args[2])) {
    ggsave(filename = "volcano.pdf",p,width = 6,height = 4)
} else {
    ggsave(paste0(args[2],"_volcano.pdf"),p,width = 6,height = 4)
}

#save.image(paste0(args[1],'.Rdata'))
