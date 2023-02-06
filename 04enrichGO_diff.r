#!/bin/bash
#################################################
#  File Name:04enrichGO_diff.r
#  Author: Minglu.Xie
#  Mail: minglu.xie.1023@student.uu.se,mingluxie@gmail.com
#  Created Time: Fri 03 Feb 2023 10:39:33 PM CET
#################################################

#Input files are from differential peak analysis results.
#Usage:
#/disk1/minglu/software/miniconda3/envs/r420/bin/Rscript enrichGO_diff.r -s mm -f knockout_up.xls -d knockout_down.xls

library("getopt")
spec<-matrix(c(
    "species", "s", "1", "character", "which species? hs, mm",
    "first", "f", "2", "character", "first differential peak file",
    "second", "d", "2", "character", "second differential peak file",
    "name", "n", "3", "chracter", "name for the output file",
    "help", "h", "0", "logical", "help"
), ncol=5, byrow=T)
opt<-getopt(spec)
# if help was asked for print a friendly message 
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}
# set some reasonable defaults for the options that are needed,
# but were not specified.
#if ( is.null(opt$species) ) {opt$species = "hs" }

if(opt$species=="hs") {
    library(ChIPseeker)
    library(ggplot2)
    library(clusterProfiler)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(EnsDb.Hsapiens.v86)
    library(org.Hs.eg.db)
    txdb=TxDb.Hsapiens.UCSC.hg38.knownGene
    aodb=org.Hs.eg.db
    aodb1="org.Hs.eg.db"
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)

  
} else if( opt$species=='mm') {
    library(ChIPseeker)
    library(clusterProfiler)
    library(ggplot2)
    library(EnsDb.Mmusculus.v79)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene
    aodb=org.Mm.eg.db
    aodb1="org.Mm.eg.db"

} else {
 print("Sorry for now, more efforts are needed to develop other species")
}

#header needs to be true, since those files have colnames.
up = read.table(opt$first,header = T)
down = read.table(opt$second,header =T)
up_name<- gsub(".xls","",opt$first)
down_name<- gsub(".xls","",opt$second)

bbb <- function(data1, name1)
{
	#Here, we need to give the colnames for the three columns in the begining postion
    # colnames(data1)[1:3]<- c("chr","start","end")
	data2=makeGRangesFromDataFrame(data1)

    #annoDb is a character!!!
	data3 = annotatePeak(data2,tssRegion = c(-3000,3000),TxDb =txdb , annoDb=aodb1)
	
	#rows of peak_annotation and orginal data
	data4 = as.data.frame(data3)

	#data5 = bitr(data4$geneId, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
	cat("rows of Annotated_peaks:",nrow(data4),"\n")
	cat("rows of Orignial_peaks:",nrow(data1),"\n")
	cat("Peak_annotated_feature:","\n")
	data6=data4[,c(1,2,3,16,6,14,15,17)]
	uniq<- paste(data6$seqnames,data6$start,data6$end,sep = "_")
	temp<- cbind(data6, uniq)
	uniq2<- paste(data1$chr,data1$start,data1$end,sep = "_")
	data1$uniq<- uniq2
	data7 <- merge(data1,temp[,4:ncol(temp)],by = "uniq")
	#enhancerRank = seq(1:nrow(data6))
	#data6$enhancer_rank = enhancerRank
	write.table(data7,file=paste0(name1,"_anno.xls"),quote=FALSE,sep="\t",row.names=FALSE)
	#save(data3,file="anno.rlt.Rdata")
	cat("Start GO enrichment: \n")
	data4 = as.data.frame(data3)
	ego_norm_BP <- enrichGO(gene = data4$geneId,
	    OrgDb = aodb,
	    ont = "BP",
	    pAdjustMethod = "BH",
	    qvalueCutoff = 0.05,
	    readable = TRUE)
	#pdf(file=paste0(name1,'.BP.annoGO.pdf'),height = 10,width = 8,paper="a4")
	pd=barplot(ego_norm_BP,showCategory=20,title=paste0("The GO_BP enrichment analysis of ", name1))
	ggsave(filename = paste0(name1,'.BP.annoGO.pdf'),plot = pd,width = 8,height = 10)
	#dev.off()
	ego_norm_BP2 = as.data.frame(ego_norm_BP)

	write.table(ego_norm_BP2,file=paste0(name1,'.BP.annoGO.xls'),sep="\t",row.names = F,quote = FALSE)
	write.table(data3@annoStat,file = paste0(name1,"_annostat.txt"),quote = F,row.names = F,sep="\t")

	return(data3)
}


up_anno = bbb(up, up_name)
down_anno = bbb(down,down_name)

all = list(up =up_anno, down = down_anno)
names(all)<- c(up_name,down_name)



# if ( !is.null(opt$help) ) {
#   cat(getopt(spec, usage=TRUE))
#   q(status=1)
# }
pdf(file='anno_diiferential_peak.pdf',width=8,height=6)
plotAnnoBar(all)
plotDistToTSS(all, title = "Distribution of peaks relative to TSS" )
dev.off()

# if(is.null(opt$name)) {
#     pdf(file='anno.pdf',width=8,height=6)
#     plotAnnoBar(all)
#     plotDistToTSS(all, title = "Distribution of peaks relative to TSS" )
#     dev.off()
    
# } else {
#     pdf(file=paste0(opt$name,'anno.pdf'),width=8,height=6)
#     plotAnnoBar(all)
#     plotDistToTSS(all, title = "Distribution of peaks relative to TSS" )
#     dev.off()
# }
