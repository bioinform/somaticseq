# R version 3.5.1
# set working directory
setwd('~/Desktop/SEQC_WG1_202103/')
library(circlize)
library(ComplexHeatmap)

dat <- read.delim('results_4-24-2017/WGSvsWES_Somatic_AF_summary_3callers_bwa_0.txt',header=T,stringsAsFactors = F)
wes <- dat[,grep('WES',colnames(dat))]
wgs <- dat[,grep('WGS',colnames(dat))]
dat.ordered <- data.frame(cbind(dat[,1],wes,wgs))


truVCF <- read.table('high-confidence_sSNV_Indels_in_HC_regions_v1.2.V6.txt',header = T,stringsAsFactors = F)
truth <- paste(truVCF$CHROM,truVCF$POS,truVCF$REF,truVCF$ALT,sep="|")


precision <- matrix(ncol=72,nrow=72,NA)
recall <- matrix(ncol=72,nrow=72,NA)

for( i in 2:73){
  oneSNVs <- dat.ordered[is.na(dat.ordered[,i]) == FALSE,1]
  for(j in 2:73){
    anotherSNVs <- dat.ordered[is.na(dat.ordered[,j]) == FALSE,1]
    overlap <- intersect(oneSNVs,anotherSNVs)
    inter <- intersect(overlap,truth)
    precision[i-1,j-1] <- length(inter)/length(overlap)
    recall[i-1,j-1] <- length(inter)/length(truth)
  }
}

colnames(precision) <- colnames(dat.ordered)[2:73]
rownames(precision) <- colnames(dat.ordered)[2:73]

colnames(recall) <- colnames(dat.ordered)[2:73]
rownames(recall) <- colnames(dat.ordered)[2:73]

caller <- sapply(strsplit(as.character(colnames(dat.ordered)[2:73]),"\\."),function(x){x[1]})
platform <- sapply(strsplit(as.character(colnames(dat.ordered)[2:73]),"\\."),function(x){x[2]})
anno_df <- data.frame(cbind(caller,platform))

ha = HeatmapAnnotation(df=anno_df,
                       col = list(platform = c("WES" = "#8c95aa", "WGS" = "#FFC857"),
                                  caller = c("muTect2"="#E41A1C","somaticSniper"="#377EB8",'strelka'='#4DAF4A')))
row_ha = rowAnnotation(df=anno_df,
                       col = list(platform = c("WES" = "#8c95aa", "WGS" = "#FFC857"),
                                  caller = c("muTect2"="#E41A1C","somaticSniper"="#377EB8",'strelka'='#4DAF4A')))

pdf('precision.pdf',height=6.2,width=8)
Heatmap(precision,top_annotation = ha,show_row_names = FALSE,left_annotation = row_ha,
        col = colorRamp2(c(0.7, 0.75, 0.85, 0.9,0.95, 1), c('dodgerblue4',"dodgerblue","white", "yellow", 'goldenrod1',"red")),
        cluster_columns = FALSE, cluster_rows = FALSE,
        rect_gp = gpar(col= "white"), 
        show_column_names = FALSE,name = 'Precision' )
dev.off()

pdf('recall.pdf',height=6.2,width=8)
Heatmap(recall,top_annotation = ha,show_row_names = FALSE,left_annotation = row_ha,
        col = colorRamp2(c(0.7, 0.75, 0.85, 0.9,0.95, 1), c('dodgerblue4',"dodgerblue","white", "yellow", 'goldenrod1',"red")),
        cluster_columns = FALSE, cluster_rows = FALSE,
        rect_gp = gpar(col= "white"), 
        show_column_names = FALSE,name = 'Recall' )
dev.off()
write.table(precision,'fig5_precision.txt',row.names=T,col.names=T,sep="\t",quote=F)
write.table(recall,'fig5_recall.txt',row.names=T,col.names=T,sep="\t",quote=F)
