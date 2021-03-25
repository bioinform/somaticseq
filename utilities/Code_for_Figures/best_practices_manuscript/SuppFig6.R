# R version 3.5.1
# set working directory
setwd('~/Desktop/SEQC_WG1_202103/')

library(VennDiagram)
library(plyr)
library(ggplot2)
library(reshape)
library(scales)
require("ggrepel")

truVCF <- read.table('high-confidence_sSNV_Indels_in_HC_regions_v1.2.V6.txt',header = T,stringsAsFactors = F)
truth <- paste(truVCF$CHROM,truVCF$POS,truVCF$REF,truVCF$ALT,sep="|")

dat <- read.delim('results_4-24-2017/FFX_Somatic_AF_summary_strelka_bwa_0.txt',header=T,stringsAsFactors = F)

ffx.bfc <- dat[is.na(dat[,'FFX_IL_24h.bfc.strelka']) == FALSE,1]
ffx.trimm <- dat[is.na(dat[,'FFX_IL_24h.trimm.strelka']) == FALSE,1]

ea1.bfc <- dat[is.na(dat[,'WES_EA_1.bfc.strelka']) == FALSE,1]
ea1.trimm <- dat[is.na(dat[,'WES_EA_1.trimm.strelka']) == FALSE,1]

fd1.bfc <- dat[is.na(dat[,'WES_FD_1.bfc.strelka']) == FALSE,1]
fd1.trimm <- dat[is.na(dat[,'WES_FD_1.trimm.strelka']) == FALSE,1]

nv1.bfc <- dat[is.na(dat[,'WES_NV_1.bfc.strelka']) == FALSE,1]
nv1.trimm <- dat[is.na(dat[,'WES_NV_1.trimm.strelka']) == FALSE,1]

df <- matrix(ncol=4,nrow=8,NA)
df[1,1] <- 'FFX'
df[1,2] <- 'BFC'
df[1,3] <- length(intersect(ffx.bfc,truth))/length(ffx.bfc)
df[1,4] <- length(intersect(ffx.bfc,truth))/length(truth)

df[2,1] <- 'FFX'
df[2,2] <- 'Trimmomatic'
df[2,3] <- length(intersect(ffx.trimm,truth))/length(ffx.trimm)
df[2,4] <- length(intersect(ffx.trimm,truth))/length(truth)

df[3,1] <- 'EA_1'
df[3,2] <- 'BFC'
df[3,3] <- length(intersect(ea1.bfc,truth))/length(ea1.bfc)
df[3,4] <- length(intersect(ea1.bfc,truth))/length(truth)

df[4,1] <- 'EA_1'
df[4,2] <- 'Trimmomatic'
df[4,3] <- length(intersect(ea1.trimm,truth))/length(ea1.trimm)
df[4,4] <- length(intersect(ea1.trimm,truth))/length(truth)

df[5,1] <- 'FD_1'
df[5,2] <- 'BFC'
df[5,3] <- length(intersect(fd1.bfc,truth))/length(fd1.bfc)
df[5,4] <- length(intersect(fd1.bfc,truth))/length(truth)

df[6,1] <- 'FD_1'
df[6,2] <- 'Trimmomatic'
df[6,3] <- length(intersect(fd1.trimm,truth))/length(fd1.trimm)
df[6,4] <- length(intersect(fd1.trimm,truth))/length(truth)

df[7,1] <- 'NV_1'
df[7,2] <- 'BFC'
df[7,3] <- length(intersect(nv1.bfc,truth))/length(nv1.bfc)
df[7,4] <- length(intersect(nv1.bfc,truth))/length(truth)

df[8,1] <- 'NV_1'
df[8,2] <- 'Trimmomatic'
df[8,3] <- length(intersect(nv1.trimm,truth))/length(nv1.trimm)
df[8,4] <- length(intersect(nv1.trimm,truth))/length(truth)

df <- data.frame(df)
colnames(df) <- c('sample','method','precision','recall')
df$precision <- as.numeric(df$precision)
df$recall <- as.numeric(df$recall)

write.table(df,'Suppl.Fig.6a.txt',sep="\t",col.names=T,row.names=F,quote=F)


p <- ggplot(df, aes(x=recall, y=precision,shape=sample, color=method)) + 
  geom_point()+geom_text_repel(aes(label=sample),
                               size = 3.5) +
  scale_shape_manual(values=c(15,16,17,8))+
  theme_classic()+
  ylab('Precision')+
  xlab('Recall')+
  theme(axis.text.x = element_text(size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 11,color='black'))+
  theme(axis.title.x = element_text(size = 11,color='black'))


pdf('Suppl.Fig.6a.pdf',height = 3,width=4.5)
plot(p)  
dev.off()


############################
# Suppl Fig6 bc
############################
dat <- read.delim('~/Desktop/BestPractice_RLY_20210322/data/Suppl.Fig.6bc.txt',header=T)
no <- dat[which(dat$type == "no-BQSR"),]
yes <- dat[which(dat$type == "BQSR"),]
p <- ggplot(dat, aes(x=Recall, y=Precision, shape=type,color=caller)) + 
  geom_point(size=2.5)+
  theme_classic()+
  ylab('Precision')+
  xlab('Recall')+
  scale_shape_manual(values=c(16,1))+
  scale_color_brewer(palette="Set1")+
  theme(axis.text.x = element_text(size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 11,color='black'))+
  theme(axis.title.x = element_text(size = 11,color='black'))
pdf('Supp.Fig.6bc.pdf',height = 3,width =4.5)
plot(p)
dev.off()

