# R version 3.5.1
# set working directory
setwd('~/Desktop/SEQC_WG1_202103/')
dat <- read.delim('results_4-24-2017/FFX_Somatic_AF_summary_strelka_bwa_0.txt',header=T,stringsAsFactors = F)
library(VennDiagram)
library(plyr)
library(ggplot2)
library(reshape)
library(scales)

ffx.bfc <- dat[is.na(dat[,'FFX_IL_24h.bfc.strelka']) == FALSE,1]
ffx.trimm <- dat[is.na(dat[,'FFX_IL_24h.trimm.strelka']) == FALSE,1]

ea1.bfc <- dat[is.na(dat[,'WES_EA_1.bfc.strelka']) == FALSE,1]
ea1.trimm <- dat[is.na(dat[,'WES_EA_1.trimm.strelka']) == FALSE,1]

fd1.bfc <- dat[is.na(dat[,'WES_FD_1.bfc.strelka']) == FALSE,1]
fd1.trimm <- dat[is.na(dat[,'WES_FD_1.trimm.strelka']) == FALSE,1]

nv1.bfc <- dat[is.na(dat[,'WES_NV_1.bfc.strelka']) == FALSE,1]
nv1.trimm <- dat[is.na(dat[,'WES_NV_1.trimm.strelka']) == FALSE,1]

bfc <- get.venn.partitions(list(ffx.bfc,ea1.bfc,fd1.bfc,nv1.bfc))
trim <- get.venn.partitions(list(ffx.trimm,ea1.trimm,fd1.trimm,nv1.trimm))

count.mut <- function(x){
  ref <- sapply(strsplit(unlist(x),"\\|"),function(x){x[3]})
  alt <- sapply(strsplit(unlist(x),"\\|"),function(x){x[4]})
  mut <- paste(ref,alt,sep='>')
  return(count(mut))
}

df <- cbind(count.mut(bfc[1,6]),
            count.mut(trim[1,6])[,2],
            count.mut(bfc[15,6])[,2],
            count.mut(bfc[14,6])[,2],
            count.mut(bfc[12,6])[,2],
            count.mut(bfc[8,6])[,2],
            count.mut(trim[15,6])[,2],
            count.mut(trim[14,6])[,2],
            count.mut(trim[12,6])[,2],
            count.mut(trim[8,6])[,2])
colnames(df) <- c('mut','overlap.BFC','overlap.trimm','FFPE.BFC','EA_1.BFC','FD_1.BFC','NV_1.BFC',
                  'FFPE.trimm','EA_1.trimm','FD_1.trimm','NV_1.trimm')

df.long <- melt(df)

p <- ggplot(df.long, aes(fill=mut, y=value, x=variable)) + 
  geom_bar( stat="identity", position="fill")+
  scale_fill_brewer(palette="Paired")+
  ylab(c('Mutation Class Frequency'))+
  xlab(c(''))+
  theme_hc()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 11))

pdf('Fig4a.pdf',height = 6,width=6)
plot(p)  
dev.off()

write.table(df,'Fig4a.txt',col.names=T,row.names = F,sep="\t",quote = F)
