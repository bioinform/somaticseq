setwd('~/Desktop/SEQC_WG1_202103/')
library(VennDiagram)
dat <- read.delim('results_4-24-2017/WGSvsWES_Somatic_AF_summary_3callers_bwa_0.txt',header=T,stringsAsFactors = F)
##########################
# SuppFig10a
##########################
wes <- dat[,grep('WES',colnames(dat))]
strelka_wes <- wes[,grep('strelka',colnames(wes))]
strelka_wes <- data.frame(cbind(dat[,1],strelka_wes))
wgs <- dat[,grep('WGS',colnames(dat))]
strelka_wgs <- wgs[,grep('strelka',colnames(wgs))]
strelka_wgs <- data.frame(cbind(dat[,1],strelka_wgs))
wes_unique <- data.frame()
wes_overlap <- data.frame()
wgs_unique <- data.frame()
wgs_overlap <- data.frame()
for(i in 2:13){
  wesSNVs <- strelka_wes[is.na(strelka_wes[,i]) == FALSE,1]
  wgsSNVs <- strelka_wgs[is.na(strelka_wgs[,i]) == FALSE,1]
  SNVs <- get.venn.partitions(list(wesSNVs,wgsSNVs))
  overlap <- unlist(SNVs[1,4])
  wes <- unlist(SNVs[3,4])
  wgs <- unlist(SNVs[2,4])
  wes_unique_af <- strelka_wes[,i][match(wes,strelka_wes[,1])]
  wes_overlap_af <- strelka_wes[,i][match(overlap,strelka_wes[,1])]
  wgs_unique_af <- strelka_wgs[,i][match(wgs,strelka_wgs[,1])]
  wgs_overlap_af <- strelka_wgs[,i][match(overlap,strelka_wgs[,1])]
  overlap_ref <- sapply(strsplit(overlap,"\\|"),function(x){x[3]})
  overlap_alt <- sapply(strsplit(overlap,"\\|"),function(x){x[4]})
  overlap_mut <- paste(overlap_ref,overlap_alt,sep='>')
  wes_ref <- sapply(strsplit(wes,"\\|"),function(x){x[3]})
  wes_alt <- sapply(strsplit(wes,"\\|"),function(x){x[4]})
  wes_mut <- paste(wes_ref,wes_alt,sep='>')
  wgs_ref <- sapply(strsplit(wgs,"\\|"),function(x){x[3]})
  wgs_alt <- sapply(strsplit(wgs,"\\|"),function(x){x[4]})
  wgs_mut <- paste(wgs_ref,wgs_alt,sep='>')
  site <- sapply(strsplit(colnames(strelka_wes)[i],"\\."),function(x){x[3]})
  overlap_name <- rep(site,length(overlap_mut))
  wes_name <- rep(site,length(wes_mut))
  wgs_name <- rep(site,length(wgs_mut))
  wes_unique_dat <- data.frame(cbind(wes_name,wes_unique_af,wes_mut))
  wgs_unique_dat <- data.frame(cbind(wgs_name,wgs_unique_af,wgs_mut))
  wes_overlap_dat <- data.frame(cbind(overlap_name, wes_overlap_af,overlap_mut))
  wgs_overlap_dat <- data.frame(cbind(overlap_name, wgs_overlap_af,overlap_mut))
  wes_unique <- rbind(wes_unique,wes_unique_dat)
  wes_overlap <- rbind(wes_overlap,wes_overlap_dat)
  wgs_unique <- rbind(wgs_unique,wgs_unique_dat)
  wgs_overlap <- rbind(wgs_overlap,wgs_overlap_dat)
}

wes_unique$type <- 'WES-only'
wes_overlap$type <- 'Shared'
wgs_unique$type <- 'WGS-only'
wgs_overlap$type <- 'Shared'
colnames(wes_unique) <- c('site','af','mut','type')
colnames(wes_overlap) <- c('site','af','mut','type')
colnames(wgs_unique) <- c('site','af','mut','type')
colnames(wgs_overlap) <- c('site','af','mut','type')

wes <- data.frame(rbind(wes_unique,wes_overlap))
wgs <- data.frame(rbind(wgs_unique,wgs_overlap))
wes$af <- as.numeric(wes$af)
wes$site_f = factor(wes$site, levels=c('FD_1','FD_2','FD_3','IL_1','IL_2','IL_3','NV_1','NV_2','NV_3','EA_1','LL_1','NC_1'))
wgs$af <- as.numeric(wgs$af)
wgs$site_f = factor(wgs$site, levels=c('FD_1','FD_2','FD_3','IL_1','IL_2','IL_3','NV_1','NV_2','NV_3','EA_1','LL_1','NC_1'))
p<-ggplot(wes, aes(x=af, fill=type)) +
  geom_density(alpha=0.6)+
  facet_wrap(~ site_f, nrow=1)+
  ylab('VAF Density')+
  xlab('')+
  theme_bw()+
  coord_cartesian(ylim=c(0, 3))+
  scale_fill_manual(values=c("#4A9FB6", "#999999"))+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(size=11, color="black"))+
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  
pdf('Suppl.Fig.10a.wes.pdf',height = 2,width = 15)
plot(p)
dev.off()
wgs_output <- wgs[,c(1,2,3,4)]
write.table(wgs_output,'Suppl.Fig.10a_wgs.txt',col.names=T,row.names=F,sep="\t",quote=F)
##########################
# SuppFig10b
##########################
wes_mut <- count(wes, c("site", "type",'mut'))
wgs_mut <- count(wgs, c("site", "type",'mut'))
p <- ggplot(wes_mut, aes(fill=mut, y=freq, x=site)) + 
  geom_bar( stat="identity", position="fill")+
  scale_fill_brewer(palette="Paired")+
  ylab(c('Mutation Class Frequency'))+
  xlim(c('FD_1','FD_2','FD_3','IL_1','IL_2','IL_3','NV_1','NV_2','NV_3','EA_1','LL_1','NC_1'))+
  facet_wrap(~ type)+
  xlab(c(''))+
  theme_hc()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 11,color='black'))+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(size=12, color="black"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
  

pdf('Suppl.Fig.10b.wes.pdf',height = 5,width = 7)
plot(p)
dev.off()
write.table(wes_mut,'Suppl.Fig.10b_wes.txt',col.names=T,row.names=F,sep="\t",quote=F)
write.table(wgs_mut,'Suppl.Fig.10b_wgs.txt',col.names=T,row.names=F,sep="\t",quote=F)
##########################
# SuppFig10c
##########################
library(reshape2)
truVCF <- read.table('high-confidence_sSNV_Indels_in_HC_regions_v1.2.V6.txt',header = T,stringsAsFactors = F)
truth <- paste(truVCF$CHROM,truVCF$POS,truVCF$REF,truVCF$ALT,sep="|")

af_list <- seq(0,1,0.05)
precision <- data.frame(matrix(ncol=72,nrow=20, dimnames=list(NULL, colnames(dat)[2:73])))
recall <- data.frame(matrix(ncol=72,nrow=20, dimnames=list(NULL, colnames(dat)[2:73])))
f1_score <- data.frame(matrix(ncol=72,nrow=20, dimnames=list(NULL, colnames(dat)[2:73])))
for(i in 2:73){
  for(j in 2:21){
    oneSNVs <- dat[which((dat[,i] <= af_list[j]) == TRUE),1]
    inter <- intersect(oneSNVs,truth)
    pr <- length(inter)/length(oneSNVs)
    re <- length(inter)/length(truth)
    precision[j-1,i-1] <- pr
    recall[j-1,i-1] <- re
    f1_score[j-1,i-1] <- 2*pr*re/(pr+re)
  }
}

precision$vaf <- as.character(af_list[2:21])
recall$vaf <- as.character(af_list[2:21])
f1_score$vaf <- as.character(af_list[2:21])

precision_long <- melt(precision)
recall_long <- melt(recall)
f1_long <- melt(f1_score)

precision_long <- na.omit(precision_long)
recall_long <- na.omit(recall_long)
f1_long <- na.omit(f1_long)

colnames(precision_long) <- c('vaf','sample','precision')
colnames(recall_long) <- c('vaf','sample','recall')
colnames(f1_long) <- c('vaf','sample','f1')

precision_long$caller <- sapply(strsplit(as.character(precision_long$sample),"\\."),function(x){x[1]})
precision_long$platform <- sapply(strsplit(as.character(precision_long$sample),"\\."),function(x){x[2]})

recall_long$caller <- sapply(strsplit(as.character(recall_long$sample),"\\."),function(x){x[1]})
recall_long$platform <- sapply(strsplit(as.character(recall_long$sample),"\\."),function(x){x[2]})

f1_long$caller <- sapply(strsplit(as.character(f1_long$sample),"\\."),function(x){x[1]})
f1_long$platform <- sapply(strsplit(as.character(f1_long$sample),"\\."),function(x){x[2]})


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

precision_df <- data_summary(precision_long, varname="precision", 
                    groupnames=c("vaf", "caller","platform"))
recall_df <- data_summary(recall_long, varname="recall", 
                             groupnames=c("vaf", "caller","platform"))
f1_df <- data_summary(f1_long, varname="f1", 
                             groupnames=c("vaf", "caller","platform"))

precision_df$vaf <- as.numeric(precision_df$vaf)
recall_df$vaf <- as.numeric(recall_df$vaf)
f1_df$vaf <- as.numeric(f1_df$vaf)

precision_df$group <- paste(precision_df$caller,precision_df$platform)
recall_df$group <- paste(recall_df$caller,recall_df$platform)
f1_df$group <- paste(f1_df$caller,f1_df$platform)


precision_df <- na.omit(precision_df)
recall_df <- na.omit(recall_df)
f1_df <- na.omit(f1_df)

p_plot<- ggplot(precision_df, aes(x=vaf, y=precision, group=group, color=caller)) + 
  geom_line(aes(linetype=platform))+
 # geom_point(size = 0.1)+
  geom_pointrange(aes(ymin=precision-sd, ymax=precision+sd),size=0.1)+
  scale_color_manual(values=c('#E41A1C','#377EB8','#4DAF4A'))+
  scale_linetype_manual(values=c("solid", "dotted"))+
  ylab('Precision')+
  xlab('VAF')+
  ylim(c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.2))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
  
pdf('Suppl.Fig.10c_precision.pdf',height=3.5,width=5.5)
plot(p_plot)
dev.off()

r_plot<- ggplot(recall_df, aes(x=vaf, y=recall, group=group, color=caller)) + 
  geom_line(aes(linetype=platform))+
  # geom_point(size = 0.1)+
  geom_pointrange(aes(ymin=recall-sd, ymax=recall+sd),size=0.1)+
  scale_color_manual(values=c('#E41A1C','#377EB8','#4DAF4A'))+
  scale_linetype_manual(values=c("solid", "dotted"))+
  ylab('Recall')+
  xlab('VAF')+
  ylim(c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.2))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))

pdf('Suppl.Fig.10c_recall.pdf',height=3.5,width=5.5)
plot(r_plot)
dev.off()

f_plot<- ggplot(f1_df, aes(x=vaf, y=f1, group=group, color=caller)) + 
  geom_line(aes(linetype=platform))+
  # geom_point(size = 0.1)+
  geom_pointrange(aes(ymin=f1-sd, ymax=f1+sd),size=0.1)+
  scale_color_manual(values=c('#E41A1C','#377EB8','#4DAF4A'))+
  scale_linetype_manual(values=c("solid", "dotted"))+
  ylab('F-score')+
  xlab('VAF')+
  ylim(c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.2))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))

pdf('Suppl.Fig.10c_f1.pdf',height=3.5,width=5.5)
plot(f_plot)
dev.off()


