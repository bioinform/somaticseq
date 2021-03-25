setwd('~/Desktop/SEQC_WG1_202103/')
truVCF <- read.table('high-confidence_sSNV_Indels_in_HC_regions_v1.2.V6.txt',header = T,stringsAsFactors = F)
truth <- paste(truVCF$CHROM,truVCF$POS,truVCF$REF,truVCF$ALT,sep="|")

library(reshape2)
af <- read.delim('results_4-24-2017/WGSvsWES_Somatic_AF_summary_3callers_bwa_0.txt',header=T,stringsAsFactors = F)
dp <- read.delim('results_4-24-2017/WGSvsWES_Somatic_DP_summary_3callers_bwa_0.txt',header=T,stringsAsFactors = F)

nv1_af <- af[,grep('NV_1',colnames(af))]
nv1_dp <- dp[,grep('NV_1',colnames(dp))]

nv1_af$snv <- 'FalsePositive'
nv1_af$snv[which((af$SNV %in% truth) == T)] <- 'TruePositive'

af_long <- melt(nv1_af)
af_long <- na.omit(af_long)
dp_long <- melt(nv1_dp)
dp_long <- na.omit(dp_long)

dat <- data.frame(cbind(af_long,dp_long))
dat <- dat[,c(1,2,3,5)]
colnames(dat) <- c('type','sample','af','dp')
dat$caller <- sapply(strsplit(as.character(dat$sample),"\\."),function(x){x[1]})
wes <- dat[grep('WES',dat$sample),]
wgs <- dat[grep('WGS',dat$sample),]
write.table(wes,'Suppl.Fig.11a_wes.txt',col.names=T,row.names=F,sep="\t",quote=F)
write.table(wgs,'Suppl.Fig.11b_wgs.txt',col.names=T,row.names=F,sep="\t",quote=F)

p <- ggplot(wgs, aes(x=af, y=dp, color=caller, shape=type)) +
  geom_point(alpha = 0.7) +
  scale_shape_manual(values=c(16,3))+ 
  ylab('Depth')+
  xlab('Allele Frequency')+
  theme_classic()+
  scale_color_brewer(palette="Set1")+
  theme(axis.text.x = element_text(size = 14))+ # text on X axis
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300))
pdf('Suppl.Fig.11b_wgs.pdf',height = 4.5,width = 8)
plot(p)
dev.off()

###############################
# Supp fig 11c
############################

dat <- read.delim('~/Desktop/data/Suppl.Fig.11c.txt',header=T)
dat_long <- melt(dat)
dat_long <- na.omit(dat_long)

p<- ggplot(dat_long, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot()+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  # ylim(c(0,2.5))+
  xlab('')+
  ylab('Depth')+
  theme_classic()+
  theme(axis.text.x = element_text(size = 11,colour = 'black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,colour = 'black'))+
  theme(axis.title.y = element_text(size = 11,colour = 'black'))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1200))

pdf('Suppl.Fig.11c.pdf',height = 3,width=4)
plot(p) 
dev.off()
