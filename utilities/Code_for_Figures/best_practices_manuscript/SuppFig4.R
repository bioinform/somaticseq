setwd('~/Desktop/SEQC_WG1_202103')
library(ggplot2)
library(reshape2)
library(ggthemes)
library(scales)

######################
# supp fig4a boxplot
######################
dat <- read.delim('tonado_plot_boxplot_caller.txt',header=T,stringsAsFactors = F)
dat <- dat[1:24,]
dat <- melt(dat)
dat$platform <- gsub('.*-','',dat$Site)
dat$tag <- paste(dat$platform,dat$variable,sep='_')

p<- ggplot(dat, aes(x=tag, y=value,fill=variable)) + 
  geom_boxplot()+
#  xlim(c("MuTect2","SomaticSniper","Strelka2"))+
  xlab("") + 
  ylab("O_score")+
  theme_classic()+ # No grey backgroud
  theme(
    strip.text.x = element_text(
      size = 12,  face = "bold"
    ))+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 13,color='black'))+
  #scale_fill_brewer(palette="Set1")
  scale_fill_manual(values=c("#4DAF4A","#E41A1C","#377EB8"))

pdf('Suppl.Fig.4a.pdf',height = 4,width=5)
p
dev.off()

write.table(dat,'Suppl.Fig.4a.txt',col.names=T,row.names=F,sep="\t",quote=F)
#################################
# suppp fig4b tornado plot
#################################
abs_formatter <- function(x) {
  return(abs(x))
}

dat <- read.delim('tonado_plot_hiseq_novaseq.txt',header=T,stringsAsFactors = F)
dat$Overlap <- as.factor(dat$Overlap)
long_dat <- melt(dat)
up <- long_dat$value/2
down <- -long_dat$value/2

dat2 <- data.frame(cbind(long_dat,up,down))
dat2 <- dat2[,-3]
dat3 <- melt(dat2)
dat3$num <- as.numeric(as.character(dat3$Overlap))

colnames(dat3) <- c('overlap','platform','variable','value','num')

p <- ggplot(dat3, aes(x=num, y=value,fill=platform)) +
  geom_bar(stat="identity",position='stack',width = 0.6) + 
  facet_grid(.~platform,scales='free_x')+
  coord_flip()+scale_y_continuous(labels = abs_formatter)+
  scale_x_continuous(breaks=c(1,5,10,15,20)) +
  scale_fill_brewer(palette="Accent")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  theme(
    strip.text.x = element_text(
      size = 12,  face = "bold"
    ))+
  labs(x = "", y = "")+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,color="black", 
                                   size=12),
        axis.text.y = element_text(color="black", 
                                   size=12))

pdf('Suppl.Fig.4b.pdf',height = 3.5,width=5)
p
dev.off()
write.table(dat,'Suppl.Fig.4b.txt',col.names=T,row.names=F,sep="\t",quote=F)
