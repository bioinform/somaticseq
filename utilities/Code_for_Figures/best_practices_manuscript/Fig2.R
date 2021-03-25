setwd('~/Desktop/BestPractice_RLY_20210322/fig/')
library(reshape2)
library(ggplot2)
#######################################
# Fig2a
#######################################
dat <- read.table('~/Desktop/BestPractice_RLY_20210322/data/Fig2a.txt',header=T)
dat$x <- as.character(dat$x)
dat_long <- melt(dat)
dat_long$platform <- sapply(strsplit(as.character(dat_long$variable),"\\_"),function(x){x[1]})
dat_long$mapper <- sapply(strsplit(as.character(dat_long$variable),"\\_"),function(x){x[2]})
dat_long$x <- as.numeric(dat_long$x)

p <- ggplot(dat_long, aes(x=x, y=value, group=variable)) +
  geom_line(aes(color=mapper,linetype=platform),size=0.8)+
  scale_linetype_manual(values=c("dashed","solid"))+
  scale_color_manual(values=c('#0066B3','#FE7F2D','#999999'))+
  ylab('Cumulative overlapping SNVs')+
  xlab('Percent of SNVs in all call sets')+
  theme_classic()+
  ylim(c(0,1))+
  theme(axis.text.x = element_text(size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))

pdf('Fig2a.pdf',height=3.5,width=5)
plot(p)
dev.off()
  #######################################
# Fig2b
#######################################

dat <- read.table('~/Desktop/BestPractice_RLY_20210322/data/Fig2b.txt',header=T)
dat$x <- as.character(dat$x)
dat_long <- melt(dat)
dat_long$platform <- sapply(strsplit(as.character(dat_long$variable),"\\_"),function(x){x[1]})
dat_long$caller <- sapply(strsplit(as.character(dat_long$variable),"\\_"),function(x){x[2]})
dat_long$x <- as.numeric(dat_long$x)

p <- ggplot(dat_long, aes(x=x, y=value, group=variable)) +
  geom_line(aes(color=caller,linetype=platform),size=0.8)+
  scale_linetype_manual(values=c("dashed","solid"))+
  scale_color_manual(values=c('#E41A1C','#377EB8','#4DAF4A'))+
  ylab('Cumulative overlapping SNVs')+
  xlab('Percent of SNVs in all call sets')+
  theme_classic()+
  ylim(c(0,1))+
  theme(axis.text.x = element_text(size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))

pdf('Fig2b.pdf',height=3.5,width=5)
plot(p)
dev.off()

#######################################
# Fig2c
#######################################

dat <- read.table('~/Desktop/BestPractice_RLY_20210322/data/Fig2c.txt',header=T)
dat$x <- as.character(dat$x)
dat_long <- melt(dat)
dat_long$platform <- sapply(strsplit(as.character(dat_long$variable),"\\_"),function(x){x[1]})
dat_long$site <- sapply(strsplit(as.character(dat_long$variable),"\\_"),function(x){x[2]})
dat_long$x <- as.numeric(dat_long$x)

p <- ggplot(dat_long, aes(x=x, y=value, group=variable)) +
  geom_line(aes(color=site,linetype=platform),size=0.8)+
  scale_linetype_manual(values=c("dashed","solid"))+
  scale_color_brewer(palette="Dark2")+
  ylab('Cumulative overlapping SNVs')+
  xlab('Percent of SNVs in all call sets')+
  theme_classic()+
  ylim(c(0,1))+
  xlim(c(0,100))+
  theme(axis.text.x = element_text(size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))

pdf('Fig2c.pdf',height=3.5,width=5)
plot(p)
dev.off()

