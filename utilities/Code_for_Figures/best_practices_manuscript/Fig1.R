setwd('~/Desktop/BestPractice_RLY_20210322/fig/')
library(ggplot2)
library(dplyr)
library(reshape2)
############################
# fig 1b
###########################
dat <- read.delim('~/Desktop/BestPractice_RLY_20210322/data/Fig1b.txt',header=T)
colnames(dat) <- c('sample','total_reads','mapped_reads','coverage')
dat$total_reads <- dat$total_reads/1000000
dat$mapped_reads <- dat$mapped_reads/1000000
dat_long <- melt(dat)

#ratio <- max(dat$total_reads)/max(dat$coverage)

p <- ggplot() +
  geom_bar(data=filter(dat_long, variable %in% c("total_reads", "mapped_reads")), aes(x = sample, y = value, fill=variable) , width = 0.7,stat ="identity", position="dodge")+
  geom_point(data=filter(dat_long, variable %in% c("coverage")),aes(x = sample, y = value*19,colour=variable)) +
  geom_line(data=filter(dat_long, variable %in% c("coverage")), aes(x = sample, y = value*19,colour=variable, group=variable)) +
  theme_classic()+
  scale_fill_manual(values=c("#177E89", "#DB3A34"))+
  scale_color_manual(values=c("#FFC857"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size=11,color='black'),legend.position="bottom")+
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 11,color='black'))+
  xlab('')+
  ylim(c(0,2500))+
  scale_y_continuous("Number of Reads (Millions)", sec.axis = sec_axis(~ . / 19, name = "Reads Coverages (X)"))

pdf('Fig1b.pdf',height = 5,width = 8)
plot(p)  
dev.off()
############################
# fig 1d
###########################

dat <- read.delim('~/Desktop/BestPractice_RLY_20210322/data/Fig1d.txt',header=T)
dat$order <- paste(dat$type,dat$hour,sep="_")
p<- ggplot(dat, aes(x=order, y=GIV,fill=order)) + 
  geom_boxplot()+
  xlim(c("FFG_1 hour","FFG_2 hour","FFG_6 hour","FFG_24 hour",'FFX_1 hour','FFX_2 hour','FFX_6 hour','FFX_24 hour'))+
  xlab("") + 
  ylab("GIV score")+
  theme_classic()+ # No grey backgroud
  theme(
    strip.text.x = element_text(
      size = 12,  face = "bold"
    ))+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 13,color='black'))+
  scale_fill_manual(values=c("#fff4d5","#fbb731","#ef7030","#f48b31","#dcdcdc","#a2a9af","#49494b","#8e8e90"))+
  theme(legend.position = "none")
  
pdf('Fig1d.pdf',height = 4,width = 4.5)
plot(p)  
dev.off()

