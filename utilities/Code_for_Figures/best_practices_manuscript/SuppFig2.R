setwd('~/Desktop/SEQC_WG1_202103/')
library(ggplot2)
library(reshape2)
########################
# suppl fig 2c
########################
dat <- read.delim('~/Desktop/data/Suppl.Fig.2c.txt')
df_long <- melt(dat)
p <- ggplot(data=df_long, aes(x=X, y=value, group=variable,color=variable)) +
  geom_line()+
  geom_point()+
  scale_color_manual(values=c("#DB3A34", "#FFC857", "#084C61","#177E89"))+
  theme_classic()+
  xlim(c('TruSeq-PCRfree\n1000 ng','TruSeq-PCRfree\n250 ng','TruSeq-Nano\n100 ng','TruSeq-Nano\n10 ng','TruSeq-Nano\n1 ng','Nextera\n100 ng','Nextera\n10 ng','Nextera\n1 ng'))+
  xlab('')+
  ylab('Precentage')+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11,colour = 'black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,colour = 'black'))+
  theme(axis.title.y = element_text(size = 11,colour = 'black'))

pdf('Suppl.Fig2c.pdf',height=4,width=7)
plot(p)
dev.off()

########################
# suppl fig2d
######################
dat <- read.table('~/Desktop/data/Suppl.Fig.2d.txt',header = T)
p<- ggplot(dat, aes(x=Platform, y=Ratio, fill=Platform)) +
  geom_boxplot()+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
 # ylim(c(0,2.5))+
  xlab('')+
  ylab('GIV')+
  theme_classic()+
  theme(axis.text.x = element_text(size = 11,colour = 'black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,colour = 'black'))+
  theme(axis.title.y = element_text(size = 11,colour = 'black'))

pdf('Suppl.Fig.2d.pdf',height = 3,width=4)
plot(p) 
dev.off()
 



