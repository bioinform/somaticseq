setwd('~/Desktop/BestPractice_RLY_20210322/fig/')
dat <- read.delim('~/Desktop/BestPractice_RLY_20210322/data/Fig3b.txt',header=T)
p <- ggplot(dat, aes(x=Recall, y=Precision,shape=type, color=caller)) + 
  geom_point(size=2.5)+
  scale_shape_manual(values=c(16,17))+
  theme_classic()+
  ylab('Precision')+
  xlab('Recall')+
  scale_color_brewer(palette="Set1")+
  theme(axis.text.x = element_text(size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 11,color='black'))+
  theme(axis.title.x = element_text(size = 11,color='black'))
pdf('Fig3b.pdf',height = 3,width =4.5)
plot(p)
dev.off()
