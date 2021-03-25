setwd('~/Desktop/BestPractice_RLY_20210322/fig/')
library(ggplot2)
dat <- read.delim('~/Desktop/BestPractice_RLY_20210322/data/Suppl.Fig.3.txt',header=T)
p <- ggplot(dat, aes(x=Type, y=Percent.Non.duplicated.Reads..Mapped.Trimmed., fill=Type)) +
  geom_boxplot(lwd=0.3)+
  scale_fill_manual(values=c("#999999", "#cccccc","#fda500","#fee289"))+
  # ylim(c(0,2.5))+
  xlab('')+
  ylab('Non-duplicated Reads (%)')+
  theme_classic()+
  theme(axis.text.x = element_text(size = 11,colour = 'black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,colour = 'black'))+
  theme(axis.title.y = element_text(size = 11,colour = 'black'))
pdf('Suppl.Fig.3c.pdf',height = 3,width=4.5)
plot(p) 
dev.off()

p <- ggplot(dat, aes(x=Type, y=Percent.GC, fill=Type)) +
  geom_boxplot(lwd=0.3)+
  scale_fill_manual(values=c("#999999", "#cccccc","#fda500","#fee289"))+
  # ylim(c(0,2.5))+
  xlab('')+
  ylab('GC (%)')+
  theme_classic()+
  theme(axis.text.x = element_text(size = 11,colour = 'black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,colour = 'black'))+
  theme(axis.title.y = element_text(size = 11,colour = 'black'))
pdf('Suppl.Fig.3b.pdf',height = 3,width=4.5)
plot(p) 
dev.off()

p <- ggplot(dat, aes(x=Type, y=Median.Insert.Size, fill=Type)) +
  geom_boxplot(lwd=0.3)+
  scale_fill_manual(values=c("#999999", "#cccccc","#fda500","#fee289"))+
  # ylim(c(0,2.5))+
  xlab('')+
  ylab('Median Insert Size (%)')+
  theme_classic()+
  theme(axis.text.x = element_text(size = 11,colour = 'black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,colour = 'black'))+
  theme(axis.title.y = element_text(size = 11,colour = 'black'))
pdf('Suppl.Fig.3a.pdf',height = 3,width=4.5)
plot(p) 
dev.off()

##########################
# Suppl Fig 3d
###########################

dat <- read.delim('~/Desktop/BestPractice_RLY_20210322/data/Suppl.Fig.3d.txt',header=T)
p <- ggplot(dat, aes(x=biosample, y=Percent.Reads.Mapped.On.Target, fill=biosample)) +
  geom_boxplot(lwd=0.3)+
  scale_fill_manual(values=c("#999999", "#cccccc"))+
  # ylim(c(0,2.5))+
  xlab('')+
  ylab('Median Insert Size (%)')+
  theme_classic()+
  theme(axis.text.x = element_text(size = 11,colour = 'black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,colour = 'black'))+
  theme(axis.title.y = element_text(size = 11,colour = 'black'))
pdf('Suppl.Fig.3d.pdf',height = 3,width=4.5)
plot(p) 
dev.off()
