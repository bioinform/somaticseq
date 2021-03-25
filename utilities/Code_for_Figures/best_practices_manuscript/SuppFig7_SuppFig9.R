# R version 3.5.1
# set working directory
setwd('~/Desktop/SEQC_WG1_202103/data/All_three_callers/')

## define function
# split into category
split.dat <- function(dat){
  tag <- dat$High.confident.call
  # overall
  overall <- dat
  # In-truth set
  intruth <- dat[which(tag == 'HighConf' | tag == 'MedConf'),]
  # not in-truth set 
  notintruth <- dat[which(tag == 'LowConf' | tag == 'Unclassified'),]
  
  notdefined <- dat[which(tag == '0'),]
  
  return(list('overall' = overall,
              'intruth' = intruth,
              'notintruth' = notintruth,
              'notdefined' = notdefined))
}

# calulate pair-wise Jaccard index
PW.Jaccard.Index <- function(dat,WES_WGS=F){
  JI.dat <- matrix(ncol= 12,nrow = 12,NA)
  inter.dat <- matrix(ncol= 12,nrow = 12,NA)
  union.dat <- matrix(ncol= 12,nrow = 12,NA)
  if(WES_WGS==F){
    for(i in 3:ncol(dat)){
      oneSNVs <- dat[is.na(dat[,i]) == FALSE,1]
      for(j in 3:ncol(dat)){
        anotherSNVs <- dat[is.na(dat[,j]) == FALSE,1]
        
        JI.dat[i-2,j-2] <- length(intersect(oneSNVs,anotherSNVs))/length(union(oneSNVs,anotherSNVs))
        inter.dat[i-2,j-2] <- length(intersect(oneSNVs,anotherSNVs))
        union.dat[i-2,j-2] <- length(union(oneSNVs,anotherSNVs))
      }
    }
  }else{
    for(i in 3:14){
      oneSNVs <- dat[is.na(dat[,i]) == FALSE,1]
      for(j in 15:26){
        anotherSNVs <- dat[is.na(dat[,j]) == FALSE,1]
        
        JI.dat[i-2,j-14] <- length(intersect(oneSNVs,anotherSNVs))/length(union(oneSNVs,anotherSNVs))
        inter.dat[i-2,j-14] <- length(intersect(oneSNVs,anotherSNVs))
        union.dat[i-2,j-14] <- length(union(oneSNVs,anotherSNVs))
      }
    }
  }
  return(list('JI.dat'=JI.dat,'inter'=inter.dat,'union'=union.dat))
}

# calulate mean Jaccrd index
mean.dat <- function(dat){
  intra.center <- c(dat[2,3],dat[2,4],dat[3,2],dat[3,4],dat[4,2],dat[4,3],#FD
                    dat[5,6],dat[5,7],dat[6,5],dat[6,7],dat[7,5],dat[7,6],#IL
                    dat[10,11],dat[10,12],dat[11,10],dat[11,12],dat[12,10],dat[12,11] #NV
  )
  inter.center <- sum(dat) - sum(diag(dat)) - sum(intra.center)
  
  intra.center.mean <- mean(intra.center)
  inter.center.mean <- inter.center/(12*12 -12-18)
  overall.mean <- (sum(intra.center) + inter.center)/(12*12 -12)
  
  return(list('intra-center'=intra.center.mean,
              'inter-center'=inter.center.mean,
              'overall'=overall.mean))
}

# concordence column
con_col <- function(file,WES_WGS=F,part){
  dat <- read.delim(file,header=T,stringsAsFactors = F)
  dat.split <- split.dat(dat)
  result <- c(unlist(mean.dat(PW.Jaccard.Index(dat.split$overall,WES_WGS)[[part]])),
              unlist(mean.dat(PW.Jaccard.Index(dat.split$intruth,WES_WGS)[[part]])),
              unlist(mean.dat(PW.Jaccard.Index(dat.split$notintruth,WES_WGS)[[part]])),
              unlist(mean.dat(PW.Jaccard.Index(dat.split$notdefined,WES_WGS)[[part]])))
  return(result)
}

strelka.JI <- round(cbind(con_col('Somatic_AF_summary_R0_Original_WES_strelka.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_downsample_WES_strelka.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_Original_WGS_Exome_strelka.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_downsample_WGS_Exome_strelka.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_Original_WGS_strelka.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_downsample_WGS_strelka.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_Original_WGS_WES_strelka.txt',WES_WGS = T,part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_downsample_WGS_WES_strelka.txt',WES_WGS = T,part='JI.dat')),2)

strelka.inter <- round(cbind(con_col('Somatic_AF_summary_R0_Original_WES_strelka.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_downsample_WES_strelka.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_Exome_strelka.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_Exome_strelka.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_strelka.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_strelka.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_WES_strelka.txt',WES_WGS = T,part='inter'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_WES_strelka.txt',WES_WGS = T,part='inter')),0)

strelka.union <- round(cbind(con_col('Somatic_AF_summary_R0_Original_WES_strelka.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_downsample_WES_strelka.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_Exome_strelka.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_Exome_strelka.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_strelka.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_strelka.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_WES_strelka.txt',WES_WGS = T,part='union'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_WES_strelka.txt',WES_WGS = T,part='union')),0)

strelka.table1 <- matrix(paste0(strelka.JI,'(',strelka.inter,'/',strelka.union,')'),nrow = 12,ncol=8)
write.table(strelka.table1,'strelka.table1.txt',row.names=F,col.names = F,sep="\t")

mutect2.JI <- round(cbind(con_col('Somatic_AF_summary_R0_Original_WES_muTect2.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_downsample_WES_muTect2.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_Original_WGS_Exome_muTect2.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_downsample_WGS_Exome_muTect2.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_Original_WGS_muTect2.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_downsample_WGS_muTect2.txt',part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_Original_WGS_WES_muTect2.txt',WES_WGS = T,part='JI.dat'),
                          con_col('Somatic_AF_summary_R0_downsample_WGS_WES_muTect2.txt',WES_WGS = T,part='JI.dat')),2)

mutect2.inter <- round(cbind(con_col('Somatic_AF_summary_R0_Original_WES_muTect2.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_downsample_WES_muTect2.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_Exome_muTect2.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_Exome_muTect2.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_muTect2.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_muTect2.txt',part='inter'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_WES_muTect2.txt',WES_WGS = T,part='inter'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_WES_muTect2.txt',WES_WGS = T,part='inter')),0)

mutect2.union <- round(cbind(con_col('Somatic_AF_summary_R0_Original_WES_muTect2.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_downsample_WES_muTect2.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_Exome_muTect2.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_Exome_muTect2.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_muTect2.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_muTect2.txt',part='union'),
                             con_col('Somatic_AF_summary_R0_Original_WGS_WES_muTect2.txt',WES_WGS = T,part='union'),
                             con_col('Somatic_AF_summary_R0_downsample_WGS_WES_muTect2.txt',WES_WGS = T,part='union')),0)

mutect2.table1 <- matrix(paste0(mutect2.JI,'(',mutect2.inter,'/',mutect2.union,')'),nrow = 12,ncol=8)
write.table(mutect2.table1,'mutect2.table1.txt',row.names=F,col.names = F,sep="\t")

somaticSniper.JI <- round(cbind(con_col('Somatic_AF_summary_R0_Original_WES_somaticSniper.txt',part='JI.dat'),
                                con_col('Somatic_AF_summary_R0_downsample_WES_somaticSniper.txt',part='JI.dat'),
                                con_col('Somatic_AF_summary_R0_Original_WGS_Exome_somaticSniper.txt',part='JI.dat'),
                                con_col('Somatic_AF_summary_R0_downsample_WGS_Exome_somaticSniper.txt',part='JI.dat'),
                                con_col('Somatic_AF_summary_R0_Original_WGS_somaticSniper.txt',part='JI.dat'),
                                con_col('Somatic_AF_summary_R0_downsample_WGS_somaticSniper.txt',part='JI.dat'),
                                con_col('Somatic_AF_summary_R0_Original_WGS_WES_somaticSniper.txt',WES_WGS = T,part='JI.dat'),
                                con_col('Somatic_AF_summary_R0_downsample_WGS_WES_somaticSniper.txt',WES_WGS = T,part='JI.dat')),2)

somaticSniper.inter <- round(cbind(con_col('Somatic_AF_summary_R0_Original_WES_somaticSniper.txt',part='inter'),
                                   con_col('Somatic_AF_summary_R0_downsample_WES_somaticSniper.txt',part='inter'),
                                   con_col('Somatic_AF_summary_R0_Original_WGS_Exome_somaticSniper.txt',part='inter'),
                                   con_col('Somatic_AF_summary_R0_downsample_WGS_Exome_somaticSniper.txt',part='inter'),
                                   con_col('Somatic_AF_summary_R0_Original_WGS_somaticSniper.txt',part='inter'),
                                   con_col('Somatic_AF_summary_R0_downsample_WGS_somaticSniper.txt',part='inter'),
                                   con_col('Somatic_AF_summary_R0_Original_WGS_WES_somaticSniper.txt',WES_WGS = T,part='inter'),
                                   con_col('Somatic_AF_summary_R0_downsample_WGS_WES_somaticSniper.txt',WES_WGS = T,part='inter')),0)

somaticSniper.union <- round(cbind(con_col('Somatic_AF_summary_R0_Original_WES_somaticSniper.txt',part='union'),
                                   con_col('Somatic_AF_summary_R0_downsample_WES_somaticSniper.txt',part='union'),
                                   con_col('Somatic_AF_summary_R0_Original_WGS_Exome_somaticSniper.txt',part='union'),
                                   con_col('Somatic_AF_summary_R0_downsample_WGS_Exome_somaticSniper.txt',part='union'),
                                   con_col('Somatic_AF_summary_R0_Original_WGS_somaticSniper.txt',part='union'),
                                   con_col('Somatic_AF_summary_R0_downsample_WGS_somaticSniper.txt',part='union'),
                                   con_col('Somatic_AF_summary_R0_Original_WGS_WES_somaticSniper.txt',WES_WGS = T,part='union'),
                                   con_col('Somatic_AF_summary_R0_downsample_WGS_WES_somaticSniper.txt',WES_WGS = T,part='union')),0)

somaticSniper.table1 <- matrix(paste0(somaticSniper.JI,'(',somaticSniper.inter,'/',somaticSniper.union,')'),nrow = 12,ncol=8)
write.table(somaticSniper.table1,'somaticsniper.table1.txt',row.names=F,col.names = F,sep="\t")

library(reshape)
dat.reform <- function(dat){
  dat <- data.frame(dat)
  colnames(dat) <- c('WES_all_read','WES_downsample','WGS_in_exome_all_read','WGS_in_exome_downsample',
                     'WGS_whole_genome_all_read','WGS_whole_genome_downsample','WES&WGS_all_read','WES&WGS_downsample')
  dat$tag <- c('intra-center','inter-center','overall',
               'intra-center','inter-center','overall',
               'intra-center','inter-center','overall',
               'intra-center','inter-center','overall')
  dat$tag2 <- c('overall_SNVs','overall_SNVs','overall_SNVs',
                'In-truth_set','In-truth_set','In-truth_set',
                'Not-in_truth_set','Not-in_truth_set','Not-in_truth_set',
                'Not_defined','Not_defined','Not_defined')
  dat.long <- melt(dat)
  return(dat.long)
}

library(ggplot2)
plot.box <- function(dat){
  p <- ggplot(dat, aes(x=variable, y=value,fill=tag2)) +
    geom_boxplot(position=position_dodge(1))+
    facet_grid(. ~ tag2)+
    xlab('')+ylab('Jaccard Index')+
    theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1))+
    theme_bw()+theme(legend.position="none")+
    theme(axis.text.y = element_text(face="bold", size=12),
          axis.text.x = element_text(face="bold",size=12, hjust=1,vjust=1,angle=60))+
    theme(
      axis.title.x = element_text(size=14, face="bold"),
      axis.title.y = element_text(size=14, face="bold")
    )
  return(p)
}

strelka.JI.long <- dat.reform(strelka.JI)
png('strelka.JI.boxplot.png',res=200,height=1200,width=1800)
plot.box(strelka.JI.long)
dev.off()

mutect.JI.long <- dat.reform(mutect2.JI)
png('mutect.JI.boxplot.png',res=200,height=1200,width=1800)
plot.box(mutect.JI.long)
dev.off()

somaticsniper.JI.long <- dat.reform(somaticSniper.JI)
png('somaticsniper.JI.boxplot.png',res=200,height=1200,width=1800)
plot.box(somaticsniper.JI.long)
dev.off()

#######################################
# supp fig7
#######################################


strelka.JI.long$caller <- 'Strelka2'
mutect.JI.long$caller <- 'MuTect2'
somaticsniper.JI.long$caller <- 'SomticSniper'
df_long <- data.frame(rbind(strelka.JI.long,mutect.JI.long,somaticsniper.JI.long))


pdf('supp_fig7.pdf',height=4,width=8)
p <- ggplot(df_long, aes(x=tag2, y=value,fill=tag)) +
  #geom_boxplot(position=position_dodge(1))+
  geom_boxplot()+
  facet_grid(. ~ caller)+
  xlab('')+ylab('Jaccard Index')+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1))+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(axis.text.y = element_text(face="bold", size=11),
        axis.text.x = element_text(face="bold",size=11, hjust=1,vjust=1,angle=60))+
  theme(
    axis.title.x = element_text(size=12, face="bold"),
    axis.title.y = element_text(size=12, face="bold")
  )
plot(p)
dev.off()

strelka.JI.long[strelka.JI.long=="In-truth_set"]<-'Repeatable'
strelka.JI.long[strelka.JI.long=="Not_defined"]<-'Non-Repeatable'
strelka.JI.long[strelka.JI.long=="Not-in_truth_set"]<-'Gray zone'
strelka.JI.long[strelka.JI.long=="overall_SNVs"]<-'Overall'

pdf('Suppl.Fig7a_Strelka2.pdf',height=4,width=4.5)
p <- ggplot(strelka.JI.long, aes(x=tag2, y=value,fill=tag)) +
  #geom_boxplot(position=position_dodge(1))+
  geom_boxplot(lwd=0.3)+
#  facet_grid(. ~ caller)+
  xlab('')+ylab('Jaccard Index')+ggtitle("Strelka2")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  scale_fill_brewer(palette="RdBu")+
  xlim(c('Repeatable','Gray zone','Non-Repeatable','Overall'))+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,color='black'))+
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text( size=11,color='black'),
        axis.text.x = element_text(size=11, hjust=1,vjust=1,angle=60,color='black'))+
  theme(
    axis.title.x = element_text(size=12, color='black'),
    axis.title.y = element_text(size=12,color='black')
  )
plot(p)
dev.off()


somaticsniper.JI.long[somaticsniper.JI.long=="In-truth_set"]<-'Repeatable'
somaticsniper.JI.long[somaticsniper.JI.long=="Not_defined"]<-'Non-Repeatable'
somaticsniper.JI.long[somaticsniper.JI.long=="Not-in_truth_set"]<-'Gray zone'
somaticsniper.JI.long[somaticsniper.JI.long=="overall_SNVs"]<-'Overall'

pdf('Suppl.Fig7a_SomaticSniper.pdf',height=4,width=4.5)
p <- ggplot(somaticsniper.JI.long, aes(x=tag2, y=value,fill=tag)) +
  #geom_boxplot(position=position_dodge(1))+
  geom_boxplot(lwd=0.3)+
  #  facet_grid(. ~ caller)+
  xlab('')+ylab('Jaccard Index')+ggtitle("SomaticSniper")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  scale_fill_brewer(palette="RdBu")+
  xlim(c('Repeatable','Gray zone','Non-Repeatable','Overall'))+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,color='black'))+
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text( size=11,color='black'),
        axis.text.x = element_text(size=11, hjust=1,vjust=1,angle=60,color='black'))+
  theme(
    axis.title.x = element_text(size=12, color='black'),
    axis.title.y = element_text(size=12,color='black')
  )
plot(p)
dev.off()

mutect.JI.long[mutect.JI.long=="In-truth_set"]<-'Repeatable'
mutect.JI.long[mutect.JI.long=="Not_defined"]<-'Non-Repeatable'
mutect.JI.long[mutect.JI.long=="Not-in_truth_set"]<-'Gray zone'
mutect.JI.long[mutect.JI.long=="overall_SNVs"]<-'Overall'

pdf('Suppl.Fig7a_MuTect2.pdf',height=4,width=4.5)
p <- ggplot(mutect.JI.long, aes(x=tag2, y=value,fill=tag)) +
  #geom_boxplot(position=position_dodge(1))+
  geom_boxplot(lwd=0.3)+
  #  facet_grid(. ~ caller)+
  scale_fill_brewer(palette="RdBu")+
  xlab('')+ylab('Jaccard Index')+ggtitle("MuTect2")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  xlim(c('Repeatable','Gray zone','Non-Repeatable','Overall'))+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,color='black'))+
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text( size=11,color='black'),
        axis.text.x = element_text(size=11, hjust=1,vjust=1,angle=60,color='black'))+
  theme(
    axis.title.x = element_text(size=12, color='black'),
    axis.title.y = element_text(size=12,color='black')
  )
plot(p)
dev.off()

write.table(strelka.JI.long,'Suppl.Fig.7b_Strelka2.txt',col.names=T,row.names=F,sep="\t",quote=F)
write.table(mutect.JI.long,'Suppl.Fig.7a_MuTect2.txt',col.names=T,row.names=F,sep="\t",quote=F)
write.table(somaticsniper.JI.long,'Suppl.Fig.7c_SomaticSniper.txt',col.names=T,row.names=F,sep="\t",quote=F)

##################################
# mutation type comparison
# supp fig9
###################################
library(plyr)
library(reshape)
library(scales)
count_mut_type <- function(dat){
  # not in-truth set 
  tag <- dat$High.confident.call
  #x <- dat[which(tag == 'LowConf' | tag == 'Unclassified'),1]
  x <- dat[which(tag == '0'),1]
  ref <- sapply(strsplit(x,"\\|"),function(x){x[3]})
  alt <- sapply(strsplit(x,"\\|"),function(x){x[4]})
  mut <- paste(ref,alt,sep='>')
  mut <- mut[nchar(mut) == 3]
  return(count(mut))
}

mut_count <- function(file){
  dat <- read.delim(file,header=T,stringsAsFactors = F)
  return(count_mut_type(dat))
}


strelka.notdefined <- cbind(mut_count('Somatic_AF_summary_R0_Original_WES_strelka.txt'),
                            mut_count('Somatic_AF_summary_R0_downsample_WES_strelka.txt')[,2],
                            mut_count('Somatic_AF_summary_R0_Original_WGS_Exome_strelka.txt')[,2],
                            mut_count('Somatic_AF_summary_R0_downsample_WGS_Exome_strelka.txt')[,2],
                            mut_count('Somatic_AF_summary_R0_Original_WGS_strelka.txt')[,2],
                            mut_count('Somatic_AF_summary_R0_downsample_WGS_strelka.txt')[,2],
                            mut_count('Somatic_AF_summary_R0_Original_WGS_WES_strelka.txt')[,2],
                            mut_count('Somatic_AF_summary_R0_downsample_WGS_WES_strelka.txt')[,2])

mutect.notdefined <- cbind(mut_count('Somatic_AF_summary_R0_Original_WES_muTect2.txt'),
                           mut_count('Somatic_AF_summary_R0_downsample_WES_muTect2.txt')[,2],
                           mut_count('Somatic_AF_summary_R0_Original_WGS_Exome_muTect2.txt')[,2],
                           mut_count('Somatic_AF_summary_R0_downsample_WGS_Exome_muTect2.txt')[,2],
                           mut_count('Somatic_AF_summary_R0_Original_WGS_muTect2.txt')[,2],
                           mut_count('Somatic_AF_summary_R0_downsample_WGS_muTect2.txt')[,2],
                           mut_count('Somatic_AF_summary_R0_Original_WGS_WES_muTect2.txt')[,2],
                           mut_count('Somatic_AF_summary_R0_downsample_WGS_WES_muTect2.txt')[,2])


somaticsniper.notdefined <- cbind(mut_count('Somatic_AF_summary_R0_Original_WES_somaticSniper.txt'),
                                  mut_count('Somatic_AF_summary_R0_downsample_WES_somaticSniper.txt')[,2],
                                  mut_count('Somatic_AF_summary_R0_Original_WGS_Exome_somaticSniper.txt')[,2],
                                  mut_count('Somatic_AF_summary_R0_downsample_WGS_Exome_somaticSniper.txt')[,2],
                                  mut_count('Somatic_AF_summary_R0_Original_WGS_somaticSniper.txt')[,2],
                                  mut_count('Somatic_AF_summary_R0_downsample_WGS_somaticSniper.txt')[,2],
                                  mut_count('Somatic_AF_summary_R0_Original_WGS_WES_somaticSniper.txt')[,2],
                                  mut_count('Somatic_AF_summary_R0_downsample_WGS_WES_somaticSniper.txt')[,2])

colnames(strelka.notdefined) <- c('mut_type','WES_all_read','WES_downsample','WGS_in_exome_all_read','WGS_in_exome_downsample',
                   'WGS_whole_genome_all_read','WGS_whole_genome_downsample','WES&WGS_all_read','WES&WGS_downsample')
strelka.notdefined$caller <- 'Strelka2'
colnames(mutect.notdefined) <- c('mut_type','WES_all_read','WES_downsample','WGS_in_exome_all_read','WGS_in_exome_downsample',
                                  'WGS_whole_genome_all_read','WGS_whole_genome_downsample','WES&WGS_all_read','WES&WGS_downsample')
mutect.notdefined$caller <- 'MuTect2'
colnames(somaticsniper.notdefined) <- c('mut_type','WES_all_read','WES_downsample','WGS_in_exome_all_read','WGS_in_exome_downsample',
                                  'WGS_whole_genome_all_read','WGS_whole_genome_downsample','WES&WGS_all_read','WES&WGS_downsample')
somaticsniper.notdefined$caller <- 'SmaticSniper'

dat <- data.frame(rbind(strelka.notdefined,mutect.notdefined,somaticsniper.notdefined))
write.table(dat,'Suppl.Fig.9.txt',col.names=T,row.names = F,sep="\t",quote=F)
df.long <- melt(dat)
p <- ggplot(df.long, aes(fill=mut_type, y=value, x=variable)) + 
  geom_bar( stat="identity", position="fill")+
  scale_fill_brewer(palette="Paired")+
  xlab('')+ylab('Percentage')+
  facet_grid(. ~ caller)+
  theme_hc()+
  theme(axis.text.y = element_text( size=12),
        axis.text.x = element_text(size=12, hjust=1,vjust=1,angle=45))+
  theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14)
  )+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(size=12, color="black"))
pdf('Suppl.Fig.9.pdf',height=7,width=11)
plot(p)
dev.off()

 ### change format for Guan Meijian
library(reshape)
library(dplyr)
reform <- function(dat,tag){
  tmp <- PW.Jaccard.Index(dat,WES_WGS=F)[['JI.dat']]
  colnames(tmp) <- colnames(dat)[3:14]
  rownames(tmp) <- colnames(dat)[3:14]
  combine <-  data.frame(t(combn(colnames(dat)[3:14],2)))
  tmp <- melt(tmp)
  tmp3 <- inner_join(combine,tmp)
  tmp3$tag <- tag
  return(tmp3)
}


dat_reform <- function(wes_file,wgs_file,downsample=T){
  wes_dat <- read.delim(wes_file,header=T,stringsAsFactors = F)
  wgs_dat <- read.delim(wgs_file,header=T,stringsAsFactors = F)
  
  wes_dat.split <- split.dat(wes_dat)
  wes_dat.long <- rbind(reform(wes_dat.split$overall,'overall'),
                        reform(wes_dat.split$intruth,'In-truth'),
                        reform(wes_dat.split$notintruth,'Not In-truth'),
                        reform(wes_dat.split$notdefined,'Not Defined'))
  
  wgs_dat.split <- split.dat(wgs_dat)
  wgs_dat.long <- rbind(reform(wgs_dat.split$overall,'overall'),
                        reform(wgs_dat.split$intruth,'In-truth'),
                        reform(wgs_dat.split$notintruth,'Not In-truth'),
                        reform(wgs_dat.split$notdefined,'Not Defined'))
  
  caller <- gsub('.*\\.','',wes_dat.long$X1)
  if(downsample==T){
    type <- sapply(strsplit(wes_dat.long$X1,"\\."),function(x){x[3]})
  }else{
    type <- 'all-read'
  }
  SNV_subset <- wes_dat.long$tag
  
  dat <-data.frame(cbind(caller,type,SNV_subset))
  
  site1 <- sapply(strsplit(wes_dat.long$X1,"_"),function(x){x[2]})
  site2 <- sapply(strsplit(wes_dat.long$X2,"_"),function(x){x[2]})
  
  dat$pair_group <- 'inter-center'
  dat$pair_group[site1 == site2] <- 'intra-center'
  
  dat$wes <- wes_dat.long$value
  dat$wgs <- wgs_dat.long$value
  return(dat)
}

reshaped_dat <- rbind(dat_reform('Somatic_AF_summary_R0_downsample_WES_muTect2.txt','Somatic_AF_summary_R0_downsample_WGS_muTect2.txt'),
                      dat_reform('Somatic_AF_summary_R0_downsample_WES_somaticSniper.txt','Somatic_AF_summary_R0_downsample_WGS_somaticSniper.txt'),
                      dat_reform('Somatic_AF_summary_R0_downsample_WES_strelka.txt','Somatic_AF_summary_R0_downsample_WGS_strelka.txt'),
                      dat_reform('Somatic_AF_summary_R0_Original_WES_muTect2.txt','Somatic_AF_summary_R0_Original_WGS_muTect2.txt',downsample = F),
                      dat_reform('Somatic_AF_summary_R0_Original_WES_somaticSniper.txt','Somatic_AF_summary_R0_Original_WGS_somaticSniper.txt',downsample = F),
                      dat_reform('Somatic_AF_summary_R0_Original_WES_strelka.txt','Somatic_AF_summary_R0_Original_WGS_strelka.txt',downsample = F))

write.table(reshaped_dat,'inter-center_VS_intra-center.txt',col.names=T,row.names=F,sep="\t",quote=F)

# MAF
library(ggplot2)

dat <- read.delim(file,header=T,stringsAsFactors = F)
dat.split <- split.dat(dat)
list_mean <- function(x)round(rowMeans(x[,3:14],na.rm=T),2)
dat.list.mean <-lapply(dat.split,list_mean)
category <- c(rep('overall_SNVs',length(dat.list.mean$overall)),rep('In-truth',length(dat.list.mean$intruth)),
              rep('Not-in-truth',length(dat.list.mean$notintruth)),rep('Not-defined',length(dat.list.mean$notdefined)))
df <- data.frame(cbind(category,unlist(dat.list.mean)))
df$V2 <- as.numeric(as.character(df$V2))
png(plotname,height = 800,width = 1000,res=200)
p <- ggplot(df, aes(x=V2, fill=category)) +
  geom_density(alpha=0.4)+
  xlab('MAF')+ylab('Density')+
  theme(axis.text.y = element_text(face="bold", size=12),
        axis.text.x = element_text(face="bold",size=12, hjust=1,vjust=1,angle=45))+
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p
dev.off()

file = 'Somatic_AF_summary_R0_Original_WES_muTect2.txt'
plotname = 'Original_WES_muTect2.png'

file = 'Somatic_AF_summary_R0_Original_WES_somaticSniper.txt'
plotname = 'Original_WES_somaticSniper.png'

file = 'Somatic_AF_summary_R0_Original_WES_strelka.txt'
plotname = 'Original_WES_strelka.png'

file = 'Somatic_AF_summary_R0_Original_WGS_muTect2.txt'
plotname = 'Original_WGS_muTect2.png'

file = 'Somatic_AF_summary_R0_Original_WGS_somaticSniper.txt'
plotname = 'Original_WGS_somaticSniper.png'

file = 'Somatic_AF_summary_R0_Original_WGS_strelka.txt'
plotname = 'Original_WGS_strelka.png'
