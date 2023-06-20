library(tidyverse)
library(ggplot2)
library(viridis)

setwd('F:/Research/scratch/crow_atac/2022_05/bootstrapped_FM/')
d <- read.table('FM_Boot-2022-05-17.txt',header=T)
head(d)

d2 <- read.table('mCG-10x_LongForm_ATAC-21Dec21.txt',header=T)
d2 <- d2[,c('ID','length','AvZ','chr','Subanalysis','Tissue')]

#merge
d3 <- merge(d %>% select(-row),d2,by=c('ID','Tissue','Subanalysis')) %>% unique()

#add coordinates
d4 <- read.table('Coordinates.bed',header=F)
names(d4) <- c('chr','start','end','ID')

d5 <- merge(d3,d4,by=c('chr','ID'))
d5 <- d5 %>% select(chr,start,end,ID,length,AvZ,Tissue,Subanalysis,median,lowerCI,upperCI)

#summary plots
d5 %>% subset(Subanalysis =='gene') %>% ggplot(aes(x=chr,y=median,fill=AvZ))+
  geom_boxplot()+
  facet_grid(Tissue~.,scales='free')+
  theme_minimal()

write.table(d5,file='5mC_FM-Bootstrapped__2022-05-17.txt',quote=F,sep='\t',row.names=F)
