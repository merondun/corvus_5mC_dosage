#Plot methylation value along read lengths
setwd('F:/Research/scratch/crow_atac/2021-11/QC')
library(ggplot2)
library(gridExtra)
library(ggsci)
library(viridis)

wgbs <- read.table('Master.ATAC.MBias.Nov.txt',header=F)
names(wgbs) <- c('Sample','Position','CountM','CountU','Methylated','Coverage','Read','Context')
wgbs$Experiment <- 'WGBS'
wgbs <- wgbs %>% mutate(Batch = ifelse(grepl('-NO',Sample) == T,'Novogene','SciLife'))

bias <- wgbs
bias <- rbind(wgbs,cg1,cgsub,hz1,hz2)

a <- ggplot(bias,aes(x=Position,col=Sample))+
  geom_line(aes(y=Methylated),stat="identity",show.legend=T,size=2)+ylab('Percent Methylation')+
  theme_classic(base_size=16)+facet_grid(Experiment~Read)+
  scale_color_viridis(discrete=T,option='plasma')+
  coord_cartesian(ylim=c(0,100))
a

#save
png('M_Bias-WGBS-ATAC.png',height=4,width=10,bg='transparent',units='in',res=600)
a
dev.off()

### and counts
library(viridis)
cwgbs <- read.table('Master.ATAC.Counts.Nov.txt',header=F)
names(cwgbs) <- c('Sample','Filter','Count')
cwgbs$Experiment <- 'WGBS'
mdat <- cwgbs

head(mdat)
mdat$Filter <- factor(mdat$Filter,levels=c('RAW','TsSNP','Chr','Repeatless','Repeats'))
mdat <- mdat %>% mutate(Filter = gsub('RAW','Raw CpGs',Filter),
                        Filter = gsub('TsSNP','Transition Filter',Filter),
                        Filter = gsub('Chr','Within Named Chromosomes',Filter))
                        
mdat$label <- round(mdat$Count/1000000,1)
mdat <- subset(mdat,Filter != 'Repeatless' & Filter != 'Repeats')
c <- ggplot(mdat,aes(x=Sample,y=label,fill=Filter,label=label))+
  scale_fill_manual(values=viridis(5))+
  geom_bar(stat='identity',position=position_dodge(width=1))+theme_classic(base_size=16)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))+xlab('')+ylab('')
c

png("Counts_AllExperiments.png", height=5, width=10,bg='transparent',res=600,units='in')
grid.arrange(c)
dev.off()
