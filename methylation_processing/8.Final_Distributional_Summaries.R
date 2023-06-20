setwd('F:/Research/scratch/crow_atac/2021-11/Ratios')
library(dplyr)
library(tidyr)
library(ggplot2)

### Begin with 5mC data
mfmeth <- read.table('mCG-10x_LongForm_ATAC-21Dec21.txt',header=TRUE)
head(mfmeth)
mfmeth <- mfmeth %>% mutate(Subanalysis = gsub('cgi','CpG Islands',Subanalysis),
                            Subanalysis = gsub('chromosome','Whole Chromosome',Subanalysis),
                            Subanalysis = gsub('gene','Genes',Subanalysis),
                            Sex = gsub('F','Female',Sex),
                            Sex = gsub('M','Male',Sex))

mp2 <- ggplot(mfmeth %>% subset(Subanalysis != 'Whole Chromosome'), aes(x=mCG,fill=AvZ))+
  geom_histogram(bins=50,show.legend=T,position='identity',alpha=0.6)+ggtitle('Raw 5mC Distributions')+
  facet_grid(Subanalysis+AvZ~Tissue+Sex,scales = 'free')+
  scale_fill_manual(values=c('steelblue1','tan'))+xlab('% 5mC')+ylab('Count')+
  theme_classic(base_size=15)
mp2

png('Raw_5mC_Distributions.png',units='in',res=300,height=7,width=11)
mp2
dev.off()

#and just do some sanity checks on raw distributions
d <- read.table('mCG-10x_LongForm_ATAC-21Dec21.txt',header=T)
d1 <- d %>% group_by(ID,chr,AvZ,Subanalysis,Sex,Tissue) %>% summarize(mc = mean(mCG))
#and also add start/end
co <- read.table('Gene_CGI_Chrom.coordinates',header=T)
d2 <- merge(d1,co,by=Reduce(intersect, list(names(d1),names(co)))) %>% select(-c(end,length))
d2$mc <- round(d2$mc,2)
d2$mc[d2$mc==0] <- .01
d3 <- d2 %>% group_by(ID,chr,start,Subanalysis,Tissue,AvZ) %>% summarize(fmRraw=log2(mc[Sex=='F']/mc[Sex=='M']),
                                                                         fmSraw=mc[Sex=='F']-mc[Sex=='M'])
d4 <- melt(d3,id.vars=c('ID','chr','start','Subanalysis','Tissue','AvZ'))
d4$start <- d4$start/1000000
d4 %>% subset(chr=='chrZ') %>% 
  ggplot(aes(x=start,y=value))+
  geom_point()+
  facet_grid(Subanalysis~Tissue,scales='free')+
  theme_minimal()
