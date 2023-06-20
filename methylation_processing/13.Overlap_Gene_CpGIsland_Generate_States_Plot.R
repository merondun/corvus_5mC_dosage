setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/ATACseq_Methylation_Analyses/2021_NOV/Clip_9_Q20/2022_12_19')
.libPaths('~/miniconda3/envs/r/lib/R/library')
library(dplyr)
library(tidyverse)
library(reshape2)
library(matrixStats)
library(viridis)
library(ggpubr)
library(rstatix)

mc = read.table('FM_Boot-2022-05-17.txt',header=TRUE)
mc = mc %>% dplyr::rename(mc = median,mc_lo = lowerCI,mc_hi = upperCI)
# d5 %>% filter(Subanalysis != 'chromosome') %>% 
#   mutate(chr = gsub('chr','',chr)) %>% 
#   ggplot(aes(x=chr,y=log2(median),fill=AvZ))+
#   geom_boxplot()+ylab('Methylation (log2(f/m))')+
#   facet_grid(Tissue~Subanalysis,scales='free')+
#   theme_minimal()
cgi = read.table('Compensation_CGI_Overlapped.bed',header=FALSE)
cgi = cgi %>% dplyr::select(c(1,2,3,5,11,12,14,15,16,4))
names(cgi) = c('chr','start','end','Tissue','state','gene','cgi_start','cgi_end','ID','rnalogfm')
cgi$Subanalysis = 'cgi'
gene = read.table('Tissue_Compensation_Coordinates.bed')
gene = gene %>% dplyr::select(-c(V6,V7,V8,V9,V10))
names(gene) = c('chr','start','end','rnalogfm','Tissue','state','gene')
gene$Subanalysis = 'gene'
gene$cgi_start = NA
gene$cgi_end = NA
gene$ID = gene$gene
comp_dat = rbind(gene,cgi)
comp_dat = comp_dat %>% mutate(Tissue = gsub('spleen','Spleen',Tissue),
                               Tissue = gsub('liver','Liver',Tissue))
f = full_join(mc,comp_dat)
write.table(f,'5mC_CGI_Gene_2022-12-19.txt',quote=F,sep='\t',row.names=F)
f = read.table('5mC_CGI_Gene_2022-12-19.txt',header=TRUE)
fr = f %>% drop_na(rnalogfm) %>% drop_na(mc)

#calculate F/M
fm = fr %>% mutate(mclogfm = log2(mc))
fm = fm %>% mutate(state = gsub('partially_compensated','compensated',state))
fm = fm %>% mutate(state = gsub('compensated','DB',state),
                   state = gsub('not_compensated','Not_DB',state))
write.table(fm %>% dplyr::select(!c(row,mc_lo,mc_hi)),file='Methylation_Over_Transcripts_CGI.txt',quote=F,sep='\t',row.names=F)

#or start from here 
fm = read.table('Methylation_Over_Transcripts_CGI.txt',header=T)
#rna and 5mc values are both log2 transformed already 

#calculate rho between each 5mc and expression state (4 total comparisons)
fmrho = fm %>% group_by(chr,state,Tissue,Subanalysis) %>% 
  summarize('log2(f/m) 5mC@log2(f/m) RNA' = cor(rnalogfm,mclogfm))

#make a label for plotting 
rholab = fmrho %>% pivot_longer(!c(chr,state,Tissue,Subanalysis))  
rholab = rholab %>% separate(name,into=c('methname','genename'),sep='@',remove = T)

#plot corrleations, here only for f-m for 5mC and log2(f/m) for RNA, change those accordingly
fmp = fm %>% 
  ggplot(aes(x=mclogfm,y=rnalogfm,col=state))+
  geom_point()+
  scale_color_manual(values=viridis(3))+
  geom_vline(xintercept=0,lty=2)+
  geom_hline(yintercept=0,lty=2)+
  xlab('Methylation (log2(f/m))')+ylab('Expression (log2(f/m)')+
  geom_text(data=rholab %>% filter(state == 'DB'), mapping =(aes(x=Inf,y=Inf,label=paste0('p: ',signif(value,2)))),size=4,vjust=6,hjust=1.25)+
  geom_text(data=rholab %>% filter(state == 'not_DB'), mapping =(aes(x=Inf,y=Inf,label=paste0('p: ',signif(value,2)))),size=4,vjust=7.5,hjust=1.25)+
  facet_grid(Tissue ~ Subanalysis,scales='free')+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
fmp
pdf('Correlations_Expression_5mC.pdf',height=4,width=6,useDingbats = F)
fmp
dev.off()


#boxplots, here just show what the plot will look like, without significance test 
fmd = fm %>%
  ggplot(aes(x=state,y=mclogfm,fill=state))+
  geom_boxplot(notch = TRUE)+ylab('Methylation (log2(f/m))')+
  facet_grid(Tissue ~ Subanalysis,scales='free')+
  geom_hline(yintercept=0,lty=2)+
  geom_vline(xintercept=0,lty=2)+
  scale_fill_manual(values=c('grey70','grey20'))+
  theme_classic()

pdf('Distributions_Expression_5mC.pdf',height=4,width=6)
fmd
dev.off()

#summarize
ft = fm %>% group_by(Subanalysis,Tissue,state) %>% summarize(median = median(mc),
                                                        mean = ci(mc)[1],
                                                        lo = ci(mc)[2],
                                                        hi=ci(mc)[3])
write.table(ft, file='FM_Summaries.txt',quote=F,sep='\t',row.names=F)

#test for significance with a t-test, implement a bonferonni correction based on the tests done (by tissue, and by analysis reg v gene)
fmt = fm %>% select(ID,chr,state,Tissue,mclogfm,Subanalysis)
fmt$Tissue = factor(fmt$Tissue,levels=c('Liver','Spleen'))
stat.test <- fmt %>% ungroup %>% 
  group_by(Tissue, Subanalysis) %>%
  t_test(data =.,mclogfm ~ state) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance('p.adj')
#t-test results 
stat.test

# Tissue Subanalysis .y.     group1 group2    n1    n2 statistic    df     p p.adj p.adj.signif
# * <fct>  <chr>       <chr>   <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <dbl> <chr>       
#   1 Liver  cgi         mclogfm DB     not_DB   217   258     1.17   417. 0.243 0.972 ns          
# 2 Liver  gene        mclogfm DB     not_DB   226   268    -1.61   492. 0.107 0.428 ns          
# 3 Spleen cgi         mclogfm DB     not_DB   221   269    -1.34   482. 0.18  0.72  ns          
# 4 Spleen gene        mclogfm DB     not_DB   237   299     0.289  534. 0.773 1     ns 

#plot, should look the same as the above boxplot
ttp <- ggboxplot(
  fmt, x = "Tissue", y = "mclogfm",
  fill = "state",notch=TRUE, 
  ggtheme = theme_pubr(border = TRUE)
) +
  scale_fill_manual(values=c('grey20','grey70'))+
  xlab('')+
  geom_hline(yintercept=0,lty=2)+
  facet_grid(.~Subanalysis)+
  ylab('Methylation (log2(f/m))')
ttp
# Add statistical test p-values to the plot 
stat.test <- stat.test %>% add_xy_position(x = "Tissue",scales = 'free')
ttpa = ttp + stat_pvalue_manual(stat.test, label = "p.adj.signif",bracket.nudge.y=0.2)
pdf('Distributions_States.pdf',height=6,width=7)
ttpa
dev.off()

