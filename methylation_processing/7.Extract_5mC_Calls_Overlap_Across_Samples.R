source('~/modules/R_Functions.R')
library(reshape2)
library(matrixStats)
library(tidyverse)

#grab all files
files <- list.files(path=".", pattern="*.bed.gz", full.names=TRUE, recursive=FALSE)
allnames <- NULL

#loop through them.. CHROMOSOME
for (samp in files) {

  cat('Reading in file: ',samp,'\n')
  #read data
  ztab <- gzfile(samp,'rt');  tab <- read.table(ztab,header=F);
  #only keep positions with more than 9x coverage
  tab <- subset(tab, (V5 + V6) > 9);  tab$site <- paste0(tab$V1,'_',tab$V2); tab$V4 <- tab$V4/100

  #extract name which will be our new variable
  name <- gsub('./','',samp); name <- gsub('.bed.gz','',name)

  #chromosome file
  tab2 <- tab %>% group_by(V1,V11) %>% summarize(mCG = sum(V5)/(sum(V5)+sum(V6)),
                                                 mC = sum(V5),
                                                 mA = sum(V6)+sum(V5))
  names(tab2) <- c('ID','length',paste0('mCG.',name),paste0('mC.',name),paste0('mCA.',name));
  tab2 <- tab2 %>% mutate(AvZ = ifelse(ID == 'chrZ', 'chrZ','Autosome')); tab2$chr <- tab2$ID; tab2$Subanalysis <- 'chromosome'

  #gene file
  tabG <- subset(tab,V7 !='.')
  tab3 <- tabG %>% group_by(V7,V8,V1) %>% summarize(mCG = sum(V5)/(sum(V5)+sum(V6)),
                                                 mC = sum(V5),
                                                 mA = sum(V6)+sum(V5))
  names(tab3) <- c('ID','length','chr',paste0('mCG.',name),paste0('mC.',name),paste0('mCA.',name));
  tab3 <- tab3 %>% mutate(AvZ = ifelse(chr == 'chrZ', 'chrZ','Autosome')); tab3$Subanalysis <- 'gene'

  #cgi file
  tabC <- subset(tab,V9 !='.')
  tab4 <- tabC %>% group_by(V9,V10,V1) %>% summarize(mCG = sum(V5)/(sum(V5)+sum(V6)),
                                                     mC = sum(V5),
                                                        mA = sum(V6)+sum(V5))
  names(tab4) <- c('ID','length','chr',paste0('mCG.',name),paste0('mC.',name),paste0('mCA.',name));
  tab4 <- tab4 %>% mutate(AvZ = ifelse(chr == 'chrZ', 'chrZ','Autosome')); tab4$Subanalysis <- 'cgi'

  #merge them
  tab5 <- rbind(tab2,tab3,tab4)
  cat('Saved to variable: ',name,'\n')
  #assign it to a new variable
  assign(paste0(name),tab5)
  #add this sample to our total list
  allnames <- unique(c(allnames,name))
}

#get names
vars <- mget(allnames)
# merge the points, keeping all
masterWGBS <- Reduce(function(x,y) merge(x = x, y = y, by=c('ID','length','AvZ','chr','Subanalysis'),all=TRUE),vars)

write.table(masterWGBS,file='Step1.txt',quote=F,sep='\t',row.names=F)

#filter according to missingness, retain sites with no missing data, according to tissue
samps <- c('S_Up_H60_L_ADL_F','S_Up_H60_M_ADL_F',
           'S_Up_H59_L_ADL_M','S_Up_H59_M_ADL_M',
           'D_Ko_C45_M_ADL_F','D_Ko_C45_L_ADL_F',
           'D_Ko_C31_M_ADL_M','D_Ko_C31_L_ADL_M')

sampM <- samps[grepl('_M_AD',samps)]
#remove sites with missing data specific by tissue
covM <- masterWGBS[rowSums(is.na(masterWGBS[grepl('^mCG..*_M_ADL', names(masterWGBS))])) == 0, ]
#convert to long form
mdat <- NULL
for (samp in sampM) {
  cat('Subsetting: ',samp,'\n')
  v1 <- covM[grepl(paste0('ID|length|AvZ|chr|Subanalysis|',samp),names(covM))]
  v1$Sample <- samp
  names(v1) <- gsub(paste0('.',samp),'',names(v1))
  sex <- ifelse(grepl('ADL_F',samp) == TRUE,'F','M')
  v1$Sex <- sex
  mdat <- rbind(mdat,v1)
}
mdat$Tissue <- 'Spleen'

#and liver
sampL <- samps[grepl('_L_AD',samps)]
#remove sites with missing data specific by tissue
covL <- masterWGBS[rowSums(is.na(masterWGBS[grepl('^mCG..*_L_ADL', names(masterWGBS))])) == 0, ]
#convert to long form
ldat <- NULL
for (samp in sampL) {
  cat('Subsetting: ',samp,'\n')
  v1 <- covL[grepl(paste0('ID|length|AvZ|chr|Subanalysis|',samp),names(covL))]
  v1$Sample <- samp
  names(v1) <- gsub(paste0('.',samp),'',names(v1))
  sex <- ifelse(grepl('ADL_F',samp) == TRUE,'F','M')
  v1$Sex <- sex
  ldat <- rbind(ldat,v1)
}
ldat$Tissue <- 'Liver'

#merge
cdat <- rbind(mdat,ldat)
cdat$Sample <- gsub('_M_ADL_M','',cdat$Sample)
cdat$Sample <- gsub('_L_ADL_M','',cdat$Sample)
cdat$Sample <- gsub('_M_ADL_F','',cdat$Sample)
cdat$Sample <- gsub('_L_ADL_F','',cdat$Sample)

write.table(cdat,file='mCG-10x_LongForm_ATAC-22May9.txt',quote=F,sep='\t',row.names=F)

