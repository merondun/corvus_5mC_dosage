library(dplyr)
library(tidyr)

### Begin with Chromosome 5mC data, calculate F/M first  
#for F/M ratios 
mfmeth <- read.table('Input.txt',header=TRUE)
mfmeth$mCG <- round(mfmeth$mCG,2)
head(mfmeth)

fmdat <- NULL
tissues <- c('Liver','Spleen')
for (organ in tissues) {
  cat('Working on tissue: ',organ,'\n')
  organ.df= subset(mfmeth,Tissue == organ)
  
  for (site in unique(organ.df$ID)) {
    cat('Working on region: ',site,'\n')
    
  fv=subset(organ.df, Sex == 'F' & ID == site)
  mv=subset(organ.df, Sex == 'M' & ID == site)
  
  vector.fm=c() # a vector per region
  
  for (r in 1:10000) {
    #sample, with replacement, the 5mC values from this region. Calculate the median within-sex
    fv2 = median(sample(fv$mCG,length(fv), replace=TRUE))
    mv2 = median(sample(mv$mCG,length(mv), replace=TRUE))
    #calculate F/M, but first change any zero values to 1%
    fv2[fv2==0] <- 0.01
    mv2[mv2==0] <- 0.01
    fmmed = fv2 / mv2         
    vector.fm=append(vector.fm, fmmed)
  }
    
    ### getting the CI
    med.fm=median(vector.fm)
    CI.fm=quantile(probs=c(0.025,0.975),x=vector.fm)
    fm.dat=data.frame(med.fm,t(CI.fm))
    fm.dat$Tissue <- organ; fm.dat$ID <- site
    names(fm.dat) <- c('FMmedian','FMlowerCI','FMupperCI','Tissue','ID')
    fmdat <- rbind(fmdat,fm.dat)
  }
}

mfmeth2 <- merge(mfmeth,fmdat,by=c('Tissue','ID'))
write.table(mfmeth2,file='Output.txt',quote=F,sep='\t',row.names=F)
