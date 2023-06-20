mkdir counts
echo -e "Sample\tFilter\tCount" > counts/header.txt

for i in $(cat Samples.list); do 

#make a file with the first 2 columns
echo -e "${i}\n${i}\n${i}\n${i}\n${i}" > counts/${i}.sample
echo -e "RAW\nTsSNP\nChr\nRepeatless\nRepeats" > counts/${i}.field


#now for the counts
zcat scratch/${i}.*bismark.cov.gz | wc -l > counts/${i}.counts
cat scratch/${i}.filt_cutsite_ctsnp.tmp | wc -l >> counts/${i}.counts
cat scratch/${i}.filt_cutsite_ctsnp_chrs.tmp | wc -l >> counts/${i}.counts
zcat 5mC_raw/${i}.CpG_5mC.cov.gz | wc -l >> counts/${i}.counts
zcat 5mC_raw/${i}.CpG_in_Repeats.cov.gz | wc -l >> counts/${i}.counts

paste counts/${i}.sample counts/${i}.field counts/${i}.counts > counts/${i}.count

done
cat counts/header.txt counts/*.count > 5mC_Counts.txt
rm -rf counts
