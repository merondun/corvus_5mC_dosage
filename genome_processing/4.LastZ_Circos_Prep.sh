for i in $(ls *fai | sed 's/.fa.fai//g'); do
awk -v s=${i} '{OFS="\t"}{print "chr","-",$1,$1,"0", $2,s}' ${i}.fa.fai > ${i}.karyo
done

#and swap the ideograms with colors
sed -i 's/NCC-ZW/black/g' NCC-ZW.karyo
sed -i 's/PAR_genes_.flycatcher/grey/g' PAR_genes_.flycatcher.karyo
sed -i 's/PAR-GENES/orange/g' PAR-GENES.karyo
sed -i 's/PAR/purple/g' PAR.karyo
sed -i 's/Wscaffs_syntenyOnly/green/g' Wscaffs_syntenyOnly.karyo
#ensure you manually check
