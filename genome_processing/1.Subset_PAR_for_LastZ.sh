#cat PAR.bed
#chrZ    1       688000
#extract PAR sequence
seqkit subseq --bed PAR.bed chrZ.fa > PAR.fa

#extract PAR genes
seqkit subseq --bed Genes-PAR.bed chrZ.fa > PAR-GENES.fa
#Clean up headers so only gene name is left
sed -i 's/chrZ.*:. //' PAR-GENES.fa

#Grab only Z+W from NCC assembly, and clean up names
seqtk subseq GCF_009650955.1_bCorMon1.pri_genomic.fna NCC-W.name > NCC-ZW.fa
sed -i 's/NC_045510.1.*sequence/NCC_chrW/' NCC-ZW.fa 
sed -i 's/NC_045511.1.*sequence/NCC_chrZ/' NCC-ZW.fa 

#clean up FLY, remove everything after the first pipe
sed -i 's/|.*//' PAR_genes_.flycatcher.fa
