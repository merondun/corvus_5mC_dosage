# Extract the PAR region based on coordinates defined in PAR.bed
seqkit subseq --bed PAR.bed chrZ.fa > PAR.fa

# Extract PAR genes based on coordinates defined in Genes-PAR.bed
seqkit subseq --bed Genes-PAR.bed chrZ.fa > PAR-GENES.fa

# Clean up headers in PAR-GENES.fa file so only gene name is left
sed -i 's/chrZ.*:. //' PAR-GENES.fa

# Extract only Z and W sequences from the NCC assembly, and clean up names
seqtk subseq GCF_009650955.1_bCorMon1.pri_genomic.fna NCC-W.name > NCC-ZW.fa
sed -i 's/NC_045510.1.*sequence/NCC_chrW/' NCC-ZW.fa
sed -i 's/NC_045511.1.*sequence/NCC_chrZ/' NCC-ZW.fa

# Clean up FLY sequences, remove everything after the first pipe character
sed -i 's/|.*//' PAR_genes_.flycatcher.fa

# Define the path to the PAR region
PAR=PAR.fa
PARG=PAR-GENES.fa
Wsyn=Wscaffs_syntenyOnly.fa
Flycatch=PAR_genes_.flycatcher.fa
NCC=NCC-ZW.fa

# A. Align PAR region to the W synteny region
lastz $PAR $Wsyn --notransition --step=20 --nogapped --format=general > lastz/PAR-Wsyn.out

# B. Align PAR region to the NCC region
lastz $PAR $NCC --notransition --step=20 --nogapped --format=general > lastz/PAR-NCC.out

# C. Align PAR region to the Flycatcher region
lastz $PAR $Flycatch --notransition --step=20 --nogapped --format=general > lastz/PAR-FLY.out

# D. Align Flycatcher region to the NCC region
lastz $Flycatch[multiple] $NCC --notransition --step=20 --nogapped --format=general > lastz/FLY-NCC.out

# E. Align PAR gene region to the W synteny region
lastz $PARG[multiple] $Wsyn --notransition --step=20 --nogapped --format=general > lastz/PARG-Wsyn.out

# F. Align PAR gene region to the Flycatcher region
lastz $PARG[multiple] $Flycatch --notransition --step=20 --nogapped --format=general > lastz/PARG-FLY.out

# G. Align PAR gene region to the NCC region
lastz $PARG[multiple] $NCC --notransition --step=20 --nogapped --format=general > lastz/PARG-NCC.out

# Extract specific columns from each output file and save them into separate bed files
for i in $(ls *out | sed 's/.out//g'); do
    # Grab name/start/end/name2/start2/end2 and add a custom link color
    # Filter only for alignment length > 1000 and sequence identity > 70
    awk -v s=${i} '{OFS="\t"}{print $2, $5, $6, $6-$5, $7,$10, $11,$11-$10,$12,$13,$14, $15, "color="s}' ${i}.out | grep -v 'zstart' | awk '$10 > 70 && $4 > 1000' | cut -f4,8,9,10,11,12 --complement > ${i}.1KB-P70.bed
    awk -v s=${i} '{OFS="\t"}{print $2, $5, $6, $6-$5, $7,$10, $11,$11-$10,$12,$13,$14, $15, "color="s}' ${i}.out | grep -v 'zstart' | cut -f4,8,9,10,11,12 --complement > ${i}.all.bed
done

# Assign specific colors to each different alignment pair
sed -i 's/FLY-NCC/black/g' FLY-NCC.bed
sed -i 's/PAR-FLY/grey/g' PAR-FLY.bed
sed -i 's/PARG-FLY/orange/g' PARG-FLY.bed
sed -i 's/PARG-NCC/purple/g' PARG-NCC.bed
sed -i 's/PARG-Wsyn/green/g' PARG-Wsyn.bed
sed -i 's/PAR-NCC/yellow/g' PAR-NCC.bed
sed -i 's/PAR-Wsyn/blue/g' PAR-Wsyn.bed

# Loop over all fai files in the current directory, create a new file with the extension .karyo which will create the circos karyotype files
for i in $(ls *fai | sed 's/.fa.fai//g'); do
    # In the new .karyo file, for each line in the .fai file,
    # print "chr", "-", the first column from the .fai file (sequence name),
    # the first column again, "0", the second column from the .fai file (sequence length), 
    # and the sequence name obtained from the file name
    awk -v s=${i} '{OFS="\t"}{print "chr","-",$1,$1,"0", $2,s}' ${i}.fa.fai > ${i}.karyo
done

# Change the sequence name in each .karyo file to a color name
# This is done in preparation for visualization in a tool like Circos, 
# which can use this color name to color the corresponding ideogram
sed -i 's/NCC-ZW/black/g' NCC-ZW.karyo
sed -i 's/PAR_genes_.flycatcher/grey/g' PAR_genes_.flycatcher.karyo
sed -i 's/PAR-GENES/orange/g' PAR-GENES.karyo
sed -i 's/PAR/purple/g' PAR.karyo
sed -i 's/Wscaffs_syntenyOnly/green/g' Wscaffs_syntenyOnly.karyo

