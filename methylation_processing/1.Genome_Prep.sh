#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

#genome-based files
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_5.6
genome=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.fasta
cgi=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.cgi.bed
repeats=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.repeats.bed
annotation=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.annotation.bed
rnaseq=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.RNAseq.bed
atacseq=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.ATACseq.bed
rnasub=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.SUB.RNAseq.bed
atacsub=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.SUB.ATACseq.bed
chrsbed=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.CHRS.bed
chrs=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.genome
ctsnp=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.CTSNPs.bed
allsnp=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.AllSNPs.bed

#sequence-based files
rawdata=/crex/proj/uppstore2019047/nobackup/wgbs/corvus/raw

#sample-id list
list=/crex/proj/uppstore2019047/nobackup/wgbs/corvus/lists/All_Samples.list

bismark_genome_preparation ${genomedir}
