#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_5.6
genes=$genomedir/Corvus_cornix__S_Up_H32_v5.6_polished.UniqueTranscripts.bed
cgis=$genomedir/Corvus_cornix__S_Up_H32_v5.6_polished.cgi.bed
chrs=$genomedir/Corvus_cornix__S_Up_H32_v5.6_polished.CHRS.bed

for samp in $(cat Samples.list); do

#overlap 5mC positions with data, and then also with CpG islands:
$chr $pos $pos $mCG $countMETH $countUNMETH $genestart $geneend $gene $cgi
zcat $samp.CpG_5mC.cov.gz | \
    bedtools intersect -a - -b $genes -loj | \
    bedtools intersect -a - -b $cgis -loj |  \
    bedtools intersect -a - -b $chrs -loj | \
    awk '{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $10,$9-$8, $14, $13-$12, $17-$16}' | \
    bgzip -c > $samp.bed.gz
done
