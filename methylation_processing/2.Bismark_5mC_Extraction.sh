#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

SCRATCH=tmp/$SLURM_JOB_ID
mkdir tmp
mkdir $SCRATCH

#genome-based files
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_5.6
genome=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.fasta
cgi=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.cgi.bed
chrsbed=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.CHRS.bed
chrs=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.genome
ctsnp=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.TS_SNPs.bed
repeats=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.repeats.bed

#working directories
rawdata=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/WGBS
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/ATACseq_Methylation_Analyses/2021_NOV/Clip_9_Q20
workdir=${basedir}/scratch
outdir=${basedir}/5mC_raw

mkdir ${basedir}
mkdir ${workdir}
mkdir ${outdir}
mkdir ${outdir}/fastqc

RUN=$1

#cat ${rawdata}/${RUN}*__R1__*fastq.gz > $SCRATCH/${RUN}.R1.fastq.gz
#cat ${rawdata}/${RUN}*__R2__*fastq.gz > $SCRATCH/${RUN}.R2.fastq.gz

cd ${workdir}

#Trim adapters
#trim_galore --fastqc -j 10 --clip_r1 9 --clip_r2 9 --paired --quality 20 --output_dir ${workdir} ${SCRATCH}/${RUN}.R1.fastq.gz ${SCRATCH}/${RUN}.R2.fastq.gz

#Re-run fastqc on final cleaned data
#fastqc -t 8 ${workdir}/${RUN}*_val*.fq.gz

#mv ${workdir}/${RUN}*.cleaned_fastqc.* ${outdir}/fastqc

#Map reads
bismark --parallel 7 --output_dir ${workdir} --genome ${genomedir} -1 ${workdir}/${RUN}.R1_val_1.fq.gz -2 ${workdir}/${RUN}.R2_val_2.fq.gz

#Deduplicate
deduplicate_bismark --paired ${workdir}/${RUN}.R1_val_1_bismark_bt2_pe.bam --output_dir ${workdir}

#Extract methylation
bismark_methylation_extractor --parallel 10 --gzip --bedgraph --ignore_3prime_r2 1 --ignore_r2 3 --buffer_size 60G --cytosine_report --genome_folder ${genomedir} --output ${workdir} ${workdir}/${RUN}.R1_val_1_bismark_bt2_pe.deduplicated.bam

#Remove any positions overlapping a transition
bedtools subtract -A -a ${workdir}/${RUN}.*bismark.cov.gz -b ${ctsnp} > ${workdir}/${RUN}.filt_cutsite_ctsnp.tmp

#Only keep sites on the major chromosomes
bedtools intersect -wb -a ${chrsbed} -b ${workdir}/${RUN}.filt_cutsite_ctsnp.tmp | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' | sort | uniq > ${workdir}/${RUN}.filt_cutsite_ctsnp_chrs.tmp

#Divide files further into repeat-region 5mC, and then our final 5mC positions for analysis
bedtools intersect -wb -a ${repeats} -b ${workdir}/${RUN}.filt_cutsite_ctsnp_chrs.tmp | awk '{OFS="\t"}{print $5, $6, $7, $8, $9, $10}' | sort | uniq | bgzip -c > ${outdir}/${RUN}.CpG_in_Repeats.cov.gz

#Subtract repeats
bedtools subtract -A -a ${workdir}/${RUN}.filt_cutsite_ctsnp_chrs.tmp -b ${repeats} | bgzip -c > ${outdir}/${RUN}.CpG_5mC.cov.gz

### Count methylation calls between raw and filtered sets
mkdir $basedir/counts

#make a file with the first 2 columns
echo -e "${RUN}\n${RUN}\n${RUN}\n${RUN}\n${RUN}" > $basedir/counts/${RUN}.sample
echo -e "RAW\nTsSNP\nChr\nRepeatless\nRepeats" > $basedir/counts/${RUN}.field

#now for the counts
zcat $basedir/scratch/${RUN}.*bismark.cov.gz | wc -l > $basedir/counts/${RUN}.counts
cat $basedir/scratch/${RUN}.filt_cutsite_ctsnp.tmp | wc -l >> $basedir/counts/${RUN}.counts
cat $basedir/scratch/${RUN}.filt_cutsite_ctsnp_chrs.tmp | wc -l >> $basedir/counts/${RUN}.counts
zcat $basedir/5mC_raw/${RUN}.CpG_5mC.cov.gz | wc -l >> $basedir/counts/${RUN}.counts
zcat $basedir/5mC_raw/${RUN}.CpG_in_Repeats.cov.gz | wc -l >> $basedir/counts/${RUN}.counts

paste $basedir/counts/${RUN}.sample $basedir/counts/${RUN}.field $basedir/counts/${RUN}.counts > $basedir/counts/${RUN}.count
