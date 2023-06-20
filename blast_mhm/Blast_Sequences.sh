#Location of Genome
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_5.6/Corvus_cornix__S_Up_H32_v5.6_polished.fasta
#Location of TARGET
mhm=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/ATACseq_Methylation_Analyses/2021_NOV/Clip_9_Q20/MHM_Sequences.fa
#Location of your desired blastDB (directory must be created)
bdb=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/modules/blastdb/Corvus5.6

echo "$(tput setaf 1) $(tput setab 7) Creating Blast DB...$(tput sgr 0)"
makeblastdb -in $genome -parse_seqids -dbtype nucl -title corvus5.6 -out $bdb

echo "$(tput setaf 1) $(tput setab 7) Blast away...$(tput sgr 0)"
blastn -query $mhm -db $bdb \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand' \
       -out target.blast

#create header
echo 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand' | tr ' ' '\t' > target.head
cat target.head target.blast > MHMchicken_5.6corvus.blast

#qseqid                                     sseqid  pident  length  mismatch  gapopen  qstart  qend  sstart     send       evalue    bitscore  sstrand
#MHM1_27.140_27.177Mb_2780bp_unit_13_times  chr2    81.013  79      11        3        2356    2434  147645094  147645168  2.13e-06  60.2      plus
#MHM1_27.237_27.270Mb_2781bp_unit_12_times  chr2    81.013  79      11        3        761     839   147645168  147645094  2.14e-06  60.2      minus
