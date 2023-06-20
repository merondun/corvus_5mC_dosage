for RUN in $(ls *M-bias.txt | sed 's/\..*//g' ); do 

sed -n '/CHG context (R1)/q;p' ${RUN}.*M-bias.txt | egrep -v 'methyl|=|context' \
| grep . | awk -v s=${RUN} '{OFS="\t"}{print s,$0,"R1","CpG"}' > ${RUN}.MBIAS1.tmp

sed -n '/CpG context (R2)/,$p' ${RUN}.*M-bias.txt | sed -n '/CHG context (R2)/q;p'  | egrep -v 'methyl|=|context' \
| grep . | awk -v s=${RUN} '{OFS="\t"}{print s,$0,"R2","CpG"}' > ${RUN}.MBIAS2.tmp

cat ${RUN}.MBIAS1.tmp ${RUN}.MBIAS2.tmp > ${RUN}.M-Bias.txt
rm ${RUN}*MBI*tmp
done
