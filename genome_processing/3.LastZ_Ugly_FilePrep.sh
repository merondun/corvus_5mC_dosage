#grab name/start/end/name2/start2/end2 and add a custom link color
for i in $(ls *out | sed 's/.out//g'); do

awk -v s=${i} '{OFS="\t"}{print $2, $5, $6, $6-$5, $7,$10, $11,$11-$10,$12,$13,$14, $15, "color="s}' ${i}.out | \
    grep -v 'zstart' | awk '$10 > 70 && $4 > 1000' | cut -f4,8,9,10,11,12 --complement \
    > ${i}.1KB-P70.bed
awk -v s=${i} '{OFS="\t"}{print $2, $5, $6, $6-$5, $7,$10, $11,$11-$10,$12,$13,$14, $15, "color="s}' ${i}.out | \
    grep -v 'zstart' | cut -f4,8,9,10,11,12 --complement \
    > ${i}.all.bed
done

#now swap our IDs with specific colors
sed -i 's/WScaff'
sed -i 's/FLY-NCC/black/g' FLY-NCC.bed
sed -i 's/PAR-FLY/grey/g' PAR-FLY.bed
sed -i 's/PARG-FLY/orange/g' PARG-FLY.bed
sed -i 's/PARG-NCC/purple/g' PARG-NCC.bed
sed -i 's/PARG-Wsyn/green/g' PARG-Wsyn.bed
sed -i 's/PAR-NCC/yellow/g' PAR-NCC.bed
sed -i 's/PAR-Wsyn/blue/g' PAR-Wsyn.bed
