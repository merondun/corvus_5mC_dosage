cat *_Compensation.txt | bedtools sort -i - > Tissue_Compensation_Coordinates.bed
bedtools closest -a Tissue_Compensation_Coordinates.bed -b Corvus_cornix__S_Up_H32_v5.6_polished.cgi.bed > Compensation_CGI_Overlapped.bed
