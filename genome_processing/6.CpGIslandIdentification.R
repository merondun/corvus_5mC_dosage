library(makeCGI);
.CGIoptions=CGIoptions();
.CGIoptions$rawdat.type="txt";
.CGIoptions$species="Corvus_cornix__S_Up_H32_v5.6_polished";
.CGIoptions;
makeCGI(.CGIoptions);
5000   #window size

#in bash afterwards
#awk '($4 > 250) && ($7 > 0.6) && ($8 > 0.75){print $1, $2, $3}' CGI-Corvus_cornix__S_Up_H32_v5.6_polished.txt | sed '1d' | tr ' ' \\t  > ../../Corvus_cornix__S_Up_H32_v5.6_polished.cgi.bed
