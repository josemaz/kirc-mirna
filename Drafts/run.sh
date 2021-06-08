[ ! -d pipeline ] && mkdir pipeline
Rscript biomart.R
Rscript RNAseq.R