source("tools.R")

datos <- download.paired(loading=TRUE)

bm <- read.csv("pipeline/biomart.csv")

mirna <- clean.mirna(datos$mir, bm) # only miRNA: (275,611)
rnaseq <- clean.rna(datos$rna, bm) # only RNAseq: (16227,599)

common.pats = intersect(mirna$sample,rnaseq$sample)

rnaseq <- rnaseq[,rnaseq$sample %in% common.pats]
mirna <- mirna[,match(rnaseq$sample, mirna$sample)]
colData(mirna) <- colData(rnaseq)

# Normalization
norm.rna(rnaseq)
norm.mir(mirna)


