source("R/tools.R")

datos <- download.paired(loading=TRUE)

bm <- read.csv("data/tables/biomart.csv")

mirna <- clean.mirna(datos$mir, bm) # only miRNA: (275,611)
rnaseq <- clean.rna(datos$rna, bm) # only RNAseq: (16227,599)

#! Patient merge
common.pats = intersect(mirna$sample,rnaseq$sample)
rnaseq <- rnaseq[,rnaseq$sample %in% common.pats]
mirna <- mirna[,match(rnaseq$sample, mirna$sample)]
#! Ajuste de informacion clinica entre rna y mir
colData(mirna) <- colData(rnaseq)

saveRDS(rnaseq, file = "data/RDS/rnaseq-clean.rds")
saveRDS(mirna, file = "data/RDS/mirna-clean.rds")