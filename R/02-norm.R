source("R/tools.R")
library(tibble)

#! INPUT
rnaseq <- readRDS("data/RDS/rnaseq-clean.rds")
mirna <- readRDS("data/RDS/mirna-clean.rds")

#! NORMALIZATION
dout <- "data/plots/QC"
dir.create(dout, recursive=TRUE)
rnaseq <- norm.rna(rnaseq,dout) # validated by PhD. JE
mirna <- norm.mir(mirna,dout) # validated by PhD. JE

#! SAVE RDS
saveRDS(rnaseq, file = "data/RDS/rnaseq-norm.rds")
saveRDS(mirna, file = "data/RDS/mirna-norm.rds")


