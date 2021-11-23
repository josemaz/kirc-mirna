source("R/tools.R")
library(tibble)


rnaseq <- readRDS("data/RDS/rnaseq-clean.rds")
mirna <- readRDS("data/RDS/rmirna-clean.rds")

#! Normalization
rnaseq <- norm.rna(rnaseq) # validated by PhD. JE
mirna <- norm.mir(mirna) # validated by PhD. JE

saveRDS(rnaseq, file = "data/RDS/rnaseq-norm.rds")
saveRDS(mirna, file = "data/RDS/mirna-norm.rds")

dout <- "data/tables/exp"
dir.create(dout)

for (i in levels(rnaseq$grupo)){
    e <- assay(rnaseq[,rnaseq$grupo == i])
    e <- cbind(gname = rownames(e), e)
    write.table(e, file=paste0(dout,"/",i,".tsv"), 
        quote = FALSE, sep = "\t", row.names = FALSE)
}