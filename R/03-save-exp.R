
rnaseq <- readRDS("data/RDS/rnaseq-norm.rds")
mirna <- readRDS("data/RDS/mirna-norm.rds")


#! SAVE EXPRESSION TABLES
dout <- "data/tables/exp"
dir.create(dout, recursive=TRUE)
for (i in levels(rnaseq$grupo)){
    rna <- assay(rnaseq[,rnaseq$grupo == i])
    rna <- cbind(gname = rownames(rna), rna) # paste col gene name
    write.table(rna, file=paste0(dout,"/",i,"-rna.tsv"), 
        quote = FALSE, sep = "\t", row.names = FALSE)
    mir <- assay(mirna[,mirna$grupo == i])
    mir <- cbind(gname = rownames(mir), mir) # paste col gene name
    write.table(mir, file=paste0(dout,"/",i,"-mirna.tsv"), 
        quote = FALSE, sep = "\t", row.names = FALSE)
    all <- rbind(rna, mir)
    write.table(all, file=paste0(dout,"/",i,"-all.tsv"), 
        quote = FALSE, sep = "\t", row.names = FALSE)
}
#! SAVE ANNOTATION
write.table(rowData(rnaseq), file=paste0(dout,"/annot-rna.tsv"), 
        quote = FALSE, sep = "\t", row.names = TRUE)
write.table(rowData(mirna), file=paste0(dout,"/annot-mir.tsv"), 
        quote = FALSE, sep = "\t", row.names = TRUE)
cols <- intersect(colnames(rowData(rnaseq)),colnames(rowData(mirna)))
all.noms <- rbind(rowData(rnaseq)[cols],rowData(mirna)[cols])
write.table(all.noms, file=paste0(dout,"/annot-all.tsv"), 
        quote = FALSE, sep = "\t", row.names = TRUE)
#! SAVE CLINICAL DATA
write.table(colData(rnaseq), file=paste0(dout,"/clinical.tsv"), 
        quote = FALSE, sep = "\t", row.names = TRUE)