source("R/tools.R")

rnaseq <- readRDS("pipeline/rnaseq-norm.rds")
mirna <- readRDS("pipeline/mirna-norm.rds")

#! DEG analysis
#! form = ~ grupo
# res.rna <- degs(rnaseq, noms=rowData(rnaseq)$gene_name, prefix = "rna") # validated by PhD. JE
deg.rna <- degs(rnaseq) # validated by PhD. JE
deg.mir <- degs(mirna)
deg.rna.an <- deg.annot.write(deg.rna, bm, prefix = "rna", c="ensembl_gene_id")
deg.mir.an <- deg.annot.write(deg.mir, bm, prefix = "mir", c="mirbase_id")
deg.all.an <- paste.degs(deg.rna.an, deg.mir.an, prefix = "all")

pval <- 1e-5

#! ctrl vs all stages Venn Diagrams
#! Cuts of LogFC & Pval
cut.rna.deg <- ups.downs(deg.rna, 1.0, pval)
cut.mir.deg <- ups.downs(deg.mir, 0.5, pval)
venn.degs(cut.rna.deg,1:4,"rna-LFC_10")
venn.degs(cut.mir.deg,1:4,"mir-LFC_10")

#! Make volcanos
source("tools.R")
volcanos(deg.rna.an, prefix = "rna", lfc=1.0, pv=pval, c="gene_name") # vaidated by PhD. JE
volcanos(deg.mir.an, prefix = "mir", lfc=0.5, pv=pval, c="mirbase_id")