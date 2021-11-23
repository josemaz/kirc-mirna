source("tools.R")

datos <- download.paired(loading=TRUE)

bm <- read.csv("pipeline/biomart.csv")

mirna <- clean.mirna(datos$mir, bm) # only miRNA: (275,611)
rnaseq <- clean.rna(datos$rna, bm) # only RNAseq: (16227,599)

#! Patient merge
common.pats = intersect(mirna$sample,rnaseq$sample)
rnaseq <- rnaseq[,rnaseq$sample %in% common.pats]
mirna <- mirna[,match(rnaseq$sample, mirna$sample)]
#! Ajuste de informacion clinica entre rna y mir
colData(mirna) <- colData(rnaseq)

#! Normalization
rnaseq <- norm.rna(rnaseq) # validated by PhD. JE
mirna <- norm.mir(mirna) # validated by PhD. JE

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




#! Paste & save Result tables of RNA & miRNA
for (r.name in names(res.rna)){
  print(r.name)
  r.all <- rbind(res.rna[[r.name]],res.mir[[r.name]])
  fout = paste0("pipeline/DEGS/res-all-",r.name,".tsv")
  r.all <- tibble::rownames_to_column(as.data.frame(r.all),"gname")
  write.table(r.all, file = fout, row.names = F, quote=FALSE, sep='\t')
  print(fout)
}

nets <- load.nets("remote/MI","*-1e5-genmirna.tsv")
enrich.target.genes(nets[["stage1"]], toupper(cut.mir.deg$ups$stage1_ctrl),
                    lrna = cut.rna.deg$downs$stage1_ctrl)
enrich.target.genes(nets[["stage2"]], toupper(cut.mir.deg$ups$stage2_ctrl),
                    lrna = cut.rna.deg$downs$stage2_ctrl)
enrich.target.genes(nets[["stage3"]], toupper(cut.mir.deg$ups$stage3_ctrl),
                    lrna = cut.rna.deg$downs$stage3_ctrl)
enrich.target.genes(nets[["stage4"]], toupper(cut.mir.deg$ups$stage4_ctrl),
                    lrna = cut.rna.deg$downs$stage4_ctrl)
#######
enrich.target.genes(nets[["stage1"]], toupper(cut.mir.deg$downs$stage1_ctrl),
                    lrna = cut.rna.deg$ups$stage1_ctrl)
enrich.target.genes(nets[["stage2"]], toupper(cut.mir.deg$downs$stage2_ctrl),
                    lrna = cut.rna.deg$ups$stage2_ctrl)[]
enrich.target.genes(nets[["stage3"]], toupper(cut.mir.deg$downs$stage3_ctrl),
                    lrna = cut.rna.deg$ups$stage3_ctrl)
enrich.target.genes(nets[["stage4"]], toupper(cut.mir.deg$downs$stage4_ctrl),
                    lrna = cut.rna.deg$ups$stage4_ctrl)
###########################
source("tools.R")
enrich.target.genes(nets[["stage1"]], toupper(cut.mir.deg$downs$stage1_ctrl),
                    lrna = cut.rna.deg$ups$stage1_ctrl, bm)
enrich.target.genes(nets[["stage2"]], toupper(cut.mir.deg$downs$stage2_stage1),
                    lrna = cut.rna.deg$ups$stage2_stage1, bm)
enrich.target.genes(nets[["stage3"]], toupper(cut.mir.deg$downs$stage3_stage2),
                    lrna = cut.rna.deg$ups$stage3_stage2, bm)
enrich.target.genes(nets[["stage4"]], toupper(cut.mir.deg$downs$stage4_stage3),
                    lrna = cut.rna.deg$ups$stage4_stage3, bm)

sets <- list()
for (n in names(nets)){
  print(n)
  inter <- paste0(nets[[n]]$Source,"_",nets[[n]]$Target)
  sets[[n]] <- inter
}
colors <- c("#6b7fff", "#c3db0f", "#7FC97F", "#BEAED4", "#FDC086")
colors <- colors[1:length(sets)]
venn.diagram(x = sets,
             # category.names = names(s),
             filename = "output/venn-inter.png",
             output=TRUE,
             imagetype="png",
             scaled = FALSE,
             col = "black",
             fill = colors,
             cat.col = colors,
             cat.cex = 1,
             cat.dist = 0.25,
             margin = 0.15,
             cex = 0.7,
             main = "Interacciones en etapas"
)
all <- Reduce(intersect, sets[2:5])
only.etapas <- Reduce(intersect, sets[2:5])
