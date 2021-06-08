source("tools.R")

# for(i in names(nets)){
#   print(i)
#   enrich.target.genes(nets[[i]], 
#                       toupper(cut.mir.deg$ups$stage1_ctrl),
#                       lrna = cut.rna.deg$downs$stage1_ctrl)
# }






#! Input for ARACNEe
#toAracne
# TODO: catch unique gene names
odir <- "output/expr/"
dir.create(odir, recursive = T)
rowData(mirna)$mirbase_id = rowData(mirna)$gene_name
rowData(mirna)$external_gene_name <- NA
rowData(mirna)$original_ensembl_gene_id <- NA
rowData(mirna)$gene_name = rownames(mirna)
annot <- rbind(rowData(rnaseq),rowData(mirna))
write.table(annot, file = paste0(odir,"gene-annot.tsv"), 
            row.names = F, quote=F, sep='\t')
#iterate over groups
m <- rbind(assay(rnaseq),assay(mirna))
for (grp in levels(rnaseq$grupo)){
  print(grp)
  d <- m[,rnaseq$grupo == grp]
  print(dim(d))
  fout <- paste0( odir, grp,"-expr.tsv")
  d <- tibble::rownames_to_column(as.data.frame(d), "gname")
  write.table(d, file = fout, row.names = F, quote=F, sep='\t')
}
write.table( as.data.frame(colData(mirna)), 
             file = paste0(odir,"sample-annot.tsv"), 
             row.names = F, quote=F, sep='\t')









# library(VennDiagram)

# #########
# #! ANALYSIS
# s <- t[1:4]
# ups <- list()
# downs <- list()
# for (nom in names(s)) {
#   print(nom)
#   res <- s[[nom]]
#   downs[[nom]] <- rownames(res[ (res$pvalue< 1e-5) & (res$log2FoldChange < -2.0) ,])
#   ups[[nom]] <- rownames(res[ (res$pvalue< 1e-5) & (res$log2FoldChange > 2.0) ,])
# }
# 
# colors <- c("#6b7fff", "#c3db0f", "#7FC97F", "#BEAED4", "#FDC086")
# colors <- colors[1:length(ups)]
# print(colors)
# 
# venn.diagram(x = ups,
#              # category.names = names(s),
#              filename = "pipeline/deg-ctr-all.png",
#              output=TRUE,
#              imagetype="png", 
#              scaled = FALSE,
#              col = "black",
#              fill = colors,
#              cat.col = colors,
#              cat.cex = 1,
#              cat.dist = 0.25,
#              margin = 0.15,
#              cex = 0.7
# )
# 
# all.downs <- Reduce(intersect, downs)
# all.ups <- Reduce(intersect, ups)
# all <- data.frame(downs = all.downs, ups = all.ups)
# 
# 
# l <- list(downs = all.downs, ups = all.ups )
# all <- data.frame(lapply(l, `length<-`, max(lengths(l))))





# ###################################################################3
# #!- Gene Mean
 # for(i in levels(rnaseq$grupo)){
 #  d <- rnaseq[,rnaseq$grupo == i]
 #  rownames(d) <- rowData(d)$gene_name
 #  d <- assay(d)
 #  print(paste0(i," ",mean(d[rownames(d) == "SAA2-SAA4",])))
 # }



# ###################################################################3
# #!- Making one contrast of combination in x
# contrastes <- function(x,dds1,px="all"){
#   print(green(c(x[1],x[2])))
#   res <- results(dds1,contrast=c("grupo",x[2],x[1]))
#   # print(head(res[order(res$log2FoldChange),],20))
#   fout <- paste0("pipeline/DEGS/deg-",px,"-",x[2],"-",x[1],".png")
#   print(fout)
#   png(fout,width=800,height=600)
#   print(EnhancedVolcano(res,
#                         lab = rownames(res),
#                         x = 'log2FoldChange',
#                         y = 'pvalue',
#                         pCutoff = 1e-5,
#                         FCcutoff = 2.0,
#                         subtitle = paste0(x[2],' vs ',x[1])
#   ))
#   dev.off()
# 
# }

# ###################################################################3
# #!- Plot Degs multigroup (todos vs. todos)
# degs <- function(dats, noms = NULL, form = ~ grupo, prefix = "rna"){
#   assay(dats) <- round(assay(dats))
#   #! DEGs multigroup
#   # noms es la etiqueta con la que se pintan los genes en los volcanos
#   stopifnot(!is.null(noms))
#   rownames(dats) = noms
#   print(green(paste0('[dats] before rownames duplicated: ',nrow(rnaseq))))
#   dats <- dats[!duplicated(rownames(dats)),]
#   print(green(paste0('[dats] after rownames duplicated: ',nrow(rnaseq))))
#   # dats$grupo <- as.factor(str_replace(dats$grupo," ",""))
#   # dats$grupo <- as.factor(str_replace(dats$grupo,"Control","ctr"))
#   ddsSE <- DESeqDataSet(dats, design = form)
#   dds <- DESeq(ddsSE)
# 
#   #! WORK FOR CONTRASTS
#   combi <- combn(levels(dats$grupo), 2)
#   dir.create("pipeline/DEGS")
#   t <- apply(combi,2,contrastes, dds1=dds, px=prefix)
# }


