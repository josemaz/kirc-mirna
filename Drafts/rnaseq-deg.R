require(DESeq2)
require(EnhancedVolcano)
library(stringr)
source("tools.R")

load(file="pipeline/rnaseq-se.rda")

# DEGs multigrupo
dats <- rnaseq
rownames(dats) = rowData(dats)$gene_name
print(green(paste0('RNAseq gene: ',nrow(dats))))
dats <- dats[!duplicated(rownames(dats)),]
print(green(paste0('RNAseq gene: ',nrow(rnaseq))))
dats$grupo <- as.factor(str_replace(dats$grupo," ",""))
dats$grupo <- as.factor(str_replace(dats$grupo,"Control","ctr"))
ddsSE <- DESeqDataSet(dats, design = ~ grupo)
dds <- DESeq(ddsSE)

# c('factorName','numeratorLevel','denominatorLevel'),
res <- results(dds,contrast=c("grupo","stagei","ctr"))

# CONTRASTS
combi <- combn(levels(dats$grupo), 2)
make.degs <- function(x){
  print(green(c(x[1],x[2])))
  res <- results(dds,contrast=c("grupo",x[2],x[1]))
  # print(head(res[order(res$log2FoldChange),],20))
  fout <- paste0("pipeline/DEGS/deg-",x[2],"-",x[1],".png")
  print(fout)
  png(fout,width=800,height=600)
  print(EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  pCutoff = 1e-2,
                  FCcutoff = 0.75,
                  subtitle = paste0(x[2],' vs ',x[1])
  ))
  dev.off()
}
dir.create("pipeline/DEGS")
t <- apply(combi,2,make.degs)





# ngrp1 <- "st_1"
# ngrp2 <- "NT"
# fout <- paste0("pipeline/DEGS/deg-",ngrp1,"-",ngrp2,".png")
# png(fout,width=800,height=600)
# EnhancedVolcano(res,
#                 lab = rownames(res),
#                 x = 'log2FoldChange',
#                 y = 'pvalue',
#                 pCutoff = 1e-2,
#                 FCcutoff = 0.75,
#                 subtitle = paste0(ngrp1,' vs ',ngrp2)
# )
# dev.off()














# ### DEG fomr TCGAbiolinks
# grp1 <- dataFilt[,rnaseq.clinic[rnaseq.clinic$grp == "NT",]$barcode]
# grp2 <- dataFilt[,rnaseq.clinic[rnaseq.clinic$grp == "st_1",]$barcode]
# dataDEGs <- TCGAanalyze_DEA(mat1 = grp1, grp2,
#                             Cond1type = "Normal",
#                             Cond2type = "Tumor",
#                             # fdr.cut = 0 ,
#                             # logFC.cut = 0,
#                             method = "glmLRT")
# dataDEGs$ids <- rownames(dataDEGs)
# rownames(dataDEGs) <- NULL
# # dataDEGs <- merge[rownames(dataDEGs) %in% bm$ensembl_gene_id, ]
# dataDEGs <- merge( dataDEGs, bm, by.x = "ids", by.y = "ensembl_gene_id")
# EnhancedVolcano(dataDEGs,
#                 lab = dataDEGs$gene_name,
#                 x = 'logFC',
#                 y = 'PValue',
#                 pCutoff = 10e-6,
#                 FCcutoff = 1.5,
#                 title = 'NT vs st_1')
# rm(grp1, grp2,dataDEGs)
