library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)
library(tidyverse)
library(EnhancedVolcano)
library(VennDiagram)
# library(ggfortify)
# library(ggbiplot)
library(NOISeq)
library(dplyr)



#####################################################################
#### miRNAs
#####################################################################

# query <- GDCquery(project = "TCGA-KIRC", 
#                       data.category = "Transcriptome Profiling",
#                       data.type = "Isoform Expression Quantification",
#                       workflow.type = "BCGSC miRNA Profiling",
#                       # barcode = c("TCGA-CJ-5680")
#                   )
# GDCdownload(query)
# miRNA.iso <- GDCprepare(query)
# save(data, file="miRNAs-isoform.rda")
# load(file="miRNAs-isoform.rda")

# query <- GDCquery(project = "TCGA-KIRC",
#                   data.category = "Transcriptome Profiling",
#                   data.type = "miRNA Expression Quantification",
#                   workflow.type = "BCGSC miRNA Profiling")
# GDCdownload(query)
# miRNA.exp <- GDCprepare(query)
# save(miRNA.exp, file="miRNAs-exp.rda")
load(file="miRNAs-exp.rda")


#!## PARING
#! miRNA
mirna <- select(miRNA.exp,starts_with("read_count_"))
mirna <- cbind(miRNA.exp$miRNA_ID,mirna)
# colnames(mirna)[1] <- "mirID"
rownames(mirna) <- mirna[,1]
mirna <- data.matrix(mirna[,-1])
colnames(mirna) <- gsub('read_count_', '', colnames(mirna), fixed=TRUE)
df1 <- data.frame(mir_sample = colnames(mirna), stringsAsFactors=FALSE)
df1$mir_id <- 
  str_extract(df1$mir_sample, "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{3}")
df1 <- df1[!duplicated(df1$mir_id),] # muestras unicas de mirnas
# RNAseq
df2 <- data.frame(rna_sample = exp.clinic$barcode, stringsAsFactors=FALSE)
df2$rna_id <-
  str_extract(df2$rna_sample, "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{3}")
df2<- df2[!duplicated(df2$rna_id),] # muestras unicas de RNAseq
# #! Merge
# venn.diagram(
#   list(mirna = df1$mir_id, rnaseq = df2$rna_id),
#   fill = c("green", "blue"),
#   "cases.tiff"
# )
cases = merge(df1,df2,by.x="mir_id",by.y="rna_id") # 583 paired cases
cases = rename(cases,'id'='mir_id')
cases <- merge( cases, exp.clinic, by.x='rna_sample', by.y='barcode')
# #! checar que todos sean casos unicos
# #! Crate data expression object
# expr <- list(
#   mir = mirna[,colnames(mirna) %in% cases$mir_sample ],
#   rna = dataFilt[,colnames(dataFilt) %in% cases$rna_sample]
# )
colData <- DataFrame(
  definition = cases$definition,
  Etapa = cases$grp,
  row.names = cases$mir_sample,
  barcode = cases$mir_sample )
#! Ajusta columnas de expresion
mirna <- mirna[,colnames(mirna) %in% cases$mir_sample]

#! Annotation mirna
dim(mirna)
# df <- merge( assay(mirna), bm, by.x=0, by.y='mirbase_id')
df <- bm[bm$mirbase_id %in% rownames(mirna),]
df <- df[!duplicated(df$mirbase_id),] # (1778) borra mirbase duplicados
rowData <- DataFrame(
  ensembl_gene_id = df$ensembl_gene_id,
  geneLength = df$geneLength,
  gcContent = df$gcContent,
  gene_name = df$gene_name,
  chr = df$chr,
  band = df$band,
  gstart = df$start_position,
  row.names = df$mirbase_id
)
#ajusta renglones en mirna
mirna <- as.matrix(mirna[rownames(mirna) %in% df$mirbase_id,])
# Create SE
mirna <- SummarizedExperiment(assays=SimpleList(counts=mirna),
                              colData=colData,
                              rowData=rowData)
rm(df1,df2,colData,rowData)



# MIRNA analysis
# Spearman correlation filter
dim(mirna)   # mirnas:1881x583
# minimo para empezar acortar muestras
dataPre <-  TCGAanalyze_Preprocessing(mirna, cor.cut = 0.7) #~.75
dim(dataPre) # 1778x138
# normalization of genes
# df<- as.data.frame(rowData(mirna))
# ggplot(df, aes(x=geneLength)) + geom_density()
# dataNorm <- TCGAanalyze_Normalization(tabDF = dataPre, geneInfo =  rowData(mirna),
#                                       method = "geneLength")
# dataNorm <- TCGAanalyze_Normalization(tabDF = dataPre, geneInfo =  rowData(mirna),
#                                       method = "gcContent")
# dim(dataNorm)
# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataPre,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
dim(dataFilt)
# Filter by minimum mean and number of zeros
threshold <- round(dim(dataFilt)[2]/2)
filtro <- dataFilt[rowSums(dataFilt== 0) <= threshold, ]
dim(filtro)
filtro <- filtro[rowMeans(filtro) >= 2, ]
dim(filtro)
dataFilt <- filtro
rm(filtro, dataNorm, dataPre)
### DEG
ngrp1 <- "NT"
ngrp2 <- "st_1"
grp1 <- dataFilt[,mirna[,mirna$Etapa == ngrp1]$barcode]
grp2 <- dataFilt[,mirna[,mirna$Etapa == ngrp2]$barcode]
dataDEGs <- TCGAanalyze_DEA(mat1 = grp1, mat2 = grp2,
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            # fdr.cut = 0 ,
                            # logFC.cut = 0,
                            method = "glmLRT")
dataDEGs$ids <- rownames(dataDEGs)
rownames(dataDEGs) <- NULL
# dataDEGs <- merge[rownames(dataDEGs) %in% bm$ensembl_gene_id, ]
dim(dataDEGs)
df <- as.data.frame(rowData(mirna))
dataDEGs <- merge( dataDEGs, df, by.x = "ids", by.y = 0)
dim(dataDEGs)
EnhancedVolcano(dataDEGs,
                lab = dataDEGs$gene_name,
                x = 'logFC',
                y = 'PValue',
                pCutoff = 1e-2,
                FCcutoff = 0.75,
                title = paste0(ngrp1,' vs ',ngrp2)
                )
## PCA
df <- as.data.frame(t(cbind(grp1,grp2)))
cancer.pca <- stats::prcomp(df)
cancer.pca<- as.data.frame(cancer.pca$x)
cancer.pca$group <- c(rep("blue", ncol(grp1)), rep("red", ncol(grp2)))
ggplot(cancer.pca,aes(x=PC1,y=PC2,color=group )) + geom_point()

# BATCH
df <- data.frame(Batch = c(rep("blue", ncol(grp1)), rep("red", ncol(grp2))) )
df$Batch <- as.factor(df$Batch)
rownames(df) <- colnames(cbind(grp1,grp2))
countsCorrected <- TCGAbatch_Correction(cbind(grp1,grp2),UnpublishedData = TRUE,AnnotationDF =df)

## PCA
grp1 <- countsCorrected[,mirna[,mirna$Etapa == ngrp1]$barcode]
grp2 <- countsCorrected[,mirna[,mirna$Etapa == ngrp2]$barcode]
df <- as.data.frame(t(cbind(grp1,grp2)))
cancer.pca <- stats::prcomp(df)
cancer.pca<- as.data.frame(cancer.pca$x)
cancer.pca$group <- c(rep("blue", ncol(grp1)), rep("red", ncol(grp2)))
ggplot(cancer.pca,aes(x=PC1,y=PC2,color=group )) + geom_point()

### ARSYN 1
myfact <- data.frame(grp = c(rep("blue", ncol(grp1)), rep("red", ncol(grp2))))
myfact$grp <- as.factor(myfact$grp)
cuentas <- cbind(grp1,grp2)
rownames(myfact) <- colnames(cuentas)
mydata = readData( cbind(grp1,grp2), 
                   factors = myfact)

myPCA = dat(mydata, type = "PCA")
explo.plot(myPCA, factor = "grp")

mydata2corr1 = ARSyNseq(mydata, factor = "grp", 
                        batch = TRUE, norm = "n",  logtransf = FALSE)

myPCA = dat(mydata2corr1, type = "PCA")
explo.plot(myPCA, factor = "grp")


### ARSYN 2
grps <- data.frame(sta = as.factor(mirna$Etapa))
grps = grps %>% mutate(colores = case_when(
  sta == "NT" ~ "red",
  sta == "st_1" ~ "blue",
  sta == "st_2" ~ "green",
  sta == "st_3" ~ "yellow",
  sta == "st_4" ~ "black",)
)
rownames(grps) = colnames(mirna)
mydata = readData( dataFilt, factors = grps)

myPCA = dat(mydata, type = "PCA")
explo.plot(myPCA, factor = "sta")

mydata2corr1 = ARSyNseq(mydata, factor = "sta", 
                        batch = TRUE, norm = "n",  logtransf = FALSE)

myPCA = dat(mydata2corr1, type = "PCA")
explo.plot(myPCA, factor = "sta")

rm(grp1, grp2,dataDEGs)




