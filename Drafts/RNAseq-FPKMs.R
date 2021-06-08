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

################################################################################
### DOWNLOAD
dir.create("pipeline")
setwd("pipeline")
query.exp <- GDCquery(project = "TCGA-KIRC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query.exp)
rnaseq <- GDCprepare(query = query.exp, summarizedExperiment = TRUE)
save(rnaseq, file="rnaseq-tcga.rda")
rm(query.exp)
setwd("../")


### CLEAN DATA
##! eliminar duplicados
dtmp <- colData(rnaseq) %>% 
  as.data.frame() %>% 
  select(barcode,definition,tumor_stage)
dtmp = dtmp %>% mutate(grp = case_when(
  definition == "Solid Tissue Normal" ~ "NT",
  tumor_stage == "stage i" ~ "st_1",
  tumor_stage == "stage ii" ~ "st_2",
  tumor_stage == "stage iii" ~ "st_3",
  tumor_stage == "stage iv" ~ "st_4",)
)
dtmp <- as.data.frame(dtmp  %>% drop_na(grp))
dtmp %>% group_by(grp) %>% dplyr::summarize(n())
rnaseq<- subset(rnaseq,select = colnames(rnaseq) %in% dtmp$barcode)
# exp.clinic <- colData(exp) %>% as.data.frame()
rnaseq.clinic <- dtmp
rownames(rnaseq.clinic) <- NULL
#! Checar lo de la longitud de los genes
bm <- read.csv("pipeline/biomart-20210517.csv")
bm <- subset(bm, gene_biotype="protein_coding")
rownames(bm) <- bm$ensembl_gene_id
nrow(rnaseq)
#! Clean Genes
rnaseq <- rnaseq[rownames(rnaseq) %in% bm$ensembl_gene_id, ]
nrow(rnaseq)
rm(dtmp)
# assay(rnaseq)


### NORMALIZATION
#! Spearman correlation filter
dataPre <-  TCGAanalyze_Preprocessing(rnaseq, cor.cut = 0.8,
                  filename = "pipeline/rnaseq-samplescorr.png") # samples:606
#! if samples numbers changes
rnaseq.clinic <- rnaseq.clinic[rnaseq.clinic$barcode %in% colnames(dataPre),]
dim(dataPre)
#! normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPre, geneInfo =  bm,
                                      method = "geneLength")
dataNorm <- TCGAanalyze_Normalization(tabDF = dataNorm, geneInfo =  bm,
                                      method = "gcContent")
dim(dataNorm)
# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
#! Filter by minimum mean and number of zeros
dim(dataFilt)
threshold <- round(dim(dataFilt)[2]/2)
filtro <- dataFilt[rowSums(dataFilt== 0) <= threshold, ]
dim(filtro)
filtro <- filtro[rowMeans(filtro) >= 2, ]
dim(filtro)
dataFilt <- filtro
rm(filtro, dataNorm, dataPre)

### SAVE
rnaseq <- rnaseq[rownames(rnaseq) %in% rownames(dataFilt),
       colnames(rnaseq) %in% colnames(dataFilt)]
rnaseq$grp <-  rnaseq.clinic$grp
save( rnaseq, bm, file = "pipeline/rnsaseq-clean.rda")
rm(dataFilt,rnaseq.clinic,threshold)
