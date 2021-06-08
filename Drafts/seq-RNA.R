library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
source("tools.R")

####################################################################
# DOWNLOAD
dir.create("pipeline")
setwd("pipeline")
query.exp <- GDCquery(project = "TCGA-KIRC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")
GDCdownload(query.exp)
rnaseq <- GDCprepare(query = query.exp, summarizedExperiment = TRUE)
save(rnaseq, file="rnaseq-exp.rda")
load(file="rnaseq-exp.rda")
rm(query.exp)
setwd("../")


####################################################################
### CLEAN DATA
# Primero las muestras
dat <- colData(rnaseq) %>% 
  as.data.frame() %>% 
  select(patient, sample, barcode,
         definition,tumor_stage,days_to_death)
dat <- seq.grouping(dat)

colDat <- DataFrame( row.names = dat$barcode,
                     barcode = dat$barcode,
                     sample = dat$sample,
                     definition = dat$definition,
                     tumor_stage = dat$tumor_stage,
                     days_to_death = dat$days_to_death,
                     grp = dat$grp
                     )

# Segundo los genes
stopifnot(sum(duplicated(rownames(rnaseq))) == 0)
m <- seq.filtro(assay(rnaseq))
# Anotacion
bm <- read.csv("pipeline/biomart.csv")
bm <- bm[bm$gene_biotype == "protein_coding",]
bm <- bm[!duplicated(bm$ensembl_gene_id),]
m <- m[rownames(m) %in% bm$ensembl_gene_id, ]
bm <- bm[bm$ensembl_gene_id %in% rownames(m),]
#checar que no haya genes duplicados

rowDat <- DataFrame(
  geneLength = bm$geneLength,
  gcContent = bm$gcContent,
  gene_name = bm$gene_name,
  chr = bm$chr,
  band = bm$band,
  gstart = bm$start_position,
  row.names = bm$ensembl_gene_id
)

m <- m[,colnames(m) %in% colDat$barcode]

# Create SE
rnaseq <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(m)),
                              colData=colDat, rowData=rowDat)
rm(m, colDat, rowDat, bm, dat)
# Salvar antes de TCGABiolinks process
save(rnaseq, file="pipeline/rnaseq-se.rda")

# ### NORMALIZATION
# #! Spearman correlation filter
# dataPre <-  TCGAanalyze_Preprocessing(rnaseq, cor.cut = 0.8,
#                                       filename = "pipeline/rnaseq-samplescorr.png") # samples:606
# #! Ajuste de las columnas si se fueron algunas muestras
# rnaseq <- rnaseq[,rnaseq$barcode %in% colnames(dataPre)]
# dim(rnaseq)
# #! normalization of genes
# dataNorm <- TCGAanalyze_Normalization(tabDF = dataPre, geneInfo = rowData(rnaseq),
#                                       method = "geneLength")
# dataNorm <- TCGAanalyze_Normalization(tabDF = dataNorm, geneInfo = rowData(rnaseq),
#                                       method = "gcContent")
# dim(dataNorm)
# # quantile filter of genes
# dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                   method = "quantile", 
#                                   qnt.cut =  0.25)
# dim(dataFilt)
# # Ajustar en los renglones, por los genes que se fueron
# rnaseq <- rnaseq[rownames(dataFilt),]









# d1 <- subset(fact.etapas,rownames(fact.etapas) %in%colnames(rnaseq))
# mydat = NOISeq::readData( assay(rnaseq) , factors = d1)
# # myPCA = dat(mydat, type = "PCA", norm=TRUE)
# myPCA = dat(mydat, type = "PCA", logtransf = F)
# explo.plot(myPCA, factor = "tumor_stage", plottype = "scores")
# d2 <- data.frame(etapa = rnaseq$grp)
# d2$etapa <-as.factor(d2$etapa)
# rownames(d2) <- colnames(rnaseq)
# # d2 <- data.frame(tumor_stage=designExp$tumor_stage, 
# #                           row.names=colnames(rnaseq))
# mydat = NOISeq::readData( assay(rnaseq) , factors = d2)
# # myPCA = dat(mydat, type = "PCA", norm=TRUE)
# myPCA = dat(mydat, type = "PCA", logtransf = F)
# explo.plot(myPCA, factor = "tumor_stage", plottype = "scores")
# rm(dataNorm,dataFilt,dataPre)

# save(rnaseq, file="pipeline/rnaseq-se.rda")
# load(file="rnaseq-se.rda")



