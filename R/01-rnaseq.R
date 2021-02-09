library(TCGAbiolinks)
library(NOISeq)
library(edgeR)
library(DESeq)
library(biomaRt)
library(crayon)
library(SummarizedExperiment)
library(VennDiagram)
library(EDASeq)

raiz <- getwd()

# sessionInfo() # show library versions
packageVersion("TCGAbiolinks") # show library versions

# Setting Global Variables
EXP.NAME <- "kirc"
OUTDIR <- "Output"
dir.create(OUTDIR)
WORKDIR <- "pipeline/TCGA"
dir.create(WORKDIR, recursive = TRUE)
setwd(WORKDIR)
PLOTSDIR <- "Plots"
dir.create(PLOTSDIR)
PREDIR <- paste(PLOTSDIR, "QC_PRE", sep = "/")
dir.create(PREDIR)
POSTDIR <- paste(PLOTSDIR, "QC_POST", sep = "/")
dir.create(POSTDIR)
w <- 1024 #Resolucion de los plots
h <- 1024  #Resolucion de los plots
p <- 24   #Resolucion de los plots



#####################################################
cat(red("0.5 Download Data"),"\n")
query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query)
#! gene annotation GRCh38.p13
data.raw <- GDCprepare(query)
save(data.raw, file=EXP.NAME%+%"-expr.RDATA", compress="xz")
load(file=EXP.NAME%+%"-expr.RDATA") 



######################################################
cat(red("1.0 Loding Data"),"\n")
print(head(data.raw)) 
#!                 TCGA-B3-3925-01A-02R-1351-07 TCGA-V9-A7HT-01A-11R-A33Z-07
#! ENSG00000000003                         3079                         2435
#! Clean Assay (Data)
data.raw <- data.raw[, data.raw$tumor_stage != "not reported"]
#! Get expression
data.expr=assay(data.raw)
head(rownames(data.expr))
designExp=colData(data.raw)
#!                              days_to_last_follow_up tumor_stage
#!                                           <integer> <character>
#! TCGA-KO-8415-01A-11R-2315-07                   2939     stage i
#! TCGA-KN-8423-01A-11R-2315-07                    788     stage i
coltumors = designExp$patient[designExp$sample_type_id=="01"] #Primary solid Tumor
colnorms  = designExp$patient[designExp$sample_type_id=="11"] #Solid tissue normal
designExp$tumor_stage[designExp$sample_type_id=="11"] <- "ctrl"
#! check the intersection between normal & tumor samples
# flog.threshold(ERROR)
venn.diagram(x = list(A=coltumors, B=colnorms), filename = "Plots/Venn.tiff",
col = "transparent", fill = c("cornflowerblue","green"), 
alpha = 0.50,cex = 1.5, fontface = "bold", label.col="black", 
cat.cex = 1.5,margin = 0.1,category.names=c("tumor","normal"))



######################################################
cat(red("2.0 Annotation"),"\n")
# ! mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",mirror ="uswest")
# ! listAttributes(mart)
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
Annot=getBM(attributes = c("ensembl_gene_id", "percentage_gene_gc_content", "gene_biotype",
  "start_position","end_position","hgnc_id","hgnc_symbol","chromosome_name"),
  filters = "ensembl_gene_id", values=rownames(data.expr), mart=mart)
Annot$lenght=abs(Annot$end_position-Annot$start_position)
save(Annot, file="Annot.RDATA", compress="xz")
load(file="Annot.RDATA") 



#####################################################
cat(red("3.0 Clean Genes"),"\n")
#! keep no duplicates gene of protein coding
# Annot=Annot[Annot$gene_biotype=="protein_coding"&Annot$hgnc_symbol!="",]
# Annot$gene_biotype=="protein_coding"&
Annot=Annot[Annot$hgnc_symbol!="",]
Annot=Annot[!duplicated(Annot$ensembl_gene_id),]
exprots_hgnc=data.expr[rownames(data.expr)%in%Annot$ensembl_gene_id,]
# nrow(exprots_hgnc) # 19222 transcripts
#! sum duplicated probes = 2 probes mapping to the same hgnc_id
Annot$hgnc_id[duplicated(Annot$hgnc_id)]
#! [1] "HGNC:30046"
Annot[Annot$hgnc_id=="HGNC:30046",]
#!       ensembl_gene_id percentage_gene_gc_content   gene_biotype start_position
#! 41657 ENSG00000254093                      43.17 protein_coding       10764961
#! 44490 ENSG00000258724                      43.86 protein_coding       10725399
#!       end_position    hgnc_id hgnc_symbol
#! 41657     10839884 HGNC:30046       PINX1
#! 44490     10839847 HGNC:30046       PINX1
#! Rows duplicated
rowdups = rownames(exprots_hgnc)%in%Annot$ensembl_gene_id[Annot$hgnc_id=="HGNC:30046"]
if(is.matrix(exprots_hgnc[rowdups,])){
  exprots_hgnc[rowdups,][1,] = colSums(exprots_hgnc[rowdups,])
  exprots_hgnc=exprots_hgnc[c(1:(which(rowdups)[2]-1),(which(rowdups)[2]+1):nrow(exprots_hgnc)),]
}
cat(green("Clean genes: ",nrow(exprots_hgnc),"\n")) # 38158 transcripts



#####################################################
cat(red("3.0 Quality Control Pre Normalization"),"\n")
raw.counts <- data.expr[rownames(data.expr)%in%Annot$ensembl_gene_id,]
exp.conditions <- designExp[,"tumor_stage", drop = FALSE]
exp.conditions$tumor_stage = as.factor(exp.conditions$tumor_stage)
mydata <- NOISeq::readData(
  data = raw.counts, 
  length = Annot[,c("ensembl_gene_id","lenght")], 
  biotype = Annot[,c("ensembl_gene_id","gene_biotype")],
  chromosome = Annot[, c("chromosome_name", "start_position", "end_position")], 
  factors = exp.conditions,
  gc = Annot[, c("ensembl_gene_id", "percentage_gene_gc_content")]
)
#! Biodetection plot. Per group.
mybiodetection <-  NOISeq::dat(mydata, type="biodetection", factor="tumor_stage", k=0)
png(filename=paste(PREDIR, "01-biodetection.Rd_%05d.png", sep="/"),  width=w, height=h, pointsize=p)
explo.plot(mybiodetection, factor="tumor_stage" )
dev.off()
cat("Biodetection plots generated\n")
#! Count distribution per biotype. Using count per million, only for one sample
mycountsbio <-  NOISeq::dat(mydata, factor = NULL, type = "countsbio")
png(filename=paste(PREDIR, "02-countsbio.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")
dev.off()
cat("Counts distribution plot per biotype and one sample generated\n")
#! Count distribution per sample
png(paste(PREDIR, "02-protein_coding_boxplot.png", sep="/"), width=w*2, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "boxplot")
dev.off()
cat("Counts distribution plot for protein coding and all samples generated\n")
png(paste(PREDIR, "02-protein_coding_barplot.png", sep="/"), width=w*2, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "barplot")
dev.off()
cat("Counts distribution barplot for protein coding biotype and all samples generated\n")
#! Count distribution per Experimental factors
mycountsbio <- NOISeq::dat(mydata, factor = "tumor_stage", type = "countsbio")
png(paste(PREDIR, "03-protein_coding_boxplot_group.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "boxplot")
dev.off()
cat("Counts distribution boxplot for protein coding biotype and group generated\n")
png(paste(PREDIR, "04-protein_coding_barplot_group.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "barplot")
dev.off()
cat("Counts distribution barplot for protein coding biotype and group generated\n")
#! BIAS
#! Length bias detection
mylengthbias <- NOISeq::dat(mydata, factor="tumor_stage", type="lengthbias")
png(paste(PREDIR, "05-Lengthbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mylengthbias, samples = NULL, toplot = "global")
dev.off()
cat("Lenght bias plot generated\n")
#! GC bias
mygcbias <- NOISeq::dat(mydata, factor = "tumor_stage", type="GCbias")
png(paste(PREDIR, "06-GCbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mygcbias, samples = NULL, toplot = "global")
dev.off()
cat("GC bias plot generated\n")
#! PCA
pca.dat <- NOISeq::dat(mydata, type = "PCA", logtransf = F)
pca.results <- pca.dat@dat$result
## Variance explained by each component
pdf(file=paste(PREDIR, "08-PCAVariance_raw.pdf", sep="/"), width = 4*2, height = 4*2)
barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance")
dev.off()
cat("PCA variance raw plot generated.\n")
## Loading plot
pdf(file=paste(PREDIR, "09-PCALoading_raw.pdf", sep="/"), width = 4*2, height = 4*2)
plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA loadings",
     xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
     ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))
dev.off()
cat("PCA loading raw plot generated.\n")
#Los colores solo funcionan para hasta cinco grupos
#SET COLORS BY GROUP
colores <- c("gold", "red2","blue1","green2","cyan","azure4")
grupos <- unique(as.character(exp.conditions$tumor_stage))
colores <- colores[1:length(grupos)]
names(colores) = grupos
cat(green(colores) %+% "\n")
mycol <- as.character(exp.conditions$tumor_stage)
for(i in grupos[!is.na(grupos)]){
  mycol[exp.conditions$tumor_stage == i] <- unname(colores[i])
}
## Score plot
pdf(file=paste(PREDIR, "10-PCAScore_raw.pdf", sep="/"), width = 5*2, height = 5)
par(mfrow = c(1,2))
# PC1 & PC2
rango <- diff(range(pca.results$scores[,1:2]))
plot(pca.results$scores[,1:2], col = "white",
    xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
    ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
    main = "PCA scores",
    xlim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1),
    ylim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,2], col = mycol, cex = 1.5)
# legend("topright", c("Ctr", "s1-", "s2-", "s3-", "s4-"),
         # col = c("gold", "red2", "blue1", "green2", "cyan"), ncol = 2, pch = 1)
legend("topright", names(colores), col =unname(colores), ncol = 2, pch = 1)
# PC1 & PC3
rango2 = diff(range(pca.results$scores[,c(1,3)]))
plot(pca.results$scores[,c(1,3)], col = "white",
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
     main = "PCA scores",
     xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
     ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
legend("topright", names(colores), col =unname(colores), ncol = 2, pch = 1)
dev.off()
cat("PCA scores raw plot generated.\n")



######################################################
cat(red("4.0 Solve Biases or Normalization"),"\n")
#! Get only protein Coding
goodgenes <- rownames(exprots_hgnc)%in%Annot$ensembl_gene_id[Annot$gene_biotype=="protein_coding"]
exprots_hgnc <- exprots_hgnc[goodgenes,]
#! 1) filter low count genes with the threshold found in BIASES(2): The CPM are used because the size 
#! of the library per sample affects the power of the experiment. CPM=(counts/fragments
#! sequenced)*one million. Filtering those genes with average CPM below 1, would be different
#! to filtering by those with average counts below 1. 
#! Drop rows with CPM >10
countMatrixFiltered = exprots_hgnc[rowMeans(exprots_hgnc)>10,]
cat(green("Genes filtered with CPM>10: ",nrow(countMatrixFiltered),"\n"))
#! 16372 features

#! METHOD 2
  #!Preparing normalization
  genesFiltered = Annot$ensembl_gene_id%in%rownames(countMatrixFiltered)
  annotFiltered = Annot[genesFiltered,]
  rownames(annotFiltered) <- annotFiltered[,1]
  fact.etapas <- data.frame(tumor_stage=designExp$tumor_stage, 
    row.names=colnames(countMatrixFiltered))

	# Full, Full, TMM
	ln.data <- withinLaneNormalization(countMatrixFiltered, annotFiltered$lenght, which = "full")
	gcn.data <- withinLaneNormalization(ln.data, 
		annotFiltered$percentage_gene_gc_content, which = "full")
	norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
  noiseqData <- NOISeq::readData( data = norm.counts,  factors = fact.etapas)
	ARSyn = ARSyNseq(noiseqData, norm = "n", logtransf = FALSE)



######################################################
cat(red("5.0 Quality Control Post Normalization"),"\n")
normData <- NOISeq::readData(
  data = exprs(ARSyn), 
  length = annotFiltered[,c("ensembl_gene_id","lenght")], 
  factors = fact.etapas,
  gc = annotFiltered[, c("ensembl_gene_id", "percentage_gene_gc_content")]
)
#! Length bias detection
mylengthbias <- NOISeq::dat(normData, factor="tumor_stage", type="lengthbias")
png(paste(POSTDIR, "01-Lengthbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mylengthbias, samples = NULL, toplot = "global")
dev.off()
cat("Lenght bias plot generated\n")
#! GC bias
mygcbias <- NOISeq::dat(normData, factor = "tumor_stage", type="GCbias")
png(paste(POSTDIR, "02-GCbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mygcbias, samples = NULL, toplot = "global")
dev.off()
cat("GC bias plot generated\n")
#! PCA
pca.dat <- NOISeq::dat(normData, type = "PCA", logtransf = F)
pca.results <- pca.dat@dat$result
## Variance explained by each component
pdf(file=paste(POSTDIR, "03-PCAVariance_raw.pdf", sep="/"), width = 4*2, height = 4*2)
barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance")
dev.off()
cat("PCA variance raw plot generated.\n")
## Loading plot
pdf(file=paste(POSTDIR, "04-PCALoading_raw.pdf", sep="/"), width = 4*2, height = 4*2)
plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA loadings",
     xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
     ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))
dev.off()
cat("PCA loading raw plot generated.\n")
#Los colores solo funcionan para hasta cinco grupos
#SET COLORS BY GROUP
colores <- c("gold", "red2","blue1","green2","cyan","azure4")
grupos <- unique(as.character(exp.conditions$tumor_stage))
colores <- colores[1:length(grupos)]
names(colores) = grupos
cat(green(colores) %+% "\n")
mycol <- as.character(exp.conditions$tumor_stage)
for(i in grupos[!is.na(grupos)]){
  mycol[exp.conditions$tumor_stage == i] <- unname(colores[i])
}
## Score plot
pdf(file=paste(POSTDIR, "05-PCAScore_raw.pdf", sep="/"), width = 5*2, height = 5)
par(mfrow = c(1,2))
# PC1 & PC2
rango <- diff(range(pca.results$scores[,1:2]))
plot(pca.results$scores[,1:2], col = "white",
    xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
    ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
    main = "PCA scores",
    xlim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1),
    ylim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,2], col = mycol, cex = 1.5)
# legend("topright", c("Ctr", "s1-", "s2-", "s3-", "s4-"),
         # col = c("gold", "red2", "blue1", "green2", "cyan"), ncol = 2, pch = 1)
legend("topright", names(colores), col =unname(colores), ncol = 2, pch = 1)
# PC1 & PC3
rango2 = diff(range(pca.results$scores[,c(1,3)]))
plot(pca.results$scores[,c(1,3)], col = "white",
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
     main = "PCA scores",
     xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
     ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
legend("topright", names(colores), col =unname(colores), ncol = 2, pch = 1)
dev.off()
cat("PCA scores raw plot generated.\n")



######################################################
cat(red("6.0 Saving Data"),"\n")
setwd(paste0(raiz,OUTDIR))
for (i in levels(fact.etapas$tumor_stage)){
  cat("Saving grupo:",green(i),"\n")
  grupo <- ARSyn[,ARSyn$tumor_stage==i]
  cat("Samples: ",green(ncol(grupo)),"\n")
  datos <- exprs(grupo)
  d <- cbind(rownames(datos), data.frame(datos, row.names=NULL))
  colnames(d)[1] <- "genes"
  fname <- paste0(EXP.NAME, "-", gsub(" ","", i), ".tsv" )
  write.table(d, file=fname, sep='\t', quote=FALSE, row.names=FALSE)
}
ids <- data.frame()
for (i in levels(fact.etapas$tumor_stage)){
  cat("Saving grupo:",green(i),"\n")
  grupo <- ARSyn[,ARSyn$tumor_stage==i]
  grupo <- ARSyn[,ARSyn$tumor_stage==i]
  df <- data.frame(IDs = sampleNames(grupo), tipo = i)
  ids <- rbind(ids,df)
}
write.table(ids, file="ids.tsv", sep='\t', quote=FALSE, row.names=FALSE)




