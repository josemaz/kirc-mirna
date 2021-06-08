require(TCGAbiolinks)
require(SummarizedExperiment)
require(NOISeq)
require(EDASeq)
require(DESeq2)
require(EnhancedVolcano)
require(dplyr)
require(crayon)
require(stringr)
require(VennDiagram)
require(gprofiler2)

###################################################################3
# FILTER
#! Filter by minimum mean and number of zeros
filtro.zeros.means <- function(m1) {
  print("Input (dim): ")
  print(dim(m1))
  threshold <- round(dim(m1)[2]/2)
  print(paste0("Threshold: ",threshold))
  m1 <- m1[rowSums(m1 == 0) <= threshold, ]
  print(paste0("Rows After Zeros (dim): ",dim(m1)[1]))
  m1 <- m1[rowMeans(m1) >= 10, ]
  print(paste0("Rows After means (dim): ",dim(m1)[1]))
  return(m1)
}

####################################################################
#!-- PCA plotting by group
pca.grupo <- function(dat, fout = NULL){
  d2 <- data.frame(etapa = colData(dat)$grupo)
  # d2$etapa <-as.factor(d2$etapa)
  rownames(d2) <- colnames(dat)
  mydat = NOISeq::readData( assay(dat) , factors = d2)
  # myPCA = dat(mydat, type = "PCA", norm=TRUE)
  myPCA = dat(mydat, type = "PCA", logtransf = F)
  if(!is.null(fout)){
    print(paste0("Writing in: ",fout))
    png(fout)
  }
  explo.plot(myPCA, factor = "etapa", plottype = "scores")
  if(!is.null(fout))dev.off()
}

####################################################################
#!-- miRNA Samples Cleaning
clean.mirna <- function(miRNA.exp, dat.bm){
  # CLEAN and ANNOT
  m <- select(miRNA.exp,starts_with("read_count_"))
  rownames(m) <- miRNA.exp$miRNA_ID
  m <- filtro.zeros.means(m) # esta en tools.R
  colnames(m) <- gsub('read_count_', '', colnames(m), fixed=TRUE)
  #! Prepare ColData  
  colData <- DataFrame( row.names = colnames(m),
                        barcode = colnames(m) )
  colData$sample <- str_extract(colnames(m), 
                                "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{3}")
  #! Annotation
  # bm <- read.csv("pipeline/biomart.csv")
  print(paste0("[BIOMART] all rows: ",nrow(bm)))
  bm <- dat.bm[dat.bm$gene_biotype == "miRNA",]
  print(paste0("[BIOMART] only miRNA : ",nrow(bm)))
  #! Merge by annotation
  tmp <- merge(m,bm,by.x=0,by.y="mirbase_id")
  tmp <- tmp[!duplicated(tmp$Row.names),]
  nsamples = ncol(m)
  m <- tmp[,c(2:(nsamples+1))]
  rownames(m) <- tmp$Row.names
  bm <- tmp[,-c(1:(nsamples+1))]
  rownames(bm) <- tmp$Row.names
  
  #! Create SE only miRNA (275,616)
  mirna <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(m)),
                                colData=colData, rowData=bm)
  
  #! Post- clean
  mirna <- mirna[,!duplicated(mirna$sample)]
  print(paste0("[mirna] after duplicated samples: ",ncol(mirna)))
  return(mirna)
}

####################################################################
#!-- RNA Samples Cleaning
# TODO: Cambiar variable rnaseq interna por rnas
clean.rna <- function(RNA.exp,dat.bm){
  # TODO: print correctly dims reports
  #!-- Filtros
  print(dim(RNA.exp))
  lstg <- c("Primary solid Tumor","Solid Tissue Normal")
  rnaseq <- RNA.exp[,RNA.exp$definition %in% lstg]
  lstg <- c("Control","stage i","stage ii","stage iii","stage iv")
  # TODO: Rename "stage i" -> stage1, etc.
  rnaseq <- rnaseq[,rnaseq$tumor_stage %in% lstg]
  print(dim(rnaseq))
  
  #!-- Mutate and Changes
  clinic <- colData(rnaseq)
  rnaseq$grupo <- rnaseq$tumor_stage
  rnaseq$grupo[clinic$sample_type_id=="11"] <- "ctrl"
  rnaseq$grupo[rnaseq$grupo=="stage i"] <- "stage1"
  rnaseq$grupo[rnaseq$grupo=="stage ii"] <- "stage2"
  rnaseq$grupo[rnaseq$grupo=="stage iii"] <- "stage3"
  rnaseq$grupo[rnaseq$grupo=="stage iv"] <- "stage4"
  rnaseq$grupo <- as.factor(rnaseq$grupo)
  
  
  #--! Eliminate duplicates by group
  #! TODO: Cambiar por un codigo iterativo
  print(dim(rnaseq))
  print(paste0("Sample dups: ",ncol(rnaseq[,duplicated(rnaseq$sample)])))
  clinic <- colData(rnaseq)
  levels(clinic$grupo)
  d1 <- clinic[clinic$definition == "Primary solid Tumor",]
  d1 <- d1[!duplicated(d1$patient),] # 527 samples
  d2 <- clinic[clinic$definition == "Solid Tissue Normal",]
  d2 <- d2[!duplicated(d2$patient),] # 72 samples
  rnaseq <- rnaseq[,c(rownames(d1),rownames(d2))]
  rm(d1,d2,lstg,clinic)
  print(dim(rnaseq))
  
  #!-- GENES Cleaning
  m <- filtro.zeros.means(assay(rnaseq))
  rnaseq <- rnaseq[rownames(m),]
  dim(rnaseq)
  
  #!-- Annotation
  # bm <- read.csv("pipeline/biomart.csv")
  bm <- dat.bm[dat.bm$gene_biotype == "protein_coding",]
  bm <- bm[!duplicated(bm$ensembl_gene_id),]
  rnaseq <- rnaseq[rownames(rnaseq) %in% bm$ensembl_gene_id, ]
  print(dim(rnaseq))
  bm <- bm[bm$ensembl_gene_id %in% rownames(rnaseq),]
  rowData(rnaseq) <- merge(rowData(rnaseq),bm,by="ensembl_gene_id")
  print(dim(rnaseq))
  return(rnaseq)
}

###################################################################3
# Normalization and Bias correct of RNAseq
norm.rna <- function(dat.rna){
  pca.grupo(dat.rna,fout = "pipeline/PCA-rna-BeforeNorm.png")
  fac <- data.frame(tumor_stage=dat.rna$grupo, 
                    row.names=colnames(dat.rna))
  #! Pre Normalization
  ln.data <- withinLaneNormalization(assay(dat.rna), 
                                     rowData(dat.rna)$geneLength, which = "full")
  gcn.data <- withinLaneNormalization(ln.data , rowData(dat.rna)$gcContent,
                                      which = "full")
  # norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
  norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
  noiseqData <- NOISeq::readData( norm.counts, factors = fac)
  mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
  #! Post Normalization
  assay(dat.rna) <- exprs(mydata2corr1)
  pca.grupo(dat.rna,fout = "pipeline/PCA-rna-AfterNorm.png")
  return(dat.rna)
}

###################################################################3
#!- Normalization and Bias correct of miRNA
norm.mir <- function(dat.mir){
  pca.grupo(dat.mir,fout = "pipeline/PCA-mir-BeforeNorm.png")
  fac <- data.frame(tumor_stage=dat.mir$grupo, 
                    row.names=colnames(dat.mir))
  #! Pre Normalization
  noiseqData <- NOISeq::readData( assay(dat.mir), factors = fac)
  mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
  #! Post Normalization
  assay(dat.mir) <- exprs(mydata2corr1)
  pca.grupo(dat.mir,fout = "pipeline/PCA-mir-AfterNorm.png")
  return(dat.mir)
}

###################################################################3
#!- Download paired
download.paired <- function(loading = F){
  dir.create("pipeline")
  setwd("pipeline")
  # Si no existen los dos objetos tambien los descarga
  fs <- !file.exists("RNAseq-tcga.rda") & !file.exists("miRNAs-tcga.rda")
  if(!loading | fs){
    #! 616 cases
    query.mir <- GDCquery(project = "TCGA-KIRC",
                      data.category = "Transcriptome Profiling",
                      data.type = "miRNA Expression Quantification",
                      workflow.type = "BCGSC miRNA Profiling")
    #! 611 cases
    query.rna <- GDCquery(project = "TCGA-KIRC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")
    #! common 512
    common.patients <- intersect(substr(getResults(query.mir, cols = "cases"), 1, 12),
                                 substr(getResults(query.rna, cols = "cases"), 1, 12))
    #! download
    query <- GDCquery(project = "TCGA-KIRC",
                      data.category = "Transcriptome Profiling",
                      data.type = "miRNA Expression Quantification",
                      workflow.type = "BCGSC miRNA Profiling",
                      barcode = common.patients)
    GDCdownload(query)
    mirs <- GDCprepare(query, summarizedExperiment = F)
    query <- GDCquery(project = "TCGA-KIRC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts",
                      barcode = common.patients)
    GDCdownload(query)
    rnas <- GDCprepare(query = query, summarizedExperiment = TRUE)
    save(mirs, file="miRNAs-tcga.rda")
    save(rnas, file="RNAseq-tcga.rda")
    rm(query,query.mir,query.rna)
  }else{
    load(file="miRNAs-tcga.rda")
    load(file="RNAseq-tcga.rda")
  }
  setwd("..")
  l <- list(mirs,rnas)
  names(l) <- c("mir", "rna")
  return(l)
}

####################################################################3
#!- Calculate all results of DEGs multigroup (todos vs. todos)
degs <- function(dats, form = ~ grupo){
    #! DEGs multigroup
    #! Ajuste de datos
    assay(dats) <- round(assay(dats)) # solo enteros
    print("[DEG] Genes duplicados en el gname: ")
    print(rownames(dats[duplicated(rownames(dats)),]))
    dats <- dats[!duplicated(rownames(dats)),]
    #! DEseq2
    ddsSE <- DESeqDataSet(dats, design = ~ grupo)
    dds <- DESeq(ddsSE)
    #! Combiantoria de los contrastes
    combi <- combn(levels(dats$grupo), 2)
    lres <- apply(combi,2, 
      function(x,dds1 = dds){
        print(green(c(x[2],x[1])))
        res <- results(dds1,contrast=c("grupo",x[2],x[1]))
        return(res)
      } 
    )
    names(lres) <- apply(combi,2, function(x) {return(paste0(x[2],"_",x[1]))} )
    # return(list(lres = lres, ups = ups, downs = downs, todos = todos))
    return(lres)
}

####################################################################3
#!- anota y escirbe lista de resultados de una lista de DEGs
deg.annot.write <- function(ldeg, bm, prefix, c){
  dir.create("output/DEGS", recursive = T)
  nres <- list()
  for (r.name in names(ldeg)){
    x1 <- unlist(strsplit(r.name,"_"))[1]
    x2 <- unlist(strsplit(r.name,"_"))[2]
    fout = paste0("output/DEGS/res-",prefix,"-",x1,"_",x2,".tsv")
    print(fout)
    d1 <- as.data.frame(ldeg[[r.name]])
    d2 <- bm[!duplicated(bm[,c]),]
    d3 <- merge(d1,d2,by.x=0, by.y=c, no.dups=F)
    # print(paste0(nrow(d3),"-",nrow(d1)))
    names(d3)[1] <- c
    stopifnot(nrow(d3)==nrow(d1))
    write.table(d3, file = fout, quote=FALSE, sep='\t', row.names = F)
    nres[[r.name]] <- d3
  }
  return(nres)
}

####################################################################3
#!- Paste DEG of diferent experiment
paste.degs <- function(ldeg1, ldeg2, prefix){
  stopifnot(length(ldeg1) == length(ldeg2))
  l.all <- list()
  for (i in names(ldeg1)){
    x1 <- unlist(strsplit(i,"_"))[1]
    x2 <- unlist(strsplit(i,"_"))[2]
    l.all[[i]] <- ldeg1[[i]] <- rbind(ldeg1[[i]],ldeg2[[i]])
    fout = paste0("output/DEGS/res-",prefix,"-",x1,"_",x2,".tsv")
    print(fout)
    write.table(l.all[[i]], file = fout, quote=FALSE, sep='\t', row.names = F)
  } 
  return(l.all)
}


####################################################################3
# #!- Cortes de LogFC y pavalue para los results
ups.downs <- function(lres, lfc = 1.0, pval = 1e-3){
  ups <- list()
  downs <- list()
  todos <- list()
  for (nom in names(lres)) {
    print(nom)
    r <- lres[[nom]]
    downs[[nom]] <- rownames(r[ (r$pvalue < pval) & (r$log2FoldChange < -lfc) ,])
    ups[[nom]] <- rownames(r[ (r$pvalue < pval) & (r$log2FoldChange > lfc) ,])
    todos[[nom]] <- union(downs[[nom]],ups[[nom]])
  }
  return(list(ups = ups, downs = downs, todos = todos))
}

####################################################################3
# #!- Plot volvano plots of list of DEGs results
volcanos <- function(lres, prefix="rna", lfc=1.0, pv=1e-3, c=NULL){
  dir.create("output/Volcanos", recursive = T)
  for (r in names(lres)) {
    fout <- paste0("output/Volcanos/vol-",
                   prefix,"-",r,".png")
    print(fout)
    # png(fout, width=800, height=600)
    print(EnhancedVolcano(lres[[r]],
                          lab = lres[[r]][,c],
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          pCutoff = pv,
                          FCcutoff = lfc,
                          subtitle = paste0(r,"-",lfc,"-",pv),
                          labSize = 5.5
    ))
    # ggsave(fout,dpi=300,device ="png", width=800, height=600, units="px")
    ggsave(fout,dpi=300,device ="png")
    dev.off()
  }
}

####################################################################3
#! Venn
venn.degs <- function(l=NULL, sets= NULL, prefix = NULL  ){
  stopifnot( (!is.null(l)) & (!is.null(prefix)) )
  print("Writing up-regultaed")
  venn.diag(l$ups[sets],paste0(prefix,"-ups"))
  print("Writing down-regultaed")
  venn.diag(l$down[sets],paste0(prefix,"-down"))
  print("Writing all-regultaed")
  venn.diag(l$todos[sets],paste0(prefix,"-all"))
}
venn.diag <- function(lsets,pfx){
  colors <- c("#6b7fff", "#c3db0f", "#7FC97F", "#BEAED4", "#FDC086")
  colors <- colors[1:length(lsets)]
  fpng <- paste0("output/venn-",pfx,".png")
  print(fpng)
  venn.diagram(x = lsets,
               # category.names = names(s),
               filename = fpng,
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
               main = pfx
  )
}

####################################################################3
#! Carga las redes de ARACNe desde una ruta
load.nets <- function(folder = "remote/MI", pat="*-1e5-genmirna.tsv"){
  redes <- list()
  for (i in list.files(path=folder, pattern=pat)){
    n <- unlist(str_split(i,'-'))[1]
    print(n)
    redes[[n]] <- read.table(
      file = paste0(folder,"/",i),
      header = T, sep = '\t')
  }
  return(redes)
}

####################################################################3
#! CEnriquece los genes targes de una lista de genes y mirna regulados
enrich.target.genes <- function(net, lmir, lrna, bm){
  filtro1 <- net[(net$Source %in% lmir) | (net$Target %in% lmir),]
  filtro2 <- filtro1[(filtro1$Source %in% lrna) | (filtro1$Target %in% lrna),]
  l <- c(
    filtro2$Source[grepl("^ENSG", filtro2$Source)],
    filtro2$Target[grepl("^ENSG", filtro2$Target)]
  )
  l <- unique(l)
  print(length(l))
  if(length(l)<30){
    print(paste0(as.character(lmir),collapse=", "))
    d <- data.frame(ensid = l)
    d <- merge(d,bm,by.x="ensid",by.y="ensembl_gene_id")
    print(paste0(as.character(d$gene_name),collapse=", "))
    # print(paste(shQuote(d$), collapse=", "))
  }
  gostres <- gost(query = l, organism = "hsapiens")
  
  print("++++++++++++++++++++++++++++")
  res <- gostres$result
  res <- res[order(res$p_value),]
  res <- res[res$source %in% c("GO:CC","GO:BP","GO:MF"),]
  print(head(res))
}










# seq.grouping <- function(d){
#   d <- d[d$definition %in% c("Primary solid Tumor", "Solid Tissue Normal"),]
#   d = d %>% mutate(grp = case_when(
#     definition == "Solid Tissue Normal" ~ "NT",
#     tumor_stage == "stage i" ~ "st_1",
#     tumor_stage == "stage ii" ~ "st_2",
#     tumor_stage == "stage iii" ~ "st_3",
#     tumor_stage == "stage iv" ~ "st_4",)
#   )
#   dim(d)
#   d <- d[!is.na(d$grp),]
#   dim(d)
#   d[duplicated(d$sample),]
#   dim(d)
#   d1 <- d[d$definition == "Primary solid Tumor",]
#   d1 <- d1[!duplicated(d1$patient),]
#   d2 <- d[d$definition == "Solid Tissue Normal",]
#   d2 <- d2[!duplicated(d2$patient),]
#   d3 <- rbind(d1,d2)
#   return(d3)
# }
