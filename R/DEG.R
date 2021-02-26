library(readr)
library(edgeR)
library(dplyr)
library(Glimma)

deg.process <- function(dir.express, conds, dir.out ) {
  
  print(conds)
  print(dir.express)
  print(dir.out)  
  all_matrices <- lapply(conds, function(cond){
    fname <- paste0(dir.express,"/expr-miRNA-", cond, ".tsv")
    print(fname)
    matrix <- read_tsv(fname)
    return(list(matrix = matrix[,-1], nrows = nrow(matrix), 
                ncols = ncol(matrix[,-1]), genes = matrix[,1], cond = cond))
  })

  # print(all_matrices)
  matrices <- lapply(all_matrices, "[[", 1)
  nrows <- unlist(lapply(all_matrices, "[[", 2))
  ncols <- unlist(lapply(all_matrices, "[[", 3))
  genes <- unlist(all_matrices[[1]][4])
  names(genes) <- NULL
  conds <- unlist(lapply(all_matrices, "[[", 5))

  M <- cpm(data.matrix(bind_cols(matrices)), log = T)
  rownames(M) <-genes
  targets <- data.frame(id = colnames(M), 
                        group = factor(unlist(mapply(rep, conds, ncols))))

  ## MUST BE TRUE ##
  stopifnot(all(nrows == length(genes)))

  ## DISEÑO ##
  dmatrix <- model.matrix(~1 + targets$group , data = targets)
   
  ## AJUSTE ##
  fit <- lmFit(M, dmatrix)
  head(fit$coefficients)

  fit <- eBayes(fit)

  fit$fdr <- apply(fit$p.value, 2, p.adjust, method="fdr")
  p <- 1-fit$fdr
  fit$B <- log(p/(1-p))

  topTable(fit, coef=ncol(dmatrix), adjust="BH")

  # Glimma plot
  dt <- decideTests(fit, lfc = 2)
  summary(dt)

  # glMDPlot(path = "Results/KIRC/DEG/", fit, counts = M, groups = targets$group, status = dt)
  glXYPlot(x=fit$coefficients[,2], y = fit$lods[,2], path =  paste0(dir.out,"/"), counts = M,
           xlab="logFC", ylab="B", groups = targets$group, status=dt[,2],
           folder = paste0("volcan-", conds[2]), launch = FALSE )

  genesFULL <- bind_cols(
    ensemblID = rownames(fit$coefficients),
    coef = fit$coefficients[, 2],
    p_value = fit$p.value[, 2],
    FDR = fit$fdr[, 2],
    B = fit$B[, 2]
  )

  write_tsv(genesFULL, paste0(dir.out,"/", conds[2], "-deg-ebayes.tsv"))
  
}

# DIRIN <- "Results/Expression"
DIRIN <- "Output"
DIROUT <- "Results/DEG"
dir.create(DIROUT, recursive = TRUE)

#! MAIN
conditions <- c("ctrl", "stagei")
deg.process( DIRIN, conditions, DIROUT)

conditions <- c("ctrl", "stageii")
deg.process( DIRIN, conditions, DIROUT)

conditions <- c("ctrl", "stageiii")
deg.process( DIRIN, conditions, DIROUT)

conditions <- c("ctrl", "stageiv")
deg.process( DIRIN, conditions, DIROUT)





#####################################################################
# home <- path.expand("~")
# rutalib <- paste(home, "Rlibs", sep = "/")
# dir.create(rutalib, showWarnings = TRUE, recursive = TRUE)
# .libPaths(rutalib)


