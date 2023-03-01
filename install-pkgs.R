r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

cat("Using LiBRARY: ",.libPaths(),"\n")

install.packages("crayon")
require("crayon")

# System Packages
packs <- c("UpSetR", "BiocManager", "DT", "tidyverse", "devtools", 
           "VennDiagram", "ggfortify", "gprofiler2")
noinst <- setdiff(packs, rownames(installed.packages()))
cat(red(noinst) %+% '\n')
install.packages(noinst)
# Bioconductor Packages
packs <- c("TCGAbiolinks", "EDASeq", "edgeR", "EnhancedVolcano",
           "sva", "NOISeq","biomaRt", "DESeq2")
noinst <- setdiff(packs, rownames(installed.packages()))
cat(red(noinst) %+% '\n')
BiocManager::install(noinst)
# From GITHUB
packs <- data.frame( p =  c("vqv/ggbiplot") )
packs$p2 <-  sapply(strsplit(packs$p, '/'), "[[", 2)
packs <- merge(packs,installed.packages(),by.x="p2",by.y=0,all.x=TRUE)
packs <- packs[is.na(packs$Package),]
cat(red(packs$p) %+% '\n')
devtools::install_github(packs$p)












# #! BiocManager::install(version = "3.13")
# BiocManager::install("TCGAbiolinks")
# BiocManager::install("EDASeq")
# BiocManager::install("edgeR")
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("sva")
# BiocManager::install("NOISeq")
# BiocManager::install("biomaRt")
# #! BiocManager::install("biomaRt", version = "2.46")


# install.packages('DT')
# install.packages('tidyverseu')
# install.packages("devtools")
# ## devtools::install_github("yanlinlin82/ggvenn")
# install.packages("VennDiagram")
# install.packages("ggfortify")
# install.packages("ggbiplot")
# # devtools::install_github("vqv/ggbiplot")
