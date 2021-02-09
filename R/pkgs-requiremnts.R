#! sudo apt install r-base-core

# Select a directory to save library
home <- path.expand("~")
rutalib <- paste(home, "Rlibs", sep = "/")
dir.create(rutalib, showWarnings = TRUE, recursive = TRUE)
.libPaths(rutalib)


install.packages("openssl")

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BiocParallel")
BiocManager::install("NOISeq")
BiocManager::install("SummarizedExperiment")
BiocManager::install("GWASTools")
BiocManager::install("cqn")
BiocManager::install("Glimma")
BiocManager::install("DESeq")
BiocManager::install("EDASeq")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
#! BiocManager::install("TCGAbiolinks")

#! Error in curl::curl_fetch_memory(url, handle = handle) :
#!   Problem with the SSL CA cert (path? access rights?)
#! remove.packages('curl')
#! detach(package:curl)
#! library(curl)
#! quite R
#! sudo apt-get remove libcurl4-nss-dev
#! sudo apt-get install libcurl4-openssl-dev
devtools::install_github('BioinformaticsFMRP/TCGAbiolinks')
#! install.packages('curl')


#! Reinstall R ubuntu
#! To know were your package are installed. delete this folder and reinstall.
# .libPaths() 
#! sudo apt-get remove r-base-core
#! sudo apt-get remove r-base
#! sudo apt-get autoremove
#! sudo apt-get install r-base-core

install.packages("VennDiagram")
install.packages("crayon")