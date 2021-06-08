require(biomaRt)
require(dplyr)
require(tidyverse)

dir.create("pipeline")

# mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# datasets = listDatasets(mart)
ensembl = useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl")
# head(searchAttributes(ensembl, pattern="(?i)mirbase"))
print("Downloading BIOMART data")
# check https://www.ensembl.org/info/website/archives/index.html
# if Error in getNodeSet subscript out of bounds
bm <- getBM(attributes=c('ensembl_gene_id',
                         'description',
                         'gene_biotype',
                         'percentage_gene_gc_content',
                         'start_position',
                         'end_position',
                         'transcript_length',
                         'chromosome_name',
                         'band',
                         'external_gene_name',
                         'mirbase_id'),
            mart = ensembl,
            verbose = FALSE)
chroms <- c(1:22,"X","Y")
dtmp <- subset(bm, chromosome_name %in% chroms)
dtmp <- subset(dtmp, gene_biotype %in% 
                 c("protein_coding","miRNA"))
dtmp <- dtmp[!duplicated(dtmp$ensembl_gene_id), ]
dtmp$geneLength <- dtmp$end_position - dtmp$start_position
dtmp <- rename(dtmp, gcContent = percentage_gene_gc_content)
dtmp <- rename(dtmp, chr = chromosome_name)
dtmp <- rename(dtmp, gene_name = external_gene_name)
dtmp <- dtmp  %>% drop_na(gene_name)
dtmp <- dtmp[dtmp$gene_name != "", ]
# 21196 lines in biomart-20210517.csv
write.csv(dtmp, file = "pipeline/biomart.csv", row.names = FALSE)
