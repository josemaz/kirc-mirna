# source("R/02-miRNAs.R")

library(dplyr)
library(tidyr)
# library(plyr)
library(crayon)
library(tibble)
library(TCGAbiolinks)
library(stringr)
library(Biostrings)

raiz <- getwd()

# Moverme a TCGA directory
setwd("pipeline/TCGA") 

#######################################################################
cat(red("0.5 Download Data"),"\n")
query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Isoform Expression Quantification",
                  workflow.type = "BCGSC miRNA Profiling")
GDCdownload(query)
#! gene annotation GRCh38.p13
data <- GDCprepare(query)
save(data, file="miRNAs-isoform.RDATA", compress="xz")
load(file="miRNAs-isoform.RDATA")


#########################################################################
cat(red("1.0 Processing miRNA Data"),"\n")
#! read_counts tiene 50% de ceros y su promedio es menor a 10
dfwide <- pivot_wider(data, names_from = barcode, values_from  = read_count) 
dfwide <- select(dfwide, -c(1,2,3,4))
dfwide[is.na(dfwide)] <- 0
df <- dfwide %>% 
  group_by(miRNA_region) %>% 
  summarise(across(everything(), sum)) # por decirdir suma o promedio
#! fix rownames and data.frame
df <- filter(df,str_detect(miRNA_region, "MIMAT"))
# out <- strsplit(as.character(before$type),'_and_') 
# do.call(rbind, out)
df <- separate(df, miRNA_region, c("X", "mimats"), ",")
df <- df[,-1]
mirnas.expr <- as.data.frame(df[,-1])
rownames(mirnas.expr)  <- df$mimats
#! Number of Zeros
# apply(mirna.expr, 1 , function(x) length(which(x == 0)))
#! Zeros filter
dim(mirnas.expr)
#! [1] 2102  617
#! df <- df[apply(df!=0, 1, all),] # por diecidir
threshold <- round(dim(mirnas.expr)[2]/2)
filtro <- mirnas.expr[rowSums(mirnas.expr == 0) <= threshold, ]
dim(filtro)
#! [1] 1509  616
mirnas.expr <- filtro[rowMeans(filtro) >= 2, ]
rm(df,dfwide,filtro)

#! Checar pacientes unicos


#########################################################################
cat(red("2.0 Processing miRBase Data"),"\n")
setwd(raiz)
s = readDNAStringSet("Data/mature-miRBase-20210114.fa")

setwd(paste0(raiz,"/Output"))
seqs.names <- names(s)
tabla <- do.call(rbind, str_split(seqs.names," "))
MIMATs <- data.frame(tabla[,c(1,2)])
colnames(MIMATs) <- c("mature","precursor")
rm(s,tabla,seqs.names)

df <- data.frame(MIMATs = row.names(mirnas.expr), mirnas.expr)
joined <- right_join(MIMATs, df,
                           by = c("precursor" = "MIMATs"))
joined<-joined[,-2]
joined <- drop_na(joined, mature)
# write.table(joined, file="miRNA-expr.tsv", sep='\t', 
# 	quote=FALSE, row.names = FALSE)
rownames(joined) <- joined$mature
joined<-joined[,-1]
colnames(joined) <- str_replace_all(colnames(joined),"\\.","-")



#########################################################################
cat(red("3.0 Loading RNAseq Data"),"\n")
mirna <- data.frame(IDs = unlist(colnames(joined), use.names = FALSE))
mirna$IDs <- str_replace_all(mirna$IDs,"\\.","-")
mirna$cases <- str_sub(mirna$IDs,1,20)
mirna[duplicated(mirna[,'cases']),]
dim(mirna) # [1] 616   2

# rnaseq <- read.table(file="ids.tsv", sep='\t', header=TRUE)
rnaseq <- read.csv("ids.tsv", sep = "\t")
rnaseq$IDs = as.character(rnaseq$IDs)
rnaseq$cases <- str_sub(rnaseq$IDs,1,20)
rnaseq$tipo <- str_replace_all(rnaseq$tipo," ","")
dim(rnaseq) # [1] 608   3

#########################################################################
cat(red("4.0 Matching with RNAseq Data"),"\n")
rnaseq <- rnaseq[!duplicated(rnaseq[c(2,3)]),]
cat(green("rnaseq unicos:",dim(rnaseq)[1]),"\n") # [1] 568   4
merged <- inner_join(mirna, rnaseq, by = "cases")
merged$tipo <- as.factor(merged$tipo)
colnames(merged)[1] <- "mirna"
colnames(merged)[3] <- "rnaseq"
cat(green("InerJoin mirnas-rnaseq:",dim(merged)[1]),"\n")
merged <- merged[!duplicated(merged$rnaseq),]
cat(green("mirnas-rnaseq unicos:",dim(merged)[1]),"\n")
print(merged[duplicated(merged$rnaseq),])
print(merged[duplicated(merged$mirna),])

# for( i in rownames(rnaseq) ){
# 	rnaseq.case = rnaseq[i, "cases"]
# 	res <- dim(mirna[mirna$cases==rnaseq.case,])[1]
# 	if(res==1){
# 		print("miRNA unico")
# 	} else if(res ==0) {
# 		print("Sin mirna")
# 	} else if(res > 1) {
# 		print("mirna doble o mas")
# 	}
# }


#########################################################################
cat(red("5.0 Saving Data"),"\n")
targets <- data.frame()
for (i in levels(merged$tipo)){
	print(i)
	grupo <- merged[merged$tipo==i,]
	grupo$seq <- 1:nrow(grupo)
	grupo$label <- paste0(grupo$tipo,"-",grupo$seq)
	grupo <- grupo[!names(grupo) %in% "seq"]
	targets <- rbind(targets,grupo)
}
write.table(targets, file="targets.tsv", sep='\t', quote=FALSE, row.names = FALSE)

for (i in levels(targets$tipo)){
	print(i)
	fname <- paste0("kirc-",i,".tsv")
	expr <- read.csv(paste0(fname), sep = "\t")
	rownames(expr) <- expr[,1]
	expr <- expr[,-1]
	colnames(expr) <- str_replace_all(colnames(expr),"\\.","-")
	grupo <- targets[targets$tipo==i,]
	cat(green("cases ",i,": ",length(colnames(expr))),"\n")
	cat(green("cases miRNA: ",dim(grupo)[1]),"\n")
	rnaseq.all <- expr[grupo$rnaseq]
	cat(green("RNAseq cases: ",dim(rnaseq.all)[2]),"\n")
	mirnas.all <- joined[grupo$mirna]
	cat(green("miRNA cases: ",dim(mirnas.all)[2]),"\n")
	colnames(rnaseq.all) <- grupo$label
	colnames(mirnas.all) <- grupo$label
	all <- rbind(rnaseq.all,mirnas.all)
	all$gene <- rownames(all)
	all <- all %>% select(gene, everything())
	fname <- paste0("paste-miRNA-",i,".tsv")
	write.table(all, file=fname, sep='\t', quote=FALSE, row.names = FALSE)
}
















# df <- data.frame(data$reads_per_million_miRNA_mapped,data$miRNA_region)
# dfr <- df %>% 
#   subset(data.reads_per_million_miRNA_mapped >= 1)  %>% 
#   filter(str_detect(data.miRNA_region, "MIMAT")) 
# ny <- ddply(dfr,.(data.miRNA_region),summarize,sum_rpm=sum(data.reads_per_million_miRNA_mapped))
# ny$sum_rpm <- log2(ny$sum_rpm)
# tmp <- ny %>%
#   separate(data.miRNA_region, c("X", "prec"), ",") 





# library(miRBaseVersions.db)
# keytypes(miRBaseVersions.db)
# 
# result = select(miRBaseVersions.db, 
#                 keys = "MIMAT0000092", 
#                 keytype = "MIMAT", 
#                 columns = "*")


# png("rplot.png", width = 350, height = 350)
# # hist(rowMeans(mirnas.expr),breaks = 1000)
# plot(density(rowMeans(mirnas.expr)), log="y")
# dev.off()

# lapply(data, function(x) { gsub("< ", "<", x) })
