library(UpSetR)

MIs <- list()

stages <- c("ctrl","stage1","stage2","stage3","stage4")
for (s in stages){
  fname <- paste0("data/tables/mi/",s,"-all-genmirna-100k.txt")
  mi <- read.table(file = fname, sep = "\t", header = TRUE)
  print(fname)
  MIs[[s]] <- paste(mi$Source,mi$Target)
}

colored <- list("stage1","stage2","stage3","stage4")
upset(fromList(MIs), order.by = "freq", sets.bar.color = "#56B4E9",
      # set.metadata = metadata,
      queries = list(list(query = intersects, params = colored, 
                          color = "orange", active = T)),
      # list(type = "matrix_rows", column = "types", 
      #      colors = c(NT = "green", T = "navy")),
      mainbar.y.label = "Number of intereactions (log10)",
      text.scale = c(3, 3, 3, 3, 3, 3),
      number.angles = 30,
      scale.intersections = "log10",
      show.numbers = "no",
      line.size = 3,
      point.size = 6
)


metadata <- as.data.frame(cbind(stages,c("NT", "T","T","T","T")))
names(metadata) <- c("sets", "types")
upset(fromList(MIs), 
      set.metadata = list(data = metadata, 
                          plots = list(list(type = "hist", column = "types", assign = 20))))




###### DRAFT
library(stringr)
bm <- read.table(file = "data/tables/biomart.csv", sep = ",", header = TRUE)
elements <- unique(unlist(MIs))
m <- fromList(MIs)

v <- as.vector(as.integer(rownames(m[(m$stage1 == 1) & (m$stage2 == 1)& (m$stage3 == 1) & (m$stage4 == 1) & (m$ctrl == 0),])))
df <- data.frame(dat = elements[v])
df[c('gene', 'miR')] <- str_split_fixed(df$dat, ' ', 2)
df1 <- merge(df,bm,by.x ="gene", by.y = "ensembl_gene_id")
df1$gene_name
write.table(df1, file = "data/tables/genes-only-stages.tsv", 
            row.names=FALSE, sep="\t")


v <- as.vector(as.integer(rownames(m[(m$stage1 == 1) & (m$stage2 == 1)& (m$stage3 == 1) & (m$stage4 == 1) & (m$ctrl == 1),])))
df <- data.frame(dat = elements[v])
df[c('gene', 'miR')] <- str_split_fixed(df$dat, ' ', 2)
df1 <- merge(df,bm,by.x ="gene", by.y = "ensembl_gene_id")

# fromList()
  # elements <- unique(unlist(MIs))
  # data <- unlist(lapply(MIs, function(x){x <- as.vector(match(elements, x))}))
  # data[is.na(data)] <- as.integer(0); data[data != 0] <- as.integer(1)
  # data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  # data <- data[which(rowSums(data) !=0), ]
  # names(data) <- names(input)
  # return(data)


# movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
#                    header = T, sep = ";")


# mi <- read.table(file = "data/tables/mi/stage4-all-genmirna-100k.txt",
#                  sep = "\t", header = TRUE)
#   mi <- read.table(file = "data/tables/mi/stage1-all-genmirna-100k.txt", 
#                    sep = "\t", header = TRUE)
# MIs$stage1 <- paste(mi$Source,mi$Target)