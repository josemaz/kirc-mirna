library(UpSetR)

MIs <- list()

for (s in c("ctrl","stage1","stage2","stage3","stage4")){
  fname <- paste0("data/tables/mi/",s,"-all-genmirna-100k.txt")
  mi <- read.table(file = fname, sep = "\t", header = TRUE)
  print(fname)
  MIs[[s]] <- paste(mi$Source,mi$Target)
}

colored <- list("stage1","stage2","stage3","stage4")
upset(fromList(MIs), order.by = "freq", sets.bar.color = "#56B4E9",
      queries = list(list(query = intersects, params = colored, 
                          color = "orange", active = T)),
      mainbar.y.label = "Intereactions"
)




###### DRAFT

mi <- read.table(file = "data/tables/mi/ctrl-all-genmirna-100k.txt", 
                 sep = "\t", header = TRUE)
MIs$ctrl <- 
  
  mi <- read.table(file = "data/tables/mi/stage1-all-genmirna-100k.txt", 
                   sep = "\t", header = TRUE)
MIs$stage1 <- paste(mi$Source,mi$Target)