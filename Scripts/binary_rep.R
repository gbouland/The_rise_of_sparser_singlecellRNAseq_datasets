library(bit)
library(Matrix)
library(magrittr)
library(Seurat)

data <- readRDS("/tudelft.net/staff-umbrella/pQTL/000scBinary/002 Preprocessed/mouseDroplet.rds")
counts <- as.matrix(data[[1]])
seur <- CreateSeuratObject(counts = counts)
seur <- SCTransform(seur)

binary <- (counts >=1)
List <- vector("list", ncol(binary))
for(i in 1:length(List)){
  message(i)
  newBit <- bit(nrow(binary))
  newBit[which(binary[,i])] <- T
  List[[i]] <- newBit
}
bitStored <- format(object.size(List), units = "Gb")%>% gsub(pattern = " Gb",replacement = "") %>% as.numeric()
scTransform <- format(object.size(seur@assays$SCT), units = "Gb")%>% gsub(pattern = " Gb",replacement = "") %>% as.numeric()
