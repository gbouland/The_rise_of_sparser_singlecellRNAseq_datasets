library(scRNAseq)
library(magrittr)

calls <- c("AztekinTailData","WuKidneyData","MarquesBrainData","GrunHSCData")

cellRsquareds <- lapply(calls,function(call){
  message(call)
  dataset <- eval(call(call))
  expr <- dataset@assays@data$counts
  rm(dataset)
  expr <- t(t(expr)/colSums(expr))*10000
  expr <- log(expr + 1)
  binary <- (expr >=1)*1
  binary <- as.matrix(binary)
  expr <- as.matrix(expr)
  cellRsquared <- lapply(1:ncol(binary),function(i){
    cor(binary[,i],expr[,i])
  }) %>% unlist()
  rm(binary)
  rm(expr)
  return(cellRsquared)
})
names(cellRsquareds) <- calls
cellsRsquared <- stack(cellRsquareds)


pct0s <- lapply(calls,function(call){
  message(call)
  dataset <- eval(call(call))
  expr <- dataset@assays@data$counts
  rm(dataset)
  binary <- (expr >=1)*1
  binary <- as.matrix(binary)
  pct0 <- lapply(1:ncol(binary),function(i){
    mean(binary[,i])
  }) %>% unlist()
  rm(binary)
  return(pct0)
})
names(pct0s) <- calls
pct0s <- stack(pct0s)
cellsRsquared$pct0 <- pct0s$values


