library(scRNAseq)
library(Matrix)
library(magrittr)
library(tidyr)

datasets <- as.data.frame(listDatasets())
datasets$year <- tidyr::extract_numeric(datasets$Reference)
datasets <- datasets[1:2,]

pcts <- rep(0,nrow(datasets))

for(i in 1:length(pcts)){
  call <- datasets[i,"Call"]
  main <- strsplit(x = call,split = "[(]") %>% sapply(FUN = function(x)x[1])
  Arg <- strsplit(x = call,split = "[(]") %>% sapply(FUN = function(x)x[2])
  Arg <- gsub(pattern = ")","",Arg)
  if(Arg == ""){
    message(main)
    fun <- get(main)
    data <- fun()
    if(names(data@assays@data) == "counts"){
      expr <- data@assays@data$counts>=1
      rm(data)
      pct0 <- sum(colSums(expr)) / (as.numeric(nrow(expr)) * as.numeric(ncol(expr)))
      pcts[i] <- pct0
      message(pct0)
      rm(data)
      gc()
    }else{
      pcts[i] <- "2"
      rm(data)
      gc()
    }
  }
}
datasets$pct <- as.numeric(pcts)
datasets_sel <- datasets[!datasets$pct %in% c(0,2),]
datasets <- datasets[datasets$pct !=2,]



