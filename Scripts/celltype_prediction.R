## Utils##
source("utils.R")
##Packages##
library(magrittr)
library(scPred)
library(caret)
library(SingleR)
##Datasets##
datasets <- list()
##PBMC Dataset##
PBMC <- scPred::pbmc_1
counts <- PBMC@assays$RNA@counts
samplesheet <- data.frame("IDs" = rownames(PBMC@meta.data),"celltype" = make.names(PBMC$cell_type))
selectedCelltypes <- names(table(samplesheet$celltype)[table(samplesheet$celltype)>=100])
samplesheet <- samplesheet[samplesheet$celltype %in% selectedCelltypes,]
counts <- counts[,samplesheet$IDs]
counts <- counts[rowSums(counts>=1) >= 100,]
shuffled_counts <- apply(counts,2,function(cell){
  cell[cell>=1] <- sample(cell[cell>=1])
  return(cell)
})
datasets[["PBMC_shuffled"]] <- list(shuffled_counts,samplesheet)
datasets[["PBMC"]] <- list(counts,samplesheet)## M1 Dataset ##
M1 <- readRDS("C:/Users/gabouland/Documents/004 PhD/zeros_review/data/M1.rds")
counts <- M1[[1]]
samplesheet <- M1[[2]]
rm(M1)
samplesheet$celltype <- sapply(samplesheet$celltype,function(x)unlist(strsplit(x,split = " "))[1])
AstroIDs <- samplesheet[samplesheet$celltype == "Astro","IDs"]
MicroIDs <- samplesheet[samplesheet$celltype == "Micro","IDs"]
OPCIDs <- samplesheet[samplesheet$celltype == "OPC","IDs"]
toSamp <- (5000 - length(c(AstroIDs,MicroIDs,OPCIDs)))/3
OligoIDs <- samplesheet[samplesheet$celltype == "Oligo","IDs"] %>% sample(size = toSamp)
ExcIDs <- samplesheet[samplesheet$celltype == "Exc","IDs"]%>% sample(size = toSamp)
InhIDs <- samplesheet[samplesheet$celltype == "Inh","IDs"]%>% sample(size = toSamp)
selectedIDs <- c(AstroIDs,MicroIDs,OPCIDs,OligoIDs,ExcIDs,InhIDs)
samplesheet <- samplesheet[match(selectedIDs,samplesheet$IDs),]
counts <- counts[,samplesheet$IDs]
datasets[["M1"]] <- list(counts,samplesheet)

## Run to shuffle ##
shuffled_counts <- apply(counts,2,function(cell){
  cell[cell>=1] <- sample(cell[cell>=1])
  return(cell)
})
datasets[["M1_shuffled"]] <- list(shuffled_counts,samplesheet)
rm(list = ls()[ls()!="datasets"])
## AD Dataset ##
AD <- readRDS("C:/Users/gabouland/Documents/004 PhD/008scBinaryAnalysis/000 Preprocess/ADset.rds")
counts <- AD[[1]]
samplesheet <- AD[[2]]
samplesheet <- samplesheet[!samplesheet$celltype %in% c("unID","doublet","endo"),]
counts <- counts[,samplesheet$IDs]
datasets[["AD"]] <- list(counts,samplesheet)
## Run to shuffle ##
shuffled_counts <- apply(counts,2,function(cell){
  cell[cell>=1] <- sample(cell[cell>=1])
  return(cell)
})
datasets[["AD_shuffled"]] <- list(shuffled_counts,samplesheet)
rm(list = ls()[ls()!="datasets"])
## Run ##
## Binary vs counts run ##
methods <- c("counts_SingleR","binary_SingleR","counts_scPred","binary_scPred")
## Counts normal vs counts shuffle run ##
methods <- c("counts_SingleR","counts_scPred")
reps <- 1
perform <- lapply(1:reps,function(i){
  message(i)
  lapply(names(datasets),function(dataset_name){
    counts <- datasets[[dataset_name]][[1]]
    samplesheet <- datasets[[dataset_name]][[2]]
    trainingIDs <- sample(samplesheet$IDs,round(nrow(samplesheet) * 0.75))
    validationIDs <- samplesheet$IDs[!samplesheet$IDs %in% trainingIDs]
    trainingData <- counts[,trainingIDs]
    validationData <- counts[,validationIDs]
    trainingLabels <- samplesheet[match(trainingIDs,samplesheet$IDs),"celltype"]
    validationLabels <- samplesheet[match(validationIDs,samplesheet$IDs),"celltype"]
    lapply(methods,function(method){
      predictions <- runPredict(trainingData, validationData, trainingLabels,method)
      stats <- confusionMatrix(predictions,as.factor(validationLabels))
      accuracy <- stats$overall[1]
      medianF1 <- median(stats$byClass[,7],na.rm = T)
      out <- data.frame("dataset" = dataset_name,
                        "method" = method,
                        "global_accuracy" = unname(accuracy),
                        "medianF1" = unname(medianF1))
      print(out)
      return(out)
    }) %>% do.call(what = "rbind")
  }) %>% do.call(what = "rbind")
}) %>% do.call(what = "rbind")

