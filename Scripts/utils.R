aggregate_zeros <- function(sceObject){
  counts <- sceObject@assays@data$counts
  gi <- metadata(sceObject)$gene_info
  cellsheet <- data.frame("colnames" = colnames(counts),
                          "group" = sceObject@colData@listData$group_id,
                          "sampleID" = as.character(sceObject@colData@listData$sample_id))
  rownames(counts) <- paste0(rownames(counts),"-",gi$category)
  samplesIDs <- unique(cellsheet$sampleID)
  zero_aggr <- lapply(samplesIDs,function(ID){
    cells <- cellsheet[cellsheet$sampleID == ID,"colnames"]
    rowSums(counts[,cells] >=1) / length(cells)
  }) %>% do.call(what = cbind)
  colnames(zero_aggr) <- samplesIDs
  #factors <- median(colSums(zero_aggr)) / colSums(zero_aggr)
  #zero_aggr<- t(t(zero_aggr) * factors)
  #zero_aggr <- atanh(zero_aggr-0.0000001)+1
  return(zero_aggr)
}

aggregate_mean <- function(sceObject){
  res <- aggregateData(sceObject,by = "sample_id",fun = "mean")
  aggregated <- res@assays@data[[1]]
  rownames(aggregated) <- paste0(rownames(aggregated),"-",metadata(res)$gene_info$category)
  return(aggregated)
}

aggregate_sum <- function(sceObject){
  res <- aggregateData(sceObject,by = "sample_id",fun = "sum")
  aggregated <- res@assays@data[[1]]
  rownames(aggregated) <- paste0(rownames(aggregated),"-",metadata(res)$gene_info$category)
  return(aggregated)
}

aggregate_data <- function(sceObject, type){
  aggregated <- switch(type,
                       "zeros" = aggregate_zeros(sceObject),
                       "mean" = aggregate_mean(sceObject),
                       "sum" = aggregate_sum(sceObject))
  return(aggregated)
}


runTests <- function(Data, test){
  res <- switch(test,
                "limmavoom" = runLimma(Data),
                "limmatrend"= runLimmaTrend(Data),
                "edgeR" = runEdgeR(Data),
                "ttest" = runTtest(Data),
                "wilcox" = runWilcox(Data))
  
  
  res$true <- as.numeric(grepl("-de",rownames(res)))
  res$detected <- ifelse(res$fdr <= 0.05,1,0)
  return(res)
}

runLimmaTrend <- function(Data){
  require(limma)
  require(edgeR)
  groups <- strsplit(colnames(Data),split = "[.]") %>% sapply(FUN = function(x)x[2])
  dge <- DGEList(Data, group = groups)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~groups)
  logCPM <- cpm(dge, log=FALSE, prior.count=0.1)
  fit <- lmFit(logCPM, design = design)
  fit <- eBayes(fit,trend = TRUE)
  res <- data.frame("pvalue" = fit$p.value[,2])
  res$fdr <- p.adjust(res$pvalue,method = "BH")
  return(res)
}

runLimma <- function(Data){
  require(limma)
  require(edgeR)
  groups <- strsplit(colnames(Data),split = "[.]") %>% sapply(FUN = function(x)x[2])
  dge <- DGEList(Data, group = groups)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~groups)
  vm <- voom(dge, design = design)
  fit <- lmFit(vm, design = design)
  fit <- eBayes(fit)
  res <- data.frame("pvalue" = fit$p.value[,2])
  res$fdr <- p.adjust(res$pvalue,method = "BH")
  return(res)
}

runEdgeR <- function(Data){
  require(limma)
  require(edgeR)
  groups <- strsplit(colnames(Data),split = "[.]") %>% sapply(FUN = function(x)x[2])
  dge <- DGEList(Data, group = groups)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~groups)
  dge <- estimateDisp(dge, design = design)
  fit <- glmFit(dge, design = design)
  lrt <- glmLRT(fit)
  res <- topTags(lrt, n = nrow(Data))$table
  res$fdr <- p.adjust(res$PValue,method = "BH")
  return(res)
}


runTtest <- function(Data){
  res <- apply(Data,1,function(row){
    if(length(unique(row)) == 1){
      return(data.frame(t = NA, p = 1))
    }else{
      res_t <- t.test(row[grepl(".A",names(row))],row[grepl(".B",names(row))])
      data.frame(t = res_t$statistic, p = res_t$p.value)
    }
  }) %>% do.call(what = rbind)
  res$fdr <- p.adjust(res$p,method = "BH")
  return(res)
}


runWilcox <- function(Data){
  index <- grep(".A",colnames(Data))
  wilcox_P <- sapply(rownames(Data),function(gene){
    min(2 * min(limma::rankSumTestWithCorrelation(index = index, statistics = Data[gene,])), 1)
  })
  res <- data.frame("P" = wilcox_P, "fdr" = p.adjust(wilcox_P,method = "BH"))
  return(res)
  
}

getMetric <- function(y, x , metric){
  tp <- sum(y == x & y == 1)
  fp <- sum(x == 1 & y == 0)
  fn <- sum(x == 0 & y == 1)
  tn <- sum(x == y & y == 0)
  switch(metric,
         "TPR" = {
           tp / (tp +fn)
         },
         "FPR" = {
           fp / (fp + tn)
         },
         "PPV" = {
           tp / (tp +fp)
         },
         "F1" = {
           tp  / (tp + 0.5*(fp+fn))
         })
}


findFeatures <- function(counts, thr = 0.2){
  require(M3Drop)
  norm <- M3DropConvertData(counts, is.counts=TRUE)
  message("Finding variable genes")
  M3Drop_genes <- rownames(M3DropFeatureSelection(norm,
                                                  mt_method="fdr",
                                                  mt_threshold=thr,
                                                  suppress.plot = T))
  return(M3Drop_genes)
}

counts_SVM <- function(reference, target, labels){
  require(Seurat)
  require(magrittr)
  data <- cbind(reference,target)
  SeurObj <- CreateSeuratObject(counts = data)
  SeurObj <- SeurObj %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()
  trainData <- cbind(data.frame("celltype" = labels),
                     SeurObj@reductions$pca@cell.embeddings[colnames(reference),1:10])
  ctrl <- trainControl(method = "cv", savePred=T, classProb=T)
  mod <- train(celltype~., data=trainData, method = "svmLinear", trControl = ctrl)
  predicted <- predict(mod, SeurObj@reductions$pca@cell.embeddings[colnames(target),1:10])
  return(predicted)
}


counts_SingleR <- function(reference, target, labels){
  require(SingleR)
  data <- cbind(reference,target)
  Feat <- findFeatures(data)
  expr.norm <- t(t(data)/colSums(data))*10000
  normalized <- log(expr.norm + 1)
  reference <- normalized[Feat,colnames(reference)]
  target <- normalized[Feat,colnames(target)]
  predictions <- SingleR(test = target,ref = reference,labels = labels)
  return(as.factor(predictions$pruned.labels))
}

binary_SingleR <- function(reference, target, labels){
  require(SingleR)
  data <- cbind(reference,target)
  Feat <- findFeatures(data)
  binary <- (data>=1)*1
  reference <- binary[Feat,colnames(reference)]
  target <- binary[Feat,colnames(target)]
  predictions <- SingleR(test = target,ref = reference,labels = labels)
  return(as.factor(predictions$pruned.labels))
}

binary_KNN <- function(reference, target, labels){
  Feat <- findFeatures(cbind(reference,target))
  reference <- as.matrix((reference[Feat,] >=1)*1)
  target <- as.matrix((target[Feat,] >=1)*1)
  overlap <- crossprod(reference,target)
  cellIDs <- data.frame("IDs" = colnames(reference),"celltype" = labels)
  predictions <- apply(overlap,2,function(sum){
    NNs <- names(sort(sum,decreasing = T))[1:10]
    predictions <- cellIDs[match(NNs,cellIDs$IDs),"celltype"]
    return(names(sort(table(predictions),decreasing = T))[1])
  })
  return(as.factor(predictions))
}

binary_scPred <- function(reference, target, labels){
  require(scPred)
  require(Seurat)
  Feat <- findFeatures(cbind(reference,target))
  reference <- (reference[Feat,]>=1)*1
  target <- (target[Feat,]>=1)*1
  reference <- CreateSeuratObject(counts = reference)
  reference <- reference %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()
  
  reference$celltype <- labels
  reference <- getFeatureSpace(reference, "celltype")
  reference <- trainModel(reference)
  target <- CreateSeuratObject(counts = target)
  target <- NormalizeData(target)
  target <- scPredict(target, reference)
  return(as.factor(target$scpred_no_rejection))
}

counts_scPred <- function(reference, target, labels){
  require(scPred)
  require(Seurat)
  Feat <- findFeatures(cbind(reference,target))
  reference <- reference[Feat,]
  target <- target[Feat,]
  reference <- CreateSeuratObject(counts = reference)
  reference <- reference %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()
  
  reference$celltype <- labels
  reference <- getFeatureSpace(reference, "celltype")
  reference <- trainModel(reference)
  target <- CreateSeuratObject(counts = target)
  target <- NormalizeData(target)
  target <- scPredict(target, reference)
  return(as.factor(target$scpred_no_rejection))
}

runPredict <- function(reference, target, labels, method){
  predictions <- switch(method,
                        "counts_SVM" = counts_SVM(reference, target, labels),
                        "counts_SingleR"= counts_SingleR(reference, target, labels),
                        "binary_SingleR" = binary_SingleR(reference, target, labels),
                        "binary_KNN" = binary_KNN(reference, target, labels),
                        "counts_scPred" = counts_scPred(reference, target, labels),
                        "binary_scPred" = binary_scPred(reference, target, labels))
  return(predictions)
}




signed_louvain <- function(g, gamma = 1, weighted = T, comm = NULL, B = NULL, seed = NULL,
                           mod = c(NULL, 'modularity', 'potts', 'neg_sym', 'neg_asym'), class = TRUE) {
  
  # Fast and accurate multi-iterative implementation of theLouvain community detection algorithm. 
  # This implementation was adapted from the MatLab function community_louvain.m available on 
  # Brain Connectivity toolbox. https://sites.google.com/site/bctnet/
  # 
  # INPUT: 
  #          g,     a directed/undirected weighted/unweighted igraph object or matrix object
  #          gamma, (optional) resolution parameter. gamma > 1 small modules, 0 < gamma > 1, larger modules 
  #          comm,  (optional) initial community affiliation vector
  #          B,     (optional) custom objective matrix (must have the same dimensions as input graph)
  #          mod,   (optional) Objective function type. it can have the following inputs
  #                 'modularity'  Modularity (default)
  #                 'potts'       Potts-model Hamiltonian (for binary networks)
  #                 'neg_sym'     symmetric treatment of negative weights
  #                 'neg_asym',    asymmetric treatment of negative weights (see Rubinov and Sporns (2011))
  #                 custom objective-function matrix (under substitute(), bquote(), expression())
  #                 bquote(B <- (A - gamma*outer(rowSums(A),colSums(A),'*')/s)/s)
  #          class  Logical. TRUE (default) gives list list with class 'communities'. FALSE gives 
  #                 avector with the community assignments
  #
  # OUTPUT:   list of class 'communities' (if class == TRUE) or just the membership vector
  #           membership  optimal community structure 
  #           modularity  modularity value
  #           names       nodes' names
  #           algorithm   chr 'signed_louvain'
  
  stopifnot(is.igraph(g) | is.matrix(g));
  
  if(is.igraph(g)){                           # get adjacendy matrices
    if(isFALSE(weighted)){
      A <- as.matrix(as_adjacency_matrix(g));
    } else {
      A <- as.matrix(g[]);
    }
  } else {
    A <- g
  }
  
  if(!is.null(seed)){
    set.seed(seed);
  }
  
  n <- ncol(A) #get number of edges
  s <- sum(A)  #get sum of edges
  
  if(mod == 'neg_sym' | mod == 'neg_asym'){
    # getting graph with positive weights
    pos <- A;
    pos[pos < 0] <- 0;
    vpos <- sum(pos);
    Bpos <- pos - gamma*outer(rowSums(pos),colSums(pos), '*')/vpos; #positive modularity
    
    # getting graph with negative weights
    neg <- -A;
    neg[neg < 0] <- 0;
    vneg <- sum(neg);
    if(vneg != 0){
      Bneg <- neg - gamma*outer(rowSums(neg),colSums(neg), '*')/vneg; #negative modularity
    }
    else{
      Bneg <- 0;
    }
  }
  else{
    if(min(A) < -1e-10){
      stop('Input graph/matrix has negative weights.\n
           Specify "neg_sym" or "neg_asym" in mod argument')
    }
  }
  
  if(mod == 'potts' && any(a != as.logical(a))){
    stop('Potts-model Hamiltonian requires a binary network')
  }
  
  if(is.null(mod)){ 
    mod <- 'modularity' 
  }
  
  if(!is.null(B)){
    if(identical(dim(B), dim(A))){
      B <- as.matrix(B); 
    }
    else{
      stop('Graph arg and B arg must have the same dimensions')
    }
  }
  else{
    if(is.language(mod)){
      eval(mod) # Ex: substitute(B <- (A - gamma*outer(rowSums(A),colSums(A),'*')/s)/s), bquote(), expression()
    }
    else{
      if(mod == 'modularity'){
        B <- as.matrix((A - gamma*outer(rowSums(A),colSums(A),'*')/s)/s);
      }
      else{
        if(mod == tolower('potts')){
          B <- as.matrix(A - gamma*(!as.logical(A)));
        }
        else{
          if(mod == 'neg_sym'){
            B <- as.matrix(Bpos/(vpos + vneg) - Bneg/(vpos+vneg));
          }
          else{
            if(mod == 'neg_asym'){
              B <- as.matrix(Bpos/vpos - Bneg/(vpos+vneg));
            }
            else{
              stop('Must choose "modularity", "potts(binary only)", "neg_sym", "neg_asym"\n 
                   or input a valid objective function for modularity matrix B')
            }
          }
        }
      }
    }
  }
  
  if(is.null(comm)){
    M0 <- as.numeric(c(1:n));
  } 
  else{
    if(length(comm) != n){ 
      stop('Length of comm arg must be the same as number of nodes')
    }
    else{
      M0 <- comm
    }
  }
  
  Mb <- match(M0, unique(M0));
  M <- Mb
  
  B <- as.matrix((B + t(B))/2);                     # symmetrize modularity matrix
  Hnm <- matrix(0,n,n);                             # node-to-module degree
  for(i in 1:max(Mb)){                              # loop over modules
    Hnm[,i] <- rowSums(as.matrix(B[,Mb == i]))
  }
  
  Q0 <- -Inf;
  Q <- sum(B[outer(M0, M0, '==')]);                 # compute modularity
  first_iteration <- TRUE;
  
  while(Q - Q0 > 1e-10) {
    flag = TRUE;                                    # flag for within-hierarchy search
    
    while(isTRUE(flag)){
      flag = FALSE;
      
      for(i in sample(1:n, replace = F)){           #loop over all nodes randomically
        ma <- Mb[i];                                # current module i
        dQ <- Hnm[i,] - Hnm[i,ma] + B[i,i];
        dQ[ma] <- 0;                                # (line above) algorithm condition
        
        max_dQ <- max(dQ);                          # maximal increase in modularity
        mb <- which(dQ %in% max_dQ);                # and corresponding module
        if(max_dQ > 1e-10){
          flag <- TRUE;
          Mb[i] <- mb;                              # reassign module
          
          Hnm[,mb] <- Hnm[,mb] + B[, i];            # change node-to-module strengths
          Hnm[,ma] <- Hnm[,ma] - B[, i];
        }
      }
    }
    
    Mb <- match(Mb,unique(Mb));                     # new module assignments
    M0 <- M;
    if(isTRUE(first_iteration)){
      M <- Mb;
      first_iteration <- FALSE;
    }
    else{
      for(i in 1:n){                                # Loop through initial module assignments
        M[M0 == i] <- Mb[i];                        # assign new modules
      }
    }
    
    n <- max(Mb);                                   # new number of modules (weird), initially n was supposed to be nodes, not modules
    Bl <- matrix(0,n,n);                            # new weighted matrix
    for(i in 1:n){
      for(z in i:n){
        bm <- sum(B[Mb == i, Mb == z]);             # pull weights of nodes in same module
        Bl[i,z] <- bm;
        Bl[z,i] <- bm;
      }
    }
    B <- Bl;
    
    Mb <- seq_len(n);                               # initial module assignments
    Hnm <- B;                                       # node-to-module strength
    
    Q0 <- Q;
    Q <- psych::tr(B);                              # compute modularity
  }
  
  if(isTRUE(class)){
    M <- list(membership = M,
              modularity = Q,
              names = colnames(A),
              algorithm = 'signed_louvain')
    M <- structure(M, class = 'communities')  # To fit in other functions from igraph
  }
  return(M)
}
