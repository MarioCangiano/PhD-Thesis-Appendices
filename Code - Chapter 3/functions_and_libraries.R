
##########################################################################################

# source("D:/mcangiano/Glasgow/RNAseq_orthografts/variousScripts/RevealerFunctions.R")

########################################################################################## libraries

library(purrr)
library(survival)
library(survminer)
library("parmigene")
library("gProfileR")
library("clusterProfiler")
library("ggplot2")
library("Rcpp")
# library("minet")
library("igraph")
library("networktools")
library("kernlab")
library("mygene")
library("kernlab")
library('knitr')
library("snow")
# library("PerPAS")
library('limma')
library('reshape2')
library('RColorBrewer')
library('WGCNA')
# library("xlsx")
# library("Homo.sapiens")
library("gplots")
library("heatmap3")
library('igraph')
# library("ssmarina")
# library("viper")
# library("edgeR")
library("diggit")
#library("cogena")
library("DESeq2")
library("limma")
# library("NOISeq")
# library("RGBM")
# library("doMC")
# library("RTN")
# library("RTNduals")
library("parallel")
# library("VennDiagram")
library("eulerr")
library("EnrichmentBrowser")

########################################################################################## functions

findSignificantAssociations <- function(vect, meta, objective="PFI", objectiveTime="PFI.time", factor=T, cv=T, 
                                        covariates=c("Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer")){
  if(factor){
    vect <- as.factor(vect)
  }
  if(length(table(vect)) > 1){
    if(cv){
      a <- b <- c()
      for(j in 1:length(vect)){
        mySamples <- names(vect[-j])
        a2 <- summary(coxph(Surv(meta[mySamples, objectiveTime], meta[mySamples, objective]) ~ vect[-j] , data = meta[mySamples, ]))$coefficients[5]
        if(!is.null(covariates)){
          b2 <- summary(coxph(Surv(meta[mySamples, objectiveTime], meta[mySamples, objective]) ~ vect[-j] + meta[mySamples, covariates], data = meta[mySamples, ]))$coefficients[1, 5]
        }else{
          b2 <- NA 
        }
        a <- c(a, a2)
        b <- c(b, b2)
      }
      a <- mean(a)
      b <- mean(b)
    }else{
      a <- summary(coxph(Surv(meta[names(vect), objectiveTime], meta[names(vect), objective]) ~ vect , data = meta))$coefficients[5]
      if(!is.null(covariates)){
        b <- summary(coxph(Surv(meta[names(vect), objectiveTime], meta[names(vect), objective]) ~ vect + meta[names(vect), covariates], data = meta))$coefficients[1, 5]
      }else{
        b <- NA 
      }
    }
    res <- c(a, b)
  }else{
    res <- c(1,1)
  }
  if(cv){
    names(res) <- c("univariate_cv", "multivariate_cv")
  }else{
    names(res) <- c("univariate", "multivariate")
  }
  return(res)
}

prepareGRNforEB <- function(x){
  df <- c()
  for(i in 1:length(x)){
    a <- cbind(rep(names(x)[i], length(x[[i]]$pos)), x[[i]]$pos, rep("+", length(x[[i]]$pos)))
    b <- cbind(rep(names(x)[i], length(x[[i]]$neg)), x[[i]]$neg, rep("-", length(x[[i]]$neg)))
    df <- rbind(df, a, b)
  }
  colnames(df) <- c("FROM", "TO", "TYPE")
  return(df)
}

makeJMJDGraph <- function(degs){
  
  load("/mnt/data/mcangiano/Transpot/7060_Glasgow/RNAseq/Data_analyis/networkAnalysis/RTN/standardisedCounts_pearson/RTNsingleRegAracne.shadowFiltered.rda")
  myGenes <- c("ENSG00000070495_JMJD6_17:74708918-74722866", "ENSG00000081692_JMJD4_1:227918125-227923112")
  netRegAttempt <- netRegShadowFiltered[myGenes] # NULL results as expected
  netRegAttempt2 <- c(netRegAttempt, netRegShadowFiltered[unlist(netRegAttempt)])
  library(igraph)
  myNodes <- unique(unname(unlist(netRegAttempt2)))
  myDegs <- degs[sapply(myNodes, function(x){unlist(strsplit(x, "_"))[1]}), ]
  rownames(myDegs) <- sapply(myNodes, function(x){unlist(strsplit(x, "_"))[2]})
  upGenes <- rownames(myDegs)[myDegs$pvalue < 0.05 & myDegs$log2FoldChange > 0 & !is.na(myDegs$pvalue)]
  downGenes <- rownames(myDegs)[myDegs$pvalue < 0.05 & myDegs$log2FoldChange < 0 & !is.na(myDegs$pvalue)]
  myGrn <- prepareGRNforEB(netRegAttempt2)
  myGraph <- graph_from_data_frame(d = myGrn, directed = F)
  V(myGraph)$degree <- round(log2(degree(myGraph)), digits = 0) + 1
  load("/mnt/data/mcangiano/Transpot/7060_Glasgow/RNAseq/Data_analyis/networkAnalysis/startData/tfListFrancesca.rda")
  newTFs <- newTFs[!is.na(newTFs)] 
  newTFs <- newTFs[!duplicated(newTFs)]
  V(myGraph)$name <- sapply(V(myGraph)$name, function(x){unlist(strsplit(x, "_"))[2]})  
  V(myGraph)$nodeColor[V(myGraph)$name %in% names(newTFs)] <- "yellow"
  V(myGraph)$nodeColor[V(myGraph)$name %in% upGenes] <- "red"
  V(myGraph)$nodeColor[V(myGraph)$name %in% downGenes] <- "blue"
  plot(myGraph, vertex.label=ifelse(V(myGraph)$name %in% names(newTFs), V(myGraph)$name, NA), 
       vertex.size=V(myGraph)$degree, vertex.color=V(myGraph)$nodeColor,layout=layout_with_kk(myGraph), 
       edge.width=2, vertex.label.cex = 1, vertex.label.font=2, vertex.label.color="black",
       edge.color=ifelse(E(myGraph)$TYPE == "+", "firebrick", "royalblue4"))
  return(myGraph)
}

findBestGroup <- function(mat, x, meta, objective="DFI", objectiveTime="DFI.time", factor=T, m=3, 
                          covariates=c("Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer")){
  allCombs <- t(combn(x = x, m = m))
  # pos1 <- which(allCombs == "ENSG00000070495_JMJD6_17:74708918-74722866", arr.ind = T)[, 1]
  # pos2 <- which(allCombs == "ENSG00000081692_JMJD4_1:227918125-227923112", arr.ind = T)[, 1]
  # pos <- intersect(pos1, pos2)
  pos <- which(allCombs == "ENSG00000070495_JMJD6_17:74708918-74722866", arr.ind = T)[, 1]
  allCombs <- allCombs[pos, ]
  
  allP <- namesAll <- c()
  for(i in 1:nrow(allCombs)){
    columns <- allCombs[i, ]
    namesAll <- c(namesAll, paste(sort(columns), collapse = "and"))
    myReg2 <- rowSums(mat[, columns])
    myReg2[myReg2 > 1] <- 1
    myMeta0 <- cbind(meta[names(myReg2), ], as.factor(myReg2))
    a <- summary(coxph(Surv(meta[names(myReg2), objectiveTime], meta[names(myReg2), objective]) ~ myReg2 , data = meta))$coefficients[5]
    b <- summary(coxph(Surv(meta[names(myReg2), objectiveTime], meta[names(myReg2), objective]) ~ myReg2 + meta[names(myReg2), covariates], data = meta))$coefficients[1, 5]
    res <- c(a, b)
    allP <- cbind(allP, res)
  }
  colnames(allP) <- namesAll
  rownames(allP) <- c("univariate", "multivariate")
  return(allP)
}

makeBromoScore <- function(fpkms, sig="F:/mcangiano/ICGC/startData/signatures_symbols.gmt"){
  library(GSVA)
  library(GSEABase)
  library(gage)  # this package for "readList" function for reading in gmt file
  library(impute)  # this package required if data matrix contains missing data to impute
  gs = readList(sig)
  data <- as.matrix(fpkms, dimnames=row.names(fpkms))
  gsva.es <- gsva(data, gs, min.sz = 1, max.sz = 500, verbose = TRUE)
  return(gsva.es)
}


findPathSignature <- function(mat, sample= NULL, cpu=8, myCut= 0.1,
                              pd= myPd, fdr= FALSE, myAlpha= 1, onlyBimodal= FALSE){
  
  library(geneSignatureFinder)
  library(survival)
  
  message("step 1: finding seeds")
  aNCPUS <- NCPUS(nchips= cpu)
  if(length(sample) == 0){
    geData <- mat
    save(geData, file="geData.rda")
    load(file="geData.rda", envir= .GlobalEnv)
  }else{
    message("selecting samples...")
    geData <- mat[sample, ]
    save(geData, file="geData.rda")
    load(file="geData.rda", envir= .GlobalEnv)
    message("done")
  }
  months <- pd$Time
  status <- pd$Endpoint
  names(months) <- names(status) <- rownames(pd)
  stData <- Surv(months[substr(x=rownames(geData), start=1, stop=12)],
                 status[substr(x=rownames(geData), start=1, stop=12)])
  save(stData, file="stData.rda")
  load(file="stData.rda", envir= .GlobalEnv)
  rownames(stData) <- rownames(geData)
  
  seeds <- seedsFinder(cutoff= 1.95, cpuCluster= aNCPUS)
  if(onlyBimodal == TRUE){
    bimodalCells <- seeds[, "bic1"] > seeds[, "bic2"]
  }else{
    bimodalCells <- seeds[, "pValue"] >= 0 
  }
  if(fdr == TRUE){
    significantCells <- BHcorrection(seeds[, "pValue"], alpha= myAlpha) < myAlpha
  }else{
    significantCells <- seeds[, "pValue"] < myAlpha
  }
  
  seedCells <- which(bimodalCells * significantCells == 1)
  
  if(length(seedCells) > 0 ){
    message("step 2: finding signature")
    K <- length(seedCells)
    signature <- vector("list", K)
    names(signature) <- names(seedCells)
    aNCPUS <- NCPUS(nchips= cpu)
    for(k in 1:K){
      sig <- signatureFinder(seedCells[k], cpuCluster=aNCPUS, stopCpuCluster=FALSE)
      imp <- importance(aSignatureFinder=sig, cpuCluster=aNCPUS, stopCpuCluster=FALSE)
      signature[[k]] <- imp 
    }
    stopCluster(cl=aNCPUS)
    attr(signature, "Creation date") <- date()
    
    message("step 3: pruning signatures")
    for(k in 1:K){
      if(length(signature[[k]]$signature) > 2){
        repeat{
          tmp <- removeGeneFrom(aSignatureFinder=signature[[k]], cutoff=myCut)
          if(is.null(tmp)){
            break
          }else{
            signature[[k]] <- tmp
          }
        }
      }
    }
    
    message("step 4: results preparation")
    if(length(signature) > 1){
      result <- searchResultsSummaryTable(aSearchResults=signature)
      ensamble <- ensembleTable(aSearchResults=signature)
    }else{
      result <- signature
      ensamble <- signature[[1]]$importance
    }
    
    samples <- rownames(geData)
    parameters <- c(fdr, onlyBimodal, myAlpha, myCut)
    names(parameters) <- c("fdr_correction", "onlyBimodal", "alpha", "importance_cut")
    allResult <- list(signature, samples, stData, result, ensamble, parameters)
    names(allResult) <- c("signature","samples", "stData", "summaryTable", "ensambleTable", "parameters_adopted")
  }else{
    allResult <- "empty"
    message("no significant and bimodal elements!")
  }
  
  file.remove("geData.rda", "stData.rda")
  rm(geData, stData)
  return(allResult)
}

signatureClustering <- function(searchResults, upperTrim = 0.8, lowerTrim = 0.2, stData){
  
  K <- length(searchResults$signature)
  n <- length(searchResults$sample)
  searchResults <- searchResults$signature
  simpleMatch <- function(x1, x2) return(1 - sum((x1 * x2) + abs(x1-1)*abs(x2-1))/length(x1))
  clusters <- matrix(0, ncol = n, nrow = K)
  rownames(clusters) <- paste(1:K)
  colnames(clusters) <- searchResults$signature
  for(k in 1:K) {
    rownames(clusters)[k] <- searchResults[[k]]$signatureName
    clusters[k, searchResults[[k]]$classification == "poor"] <- 1
  }
  dissimilarity <- matrix(0, K, K)
  rownames(dissimilarity) <- rownames(clusters)
  colnames(dissimilarity) <- rownames(clusters)
  for(k in 1:(K-1))
    for(kk in (k+1):K)
      dissimilarity[k, kk] <- simpleMatch(clusters[k,], clusters[kk,])
  dissimilarity <- dissimilarity + t(dissimilarity)
  dissimilarity <- as.dist(dissimilarity)
  
  plot(signatureHC <- hclust(dissimilarity),
       main ="Signatures Dendrogram", sub = "")
  tmp <- sort(bars <- apply(clusters, 2, mean))
  ccol <- colorRampPalette(c("red", "yellow"))(101)[round(tmp*100)+1]
  barplot(tmp, xlab = "good <------------------> poor", main = "samples' classification",
          space = 0, col = ccol, border = NA, names.arg = "")
  abline(h = lowerTrim, col = "grey50", lty = "dashed")
  abline(h = upperTrim, col = "grey50", lty = "dashed")
  lines(c(0, n), c(0, 1), col = "grey50", lty = "dashed")
  
  samplesClassification <- rep("uncertain", n)
  samplesClassification[bars > upperTrim] <- "poor"
  samplesClassification[bars < lowerTrim] <- "good"
  sf <- survfit(stData ~ samplesClassification)
  # tValue <- survdiff(stData[samplesClassification != "uncertain"] ~ samplesClassification[samplesClassification != "uncertain"])$chisq
  # pValue <- 1 - pchisq(tValue, 1)
  plot(sf, col = c("deepskyblue4", "darkmagenta", "grey50"), xlab = "months", ylab = "survival curves", 
       main = "samples classification")
  legend("bottomright", legend = c(paste("good =", table(samplesClassification)["good"]),
                                   paste("poor = ", table(samplesClassification)["poor"]), 
                                   paste("uncertain = ", table(samplesClassification)["uncertain"])), 
         fill = c("deepskyblue4", "darkmagenta", "grey50"))
  
}

extractSignatures <- function(x, stData){
  
  library(survival)
  pos <- x$ensambleTable$wImportance != 0 & x$ensambleTable$wImportance != "NaN"
  weightedSig <- as.numeric(as.character(x$ensambleTable$wImportance[pos]))
  names(weightedSig) <- rownames(x$ensambleTable)[pos]
  bestSig <- rownames(x$summaryTable)[order(x$summaryTable$`log(pValue)`, x$summaryTable$`tValue improvement`, 
                                            x$summaryTable$length, decreasing = T)][1]
  # bestSig <- rownames(x$summaryTable)[which.min(as.numeric(as.character(x$summaryTable$`log(pValue)`)))]
  classes <- as.character(x$signature[[bestSig]]$classification)
  names(classes) <- x$samples
  fullSig <- x$signature[[bestSig]]$signature
  #stData <- x$stData
  aSurv <- survfit(stData[names(classes),] ~ classes)
  library(survminer)
  ggsurvplot(cumevents = T, cumcensor = T, aSurv, data = stData, conf.int = F, pval=T, title=paste(bestSig, 
                                                                                                   " - len: ", length(fullSig), " - good: ", table(classes)[1], " poor: ", table(classes)[2], sep=""))
  res <- list(weightedSig=weightedSig, bestSig=bestSig, fullSig=fullSig, aSurv=aSurv)
  return(res)
}

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

correctForAliases <- function(x){
  require(limma)
  aliases <- alias2SymbolTable(alias = x, species = "Hs")
  aliases[is.na(aliases)] <- x[is.na(aliases)]
  message(sum(!x %in% aliases), " aliases corrected!")
  return(aliases)
}

prepareEnrichments <- function(x, threshold = 0.05, foldChange=1.5, type="rankedList", collectionId="symbol", regulons=T, 
                               ncores=8, eset=NULL, convert=F, removePoint=T, fullList=F){
  
  if(removePoint){
    message("removing points")
    if(length(dim(x)) == 2){
    a <- sapply(rownames(x), function(x){unlist(strsplit(x, "\\."))[1]})
    pos <- which(!duplicated(a))
    rownames(x)[pos] <- a[pos]
    x <- x[pos, ]
    }else{
      a <- sapply(x, function(x){unlist(strsplit(x, "\\."))[1]})
      x <- a[!duplicated(a)]
    }
    }
  
  if(convert){
    message("converting names")
    load("/mnt/data/mcangiano/7060_Glasgow/RNAseq/Data_analyis/networkAnalysis/startData/convertedGenesFromRawCountsFromMygene.rda")
    conv <- ensToComplete[rownames(x)]
    if(length(grep("more", conv)) > 0){
      x <- x[-1* grep("more", conv), ]
      conv <- conv[-1* grep("more", conv)]
    }
    pos <- which(!is.na(conv))
    rownames(x)[pos] <- conv[pos]
    x <- x[pos, ]  
  }
  
  cleanGmt <- function(x){
    x$ont <- as.factor(sapply(as.character(x$ont), function(y){unlist(strsplit(y, split = "_"))[1]}))
    x$gene <- sapply(as.character(x$gene), function(y){unlist(strsplit(y, split = "_"))[1]})
    return(x)
  }
  
  library(clusterProfiler)
  if(collectionId == "entrez"){
    h <- read.gmt(gmtfile = "/mnt/data/mcangiano/7060_Glasgow/RNAseq/Data_analyis/networkAnalysis/startData/h.all.v6.1.entrez.gmt")
    c3 <- read.gmt(gmtfile = "/mnt/data/mcangiano/7060_Glasgow/RNAseq/Data_analyis/networkAnalysis/startData/c3.tft.v6.1.entrez.gmt")
  }else{
    h <- read.gmt(gmtfile = "D:/mcangiano/Glasgow/RNAseq_orthografts/startData/h.all.v6.2.symbols.gmt")
    c3 <- read.gmt(gmtfile = "D:/mcangiano/Glasgow/RNAseq_orthografts/startData/c3.all.v6.2.symbols.gmt") 
  }
  if(type == "rankedList"){
    pos <- which(x$padj < threshold & abs(x$log2FoldChange) > foldChange)
    if(fullList){
    rankedList <- -log10(x[, "padj"])
    sign <- ifelse(x[, "log2FoldChange"] > 0, "+", "-")
    names(rankedList) <- rownames(x)
    }else{
    rankedList <- -log10(x[pos, "padj"])
    sign <- ifelse(x[pos, "log2FoldChange"] > 0, "+", "-")
    names(rankedList) <- rownames(x)[pos]
    }
    rankedPos <- rankedList[sign=="+"]
    rankedNeg <- rankedList[sign=="-"]
    rankedPos <- rankedPos[!is.na(rankedPos)]
    rankedNeg <- rankedNeg[!is.na(rankedNeg)]
    if(regulons){
      
      # rgbmPos <- cleanGmt(read.gmt(gmtfile = "RGBM/normCounts/netRegPos.gmt"))
      # rgbmNeg <- cleanGmt(read.gmt(gmtfile = "RGBM/normCounts/netRegNeg.gmt"))
      allRgbm <- cleanGmt(read.gmt(gmtfile = "RGBM/normCounts/netRegGeneral.gmt"))
      allAracne <- cleanGmt(read.gmt(gmtfile = "aracne/filteredNormCountsMutualInformation/additive.generalAracne.gmt"))
      # allAracne2 <- cleanGmt(read.gmt(gmtfile = "aracne/filteredNormCountsMutualInformation/additiveGeneral.gmt"))
      rtnAll <- cleanGmt(read.gmt(gmtfile = "RTN/standardisedCounts_pearson/RTNsingleReg.generalAracne.gmt"))
      rtnAll2 <- cleanGmt(read.gmt(gmtfile = "RTN/standardisedCounts_pearson/RTNsingleGeneral.gmt"))
      
      overRepAnalysis <- enricher(gene = names(rankedList)[pos], TERM2GENE = allRgbm, pvalueCutoff = 1)
      overRepAnalysis2 <- enricher(gene = names(rankedList)[pos], TERM2GENE = allAracne, pvalueCutoff = 1)
      overRepAnalysis3 <- enricher(gene = names(rankedList)[pos], TERM2GENE = rtnAll, pvalueCutoff = 1)
      overRepAnalysis4 <- enricher(gene = names(rankedList)[pos], TERM2GENE = rtnAll2, pvalueCutoff = 1)
      # gseaAnalysis3 <- GSEA(geneList = sort(rankedNeg, decreasing = T), TERM2GENE = rgbmNeg, verbose = T, pvalueCutoff = 1)
      
      results <- list(enirchRtn= overRepAnalysis3, enirchRtn2= overRepAnalysis4, enrichAracne= overRepAnalysis2,
                      enrichRgbm= overRepAnalysis, thresholdEnrichRankedList=threshold, thresholdFoldChange=foldChange)
    }else{
      if(collectionId == "symbol"){
        names(rankedPos) <- sapply(names(rankedPos), function(x){unlist(strsplit(x, "_"))[2]})
        names(rankedNeg) <- sapply(names(rankedNeg), function(x){unlist(strsplit(x, "_"))[2]})
        # rankedPos <- rankedPos[!is.na(names(rankedPos))]
        # rankedNeg <- rankedNeg[!is.na(names(rankedNeg))]
        gseaAnalysis <- GSEA(geneList = sort(rankedPos, decreasing = T), TERM2GENE = c3, verbose = T, pvalueCutoff = 1)
        gseaAnalysis2 <- GSEA(geneList = sort(rankedNeg, decreasing = T), TERM2GENE = c3, verbose = T, pvalueCutoff = 1)
        gseaAnalysis3 <- GSEA(geneList = sort(rankedPos, decreasing = T), TERM2GENE = h, verbose = T, pvalueCutoff = 1)
        gseaAnalysis4 <- GSEA(geneList = sort(rankedNeg, decreasing = T), TERM2GENE = h, verbose = T, pvalueCutoff = 1)
        overRepAnalysis <- enricher(gene = sapply(names(rankedList)[pos], function(x){unlist(strsplit(x, "_"))[2]}), TERM2GENE = h, pvalueCutoff = 1)
        results <- list(c3Up=gseaAnalysis, c3Down=gseaAnalysis2, Hup=gseaAnalysis3, Hdown=gseaAnalysis4, enrichH=overRepAnalysis, threshold=threshold)
      }else{
        translated <- bitr(geneID = names(rankedList), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
        myDup <- duplicated(translated$ENSEMBL)
        message("duplicated: ", sum(myDup))
        translated <- translated[!myDup, ]
        rownames(translated) <- translated$ENSEMBL
        rankedList <- rankedList[names(rankedList) %in% translated$ENSEMBL]
        names(rankedList) <- translated[names(rankedList), "ENTREZID"]
        logNegRanked <- -log(rankedList)
        gseaAnalysis <- GSEA(geneList = sort(logNegRanked, decreasing = T), TERM2GENE = c3, verbose = T, pvalueCutoff = 1)
        overRepAnalysis <- enricher(gene = names(logNegRanked), TERM2GENE = h, pvalueCutoff = 1)
        overRepAnalysis2 <- enricher(gene = names(logNegRanked), TERM2GENE = c3, pvalueCutoff = 1)
        NES <- gseaAnalysis@result$NES
        names(NES) <- rownames(gseaAnalysis@result)
        results <- list(gsea=gseaAnalysis, enrichH=overRepAnalysis, enrichC3=overRepAnalysis2, threshold=threshold, negLogRanked=logNegRanked, NES=NES)
      }
    }
  }
  if(type == "enrichment"){
    threshold <- logNegRanked <- NES <- NA
    gseaAnalysis <- NA
    overRepAnalysis <- enricher(gene = x, TERM2GENE = h, pvalueCutoff = 1)
    overRepAnalysis2 <- enricher(gene = x, TERM2GENE = c3, pvalueCutoff = 1)
    results <- list(gsea=gseaAnalysis, enrichH=overRepAnalysis, enrichC3=overRepAnalysis2, threshold=threshold, negLogRanked=logNegRanked, NES=NES)
  }
  if(type == "ssgsea"){
    if(regulons){
      rgbmPos <- gmt2list(annofile = "RGBM/normCounts/netRegPos.gmt")
      rgbmNeg <- gmt2list(annofile = "RGBM/normCounts/netRegNeg.gmt")
      rtnAll <- gmt2list(annofile = "RTN/standardisedCounts_pearson/RTNsingleRegGeneral.gmt")
      aracnePos <- gmt2list(annofile = "aracne/filteredNormCountsMutualInformation/additivePos.gmt")
      aracnePos <- aracnePos[!duplicated(names(aracnePos))]
      names(aracnePos) <- paste(sapply(names(aracnePos), function(x){paste(unlist(strsplit(x, "_"))[2:3], collapse="_")}), "_pos", sep="")
      aracneNeg <- gmt2list(annofile = "aracne/filteredNormCountsMutualInformation/additiveNeg.gmt")
      aracneNeg <- aracneNeg[!duplicated(names(aracneNeg))]
      aracneAll <- c(aracnePos, aracneNeg)
      names(aracneNeg) <- paste(sapply(names(aracneNeg), function(x){paste(unlist(strsplit(x, "_"))[2:3], collapse="_")}), "_neg", sep="")
      names(rgbmPos) <- paste(names(rgbmPos), "_pos", sep="")
      names(rgbmNeg) <- paste(names(rgbmNeg), "_neg", sep="")
      rgbmAll <- c(rgbmPos, rgbmNeg)
      eset2 <- eset
      rownames(eset2) <- sapply(rownames(eset2), FUN = function(x){unlist(strsplit(x, "_"))[1]})
      h <- GSVA::gsva(expr = eset2, gset.idx.list = rgbmAll, method="ssgsea", parallel.sz = ncores)
      h2 <- GSVA::gsva(expr = eset, gset.idx.list = rtnAll, method="ssgsea", parallel.sz = ncores)
      h3 <- GSVA::gsva(expr = eset, gset.idx.list = aracneAll, method="ssgsea", parallel.sz = ncores)
      results <- list(rgbmSsgsea= h, rtnSsgsea= h2, aracneSsgsea= h3)
      }
  }
  return(results)
}

prepareDataForVenn <- function(x, myCon="prova.txt", names=F){
  myNames <- c()
  for(i in 1:length(x)){
    if(names){
    myNames <- c(myNames, paste(names(x)[i], ":", paste(names(x[[i]]), collapse = ","), ";", sep=""))
    }else{
    myNames <- c(myNames, paste(names(x)[i], ":", gsub(":", "-", paste(x[[i]], collapse = ",")), ";", sep=""))
    }
  }
  writeLines(text = myNames, con = myCon)
}

prepareDataForVenn2 <- function(x){
  
  # x[[6]] <- NA
  # fit <- euler(c("A" = length(x[[1]]), "B" = length(x[[2]]), "C" = length(x[[3]]),
  #                "D" = length(x[[4]]), "E" = length(x[[5]]), "F" = length(x[[6]]),
  #                "A&B" = length(intersect(x[[1]], x[[2]])),
  #                "A&C" = length(intersect(x[[1]], x[[3]])),
  #                "A&D" = length(intersect(x[[1]], x[[4]])),
  #                "A&E" = length(intersect(x[[1]], x[[5]])),
  #                "A&F" = length(intersect(x[[1]], x[[6]])),
  #                "B&C" = length(intersect(x[[2]], x[[3]])),
  #                "B&D" = length(intersect(x[[2]], x[[4]])),
  #                "B&E" = length(intersect(x[[2]], x[[5]])),
  #                "B&F" = length(intersect(x[[2]], x[[6]])),
  #                "C&D" = length(intersect(x[[3]], x[[4]])),
  #                "C&E" = length(intersect(x[[3]], x[[5]])),
  #                "C&F" = length(intersect(x[[3]], x[[6]])),
  #                "D&E" = length(intersect(x[[4]], x[[5]])),
  #                "D&F" = length(intersect(x[[4]], x[[6]])),
  #                "F&E" = length(intersect(x[[5]], x[[6]]))),
  #              shape = "ellipse", input = "union")
  # plot(fit)
} # to be completed

unregister <- function(){
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

findFullNames <- function(x, myNames, ncores=16){
  
  fullNames <- c()
  pb <- txtProgressBar(min = 0, max = length(x), style = 3)
  
  if(ncores > 1){
    library(foreach)
    library(doSNOW)
    C <- ncores
    cl <- makeCluster(C)
    registerDoSNOW(cl)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    fullNames <- foreach(i= 1:length(x), .combine=c, .options.snow = opts)%dopar%{
      d <- grep(paste("^", x[i], "_", sep = ""), myNames)
      if(length(d) == 1){
        r <- myNames[d]
      }else{
        if(length(d) > 1){
          r <- "more_than_1_occurrence"
        }else{
          r <- NA
        }
      }
      return(r)
    }
    close(pb)
    stopCluster(cl)
  }else{
  for(i in 1:length(x)){
    d <- grep(paste("^", x[i], "_", sep = ""), myNames)
    if(length(d) == 1){
      fullNames <- c(fullNames, myNames[d])
    }else{
      if(length(d) > 1){
        fullNames <- c(fullNames, "more_than_1_occurrence")
      }else{
        fullNames <- c(fullNames, NA)
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  }
  
  names(fullNames) <- x
  return(fullNames)
}

plotCountsRange <- function(x, myRange, ylim=c(0,1.5), logAlready=F){
  if(!logAlready){
    x <- log10(x)
  }
  plot(density(x[x[, 1] > myRange[1] & x[, 1] < myRange[2], 1]), 
       main= paste("log non-zero x - ", myRange[1], "-", myRange[2], sep = ""),
       ylim=ylim)
  for(i in 2:ncol(x)){
    lines(density(x[x[, i] > myRange[1] & x[, i] < myRange[2], i]), col=i)
  }
}

plotCountsStatistics <- function(x, myRange, myGroups=NA, coloured=F, myTitle="counts"){
  allMeans <- apply(x, 1, mean)
  allSds <- apply(x, 1, sd)
  pos <- which(allMeans > myRange[1] & allMeans < myRange[2])
  if(coloured){
    load("startData/tfListFrancesca.rda")
    load("/opt/shared/hercules/genomes/GRCh37.75/Human_GRCh37.75.biotype.rda")
    newTFs <- newTFs[!is.na(newTFs)]
    tfEnsembl <- sapply(newTFs, FUN = function(x){unlist(strsplit(x, "_"))[1]})
    tfEnsembl <- tfEnsembl[tfEnsembl!="more"]
    biotype[tfEnsembl] <- "tf"
    biotypeColours <- rep("black", length(table(biotype)))
    names(biotypeColours) <- names(table(biotype))
    biotypeColours["tf"] <- "blue"
    biotypeColours["antisense"] <- "red"
    biotypeColours["lincRNA"] <- "yellow"
    biotypeColours["miRNA"] <- "green"
    biotypeColours["rRNA"] <- "purple"
    biotypeColours["processed_transcript"] <- "orange"
    myCols <- biotypeColours[biotype[names(pos)]]
  }else{
    myCols <- "black"
  }
  plot(allMeans[names(pos)], allSds[names(pos)], main = paste("means vs sds - ", myTitle," - ", myRange[1], "-", myRange[2], sep = ""),
       xlab="means", ylab="standard deviation", sub=paste("number of genes - ", length(pos)), lwd=3, col=myCols)
  
  if(!is.na(myGroups)){
    colnames(x) <- myGroups[colnames(x)]
    a <- aggregate(t(x), by = list(myGroups), FUN=mean) # find a way to speed it up
    b <- aggregate(t(x), by = list(myGroups), FUN=sd)
    
    # classes <- unique(names(myGroups))
    # for(i in 1:length(classes)){
    # allMeans <- apply(normCounts[, myGroups[names(myGroups) %in% classes[i]]], 1, mean)
    # allSds <- apply(normCounts[, myGroups[names(myGroups) %in% classes[i]]], 1, sd)
    # points(allMeans[names(pos)], allSds[names(pos)], col=i+1, pch=16)
    # }
    #   myAnnot <- cbind(cellLine=as.factor(names(myGroups)))
    #   rownames(myAnnot) <- colnames(x)
    #   heatmap3(cor(x[pos, ]), ColSideColors=myAnnot, trace='none', 
    #            main=paste('Sample correlations - normCounts', myRange[1], "-", myRange[2], sep = ""), scale = "none", showRowDendro = F, labCol = NA)  
    #   legend("bottomleft", legend = unique(names(myGroups)), fill = unique(as.numeric(as.factor(names(myGroups)))))
    # 
  }
} # for groups is better to show LOWESS trends

constructRegulon2 <- function(net, corNet, minReg = 20){
  
  net <- net[rowSums(net != 0) > minReg, ]
  if(nrow(net) > 0){
    netReg2 <- vector("list", nrow(net))
    names(netReg2) <- rownames(net)
    toKeep <- NULL
    for(k in 1:nrow(net)){
      tf <- rownames(net)[k]
      reg <- colnames(net[tf, net[tf, ] != 0, drop = F])
      
      tmp <- corNet[tf, reg, drop = F]
      
      pos <- colnames(tmp[1, tmp[1, ] > 0, drop = F])
      neg <- colnames(tmp[1, tmp[1, ] < 0, drop = F])
      
      if(length(pos) < minReg & length(neg) < minReg) next
      if(length(pos) <= 1) pos <- c("dummy","dummy")
      if(length(neg) <= 1) neg <- c("dummy","dummy")
      
      toKeep <- c(toKeep, k)
      netReg2[[k]] <- vector("list", 2)
      names(netReg2[[k]]) <- c("pos", "neg")
      netReg2[[k]]$pos <- pos
      netReg2[[k]]$neg <- neg
    }
    netReg2 <- netReg2[toKeep]
    print(length(netReg2))
  }else{
    metReg2 <- NA
    print(NA)
  }
  return(netReg2)
}

prepareFilteredRegulons <- function(out, z, tf_genes, minReg=20){
  
  out <- out[rownames(out) %in% tf_genes, ] 
  i <- 0
  while((sum(out>0)/length(tf_genes))>100 | i < 10)
  {
    i <- i+1
    message(i)
    median_weight <- median(out[out>0]);
    out[out<median_weight] <- 0;
  }
  corNet <- cor(t(z[rownames(out), ]), t(z), method = "spearman")
  netReg <- constructRegulon2(net = out, corNet, minReg = minReg)
  return(netReg)
}

extractDegs <- function(x, threshold=0.05){
  out <- x$log2FoldChange[x$padj < threshold]
  names(out) <- rownames(x)[x$padj < threshold]
  out <- out[!is.na(out)]
  return(out)
}

calculateIntersections <- function(x){
  M <- M2 <- D <- D2 <- matrix(0, nrow=length(x), ncol=length(x))
  for(i in 1:length(x)){
    for(j in 1:length(x)){
      common <- intersect(names(x[[i]]), names(x[[j]]))
      M[i, j] <- length(common)/length(x[[i]])
      M2[i, j] <- length(common)
      D[i, j] <- dist(x = rbind(x[[i]][common], x[[j]][common]))
      D2[i, j] <- sum(abs(x[[i]][common] - x[[j]][common]))
    }
  }
  rownames(M) <- colnames(M) <- rownames(M2) <- colnames(M2) <- rownames(D) <- colnames(D) <- rownames(D2) <- colnames(D2) <- names(x) 
  return(list(relativeNumbers=M, absoluteNumbers=M2, distance=D, difference=D2))
}

prepareIntersection <- function(x, y, i){
  # library(Hmisc)
  # a <- list()
  a$pos <- intersect(x[[i]]$pos, y[[i]]$pos)/sum(c(length(x[[i]]$pos), lenght(y[[i]]$pos)))
  a$neg <- intersect(x[[i]]$neg, y[[i]]$neg)/sum(c(length(x[[i]]$neg), lenght(y[[i]]$neg)))
  # label(a) <- paste(names(x)[i], names(y)[i], sep="_")
  return(a)
  }

makeGmt <- function(x, outFile, type="RGBM"){
  
  library(cogena)
  if(type == "RTN"){
    x <- x[lapply(x, length)>0]
    gmtlist2file(gmtlist = x,filename = paste(outFile, ".RTN.gmt", sep=""))
  }else{
    if(type == "all"){
      x <- lapply(x, FUN = function(x){unname(unlist(x))})
      gmtlist2file(gmtlist = x,filename = paste(outFile, "General.gmt", sep=""))
    }else{
      if(type == "aracne"){
        x2 <- lapply(x, function(y){names(y$tfmode)})
        x3 <- lapply(x, function(y){names(y$tfmode[y$tfmode > 0])})
        x4 <- lapply(x, function(y){names(y$tfmode[y$tfmode < 0])})
        names(x2) <- names(x3) <- names(x4) <- names(x)
        gmtlist2file(gmtlist = x2, filename = paste(outFile, ".generalAracne.gmt", sep=""))
        if(length(x3) > 0){
          gmtlist2file(gmtlist = x3, filename = paste(outFile, ".generalPos.gmt", sep=""))
        }
        if(length(x4) > 0){
          gmtlist2file(gmtlist = x4,filename = paste(outFile, ".generalNeg.gmt", sep=""))
        }
        netReg <- list()
        for(i in 1:length(x2)){
          netReg[[i]] <- list(pos=x3[[i]], neg=x4[[i]])
        }
        names(netReg) <- names(x2)
        return(netReg)
      }else{
    xPos <- xNeg <- list()
    for(i in 1:length(x)){
      xPos[[i]] <- x[[i]]$pos
      xNeg[[i]] <- x[[i]]$neg
    }
    names(xPos) <- names(xNeg) <- names(x)
    if(length(grep("dummy", xPos)) > 0){
      xPos <- xPos[-1*grep("dummy", xPos)]
    }
    if(length(xPos) > 0){
      gmtlist2file(gmtlist = xPos, filename = paste(outFile, "Pos.gmt", sep=""))
    }
    if(length(grep("dummy", xNeg)) > 0){
      xNeg <- xNeg[-1*grep("dummy", xNeg)]
    }
    if(length(xNeg) > 0){
      gmtlist2file(gmtlist = xNeg,filename = paste(outFile, "Neg.gmt", sep=""))
    }
    }
    }
  }
}

extractNES <- function(x, slotN, threshold=0.01, convertNames=F, ora=T){
  allResults <- list()
  allPaths <- c()
  for(i in 1:length(x)){
    pos <- x[[i]][[slotN]]@result$qvalue < threshold
    if(ora){
      finalScore <- -log10(x[[i]][[slotN]]@result$qvalue)[pos]
    }else{
      finalScore <- x[[i]][[slotN]]@result$NES[pos]
    }
    names(finalScore) <- rownames(x[[i]][[slotN]]@result)[pos]
    allResults[[i]] <- finalScore
    allPaths <- union(allPaths, names(finalScore))
  }
  names(allResults) <- names(x)
  M <- matrix(0, nrow = length(x), ncol = length(allPaths))
  rownames(M) <- names(x)
  colnames(M) <- allPaths
  for(i in 1:length(allResults)){
    M[names(allResults)[i], names(allResults[[i]])] <- allResults[[i]] 
  }
  if(convertNames){
    load("startData/allConvertedNames.rda")
    mySymbols <- convertedGenes$symbol
    names(mySymbols)<- convertedGenes$query
    originalNames <- colnames(M)
    colnames(M) <- mySymbols[colnames(M)]
    res <- list(mat=M, originalNames=originalNames)
  }else{
    res <- M
  }
  return(res)
}

calculateNetworkMetrics <- function(x, tfNames=NA, convert=F, outFile="graph"){
  
  ego.effective.size <- function(g, ego, ...) {
    egonet <- induced.subgraph(g, neighbors(g, ego, ...))
    n <- vcount(egonet)
    t <- ecount(egonet)
    return(n - (2 * t) / n)
  }
  
  effective.size <- function(g, ego=NULL, ...) {
    if(!is.null(ego)) {
      return(ego.effective.size(g, ego, ...))
    }
    return (sapply(V(g), function(x) {ego.effective.size(g,x, ...)}))
  }

  tfNumber <- function(x, tfNames){
    n <- length(x[tfNames], )
  }
  
  message("adjusting matrix")
  x <- x[, apply(x, 2, sum) != 0]
  if(dim(x)[1] != dim(x)[2]){
    dif <- setdiff(colnames(x), rownames(x))
    newMat <- matrix(0, nrow = length(dif), ncol=ncol(x))
    rownames(newMat) <- dif
    colnames(newMat) <- colnames(x)
    x <- rbind(x, newMat)
    rm(newMat)
  }else{
    dif <- rownames(x)[!rownames(x) %in% tfNames]
    x[dif, ] <- 0
  }
  common <- intersect(rownames(x), colnames(x))
  x <- x[common, common]
  
  message("preparing graph")
  y <- graph_from_adjacency_matrix(adjmatrix = x, mode = "directed", weighted = T, diag = T, add.rownames = T)
  V(y)$originalNames <- V(y)$name
  x <- x[V(y)$originalNames %in% tfNames, ]
  # Clustering
  com <- walktrap.community(y)
  V(y)$membAll        <- com$membership
  modules <- com$membership
  names(modules) <- V(y)$name
  V(y)$newNames <- V(y)$name
  if(convert){
    load("startData/convertedGenesFromRawCountsFromMygene.rda")
    V(y)$convertedNames <- ensToComplete[V(y)$originalNames]
    ensToComplete <- ensToComplete[names(ensToComplete) %in% tfNames]
    V(y)$newNames       <- sapply(ensToComplete[V(y)$name], function(x){unlist(strsplit(x, "_"))[2]})
  }else{
    V(y)$newNames       <- ifelse(V(y)$name %in% tfNames, unlist(strsplit(V(y)$name, "_"))[2], NA)
  }
  
  message("calculating metrics")
  
  # vertex and edges metrics
  V(y)$degree      <- igraph::degree(y, mode = "total")
  V(y)$indegree    <- igraph::degree(y, mode = "in")
  V(y)$outdegree   <- igraph::degree(y, mode = "out")
  V(y)$betweenness <- betweenness(y, directed = T, normalized = T)
  V(y)$evcent      <- evcent(y, directed = T)$vector
  V(y)$closeness   <- closeness(y, normalized = T)
  # V(y)$flowbet     <- sna::flowbet(as.matrix(get.adjacency(y, attr="weight")))
  E(y)$betweenness <- edge.betweenness(y, directed = T)
  
  # Local position
  V(y)$effsize     <- effective.size(y, mode = "all")
  V(y)$constraint  <- constraint(y)
  
  # Whole network
  set.graph.attribute(y, "cumulativeDegreeDistribution", degree.distribution(graph = y, cumulative = T))
  set.graph.attribute(y, "density", graph.density(y))
  set.graph.attribute(y, "avgpathlength", average.path.length(y))
  # set.graph.attribute(y, "modularity", modularity(com))
  set.graph.attribute(y, "betcentralization", centralization.betweenness(y)$centralization)
  set.graph.attribute(y, "degcentralization", centralization.degree(y, mode = "total")$centralization)
  set.graph.attribute(y, "size", vcount(y))
  set.graph.attribute(y, "edgecount", ecount(y))
  
  message("writing output")
  # z <- networktools::impact(input = x, nodes = "ENSG00000001167")
  write.graph(graph = y, file = paste(outFile, ".complete.gml", sep=""), format = "gml")
  z <- induced_subgraph(graph = y, vids = which(!is.na(V(y)$newNames)))
  com <- walktrap.community(z)
  V(z)$membTf        <- com$membership
  x <- x[V(z)$originalNames, ]
  V(z)$tfNumber <- sapply(1:nrow(x), function(i){sum(x[i, colnames(x) %in% tfNames] != 0)})
  rm(x)
  # plot(com, z)
  # plot_dendrogram(com)
  plot(z, vertex.label=NA, vertex.size=rangeNormalisation(V(z)$tfNumber,
   c(1,10)), vertex.color=membership(com), edge.width=1)
  # layout <- layout.reingold.tilford(z, circular=T)
  # plot.igraph(y, vertex.label=NA, vertex.size=rangeNormalisation(V(y)$outdegree, c(1,10)), 
  #             edge.color="black",edge.width=E(y)$weight/0.00005)
  write.graph(graph = z, file = paste(outFile, ".onlyTf.gml", sep=""), format = "gml")
  save(modules, y, z, file = paste(outFile, ".graphData.rda", sep=""))
  return(list(complete=y, tfOnly=z, modules=modules))
  }

extractGraphAttributes <- function(y, z, DEGs, threshold = 0.05, foldChange=1.5, finalPlot=F){
  
  pos <- which(DEGs$padj < threshold & abs(DEGs$log2FoldChange) > foldChange)
  rankedList <- -log10(DEGs[pos, "padj"]) # we should add pos again
  sign <- ifelse(DEGs[pos, "log2FoldChange"] > 0, "+", "-")
  names(rankedList) <- rownames(DEGs)[pos]
  rankedPos <- rankedList[sign=="+"]
  rankedNeg <- rankedList[sign=="-"]
  rankedList <- rankedList[!is.na(rankedList)]
  rankedPos <- rankedPos[!is.na(rankedPos)]
  rankedNeg <- rankedNeg[!is.na(rankedNeg)]
  
  allAttributes <- lapply(list.vertex.attributes(y),function(x) get.vertex.attribute(y,x))
  names(allAttributes) <- list.vertex.attributes(y)
  allAttributes2 <- lapply(list.vertex.attributes(z),function(x) get.vertex.attribute(z,x))
  names(allAttributes2) <- list.vertex.attributes(z)
  common <- intersect(V(y)$originalNames,names(rankedList))
  y2 <- induced_subgraph(graph = y, vids = V(y)[V(y)$originalNames %in% names(rankedList)])
  allAttributes3 <- lapply(list.vertex.attributes(y2),function(x) get.vertex.attribute(y2,x))
  names(allAttributes3) <- list.vertex.attributes(y2)
  
  if(finalPlot){
    V(y2)$name <- V(y2)$originalNames
    present <- intersect(common, names(rankedList))
    V(y2)[present]$rankedList <- rankedList[present]
    V(y2)[present]$color <- ifelse(present %in% names(rankedPos), yes = 1, no = 2)
    plot(y2, vertex.label=NA, vertex.size=rangeNormalisation(V(y2)$rankedList,c(1,10)), vertex.color=V(y2)$color, edge.width=1)
  }
  
  metrics <- rbind(unlist(lapply(allAttributes[6:13], FUN = mean, na.rm=T)),
                     unlist(lapply(allAttributes2[6:13], FUN = mean, na.rm=T)),
                     unlist(lapply(allAttributes3[6:13], FUN = mean, na.rm=T)))
  rownames(metrics) <- c("all", "tfs", "degs")
  return(metrics)
}

rangeNormalisation <- function(x, range){ 
  normData <- (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
  diff <- (range[2] - range[1])
  normData2 <- (normData * diff) + range[1] 
  normData2[is.na(normData2)] <- range[1]
  return(normData2)
} 

rangeNormalisationMat <- function(x, range){
  a <- rangeNormalisation(as.vector(x), range=range)
  b <- matrix(a, nrow = nrow(x), ncol=ncol(x), byrow = F)
  rownames(b) <- rownames(x)
  colnames(b) <- colnames(x)
  return(b)
}

gaussianKernel <- function(x1, x2, alpha=1) {
  exp(- alpha * norm(as.matrix(x1-x2), type="F"))
}

makeSimilarity <- function(my.data, similarity) {
  N <- nrow(my.data)
  S <- matrix(rep(NA,N^2), ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      S[i,j] <- similarity(my.data[i,], my.data[j,])
    }
  }
  S
}

makeAffinity <- function(S, n.neighboors=2) {
  N <- length(S[,1])
  
  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  A  
}

masterRegulatorAnalysis <- function(adj=NULL, z, treatments= infoTable$GSid[which(infoTable$CRPC == T)], tf = newTFs,  
                                    controls = infoTable$GSid[which(infoTable$CRPC == F)], adjfile="", outFile, ncores=32){
  if(!is.null(adj)){
  prepareDfAracne <- function(adj, outFile){
    
    make3Col <- function(x1, tf){
      x2 <- colnames(adj)[which(adj[x1, ] != 0)]
      if(length(x2) > 0){
        x3 <- adj[x1, which(adj[x1, ] != 0)]
        x1 <- unname(rownames(adj)[x1])
        a <- cbind(rep(x1, length(x2)), x2, x3)
      }else{
        a <- c()
      }
      return(a)
    }
    
    l <- c()
    l <- do.call("rbind", sapply(X = 1:nrow(adj), FUN = function(x){make3Col(x, tf)}))
    if(!is.na(tf[1])){
      l <- l[l[, 1] %in% tf, ]
    }
    write.table(x = l, file = outFile, sep = "\t", quote = F, row.names = F, col.names = F)
  }  
  }
  
  message("regulon object preparation")
  if(!file.exists(adjfile)){
  prepareDfAracne(adj = adj, outFile = adjfile)
  }
  regul <- viper::aracne2regulon(afile = adjfile, verbose = F, eset = z)
  message("msviper")
  signature <- ssmarina::rowTtest(x = z[, treatments], y = z[, controls])
  signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
  nullmodel <- viper::ttestNull(x = z[, treatments], y = z[, controls], per = 1000, repos = TRUE, verbose = FALSE, cores=ncores)
  msviper <- viper::msviper(signature, regul, nullmodel, verbose = F, pleiotropy = T, cores=ncores)
  msviper <- viper::ledge(msviper)
  combinatorial <- viper::msviperCombinatorial(msviper, regulators = 25, verbose = FALSE, cores=ncores) # we can set a pvalue threshold
  shadowViper <- viper::shadow(msviper, regulators = 25, verbose = F, cores=ncores) # we can also set a pvalue treshold for this
  # sinergyViper <- viper::msviperSynergy(combinatorial, verbose = F)
  message("viper")
  vpres <- viper::viper(z, regul, dnull = nullmodel, verbose = FALSE, cores=ncores)
  message("marina")
  msmarina <- ssmarina::marina(signature, regul, nullmodel, verbose = F)
  msssmarina <- ssmarina::ssmarina(eset = z, regulon = regul, verbose = F)
  message("done")
  results <- list(reguolons=regul, signature=signature, nullmodel=nullmodel, ssviper=vpres, viperSig=msviper, combinatorialViper=combinatorial,
                  shadow=shadowViper, marinaSig=msmarina, ssmarina=msssmarina)
  save(results, file = paste(outFile, ".mra.rda", sep=""))
  return(results)
} 

prepareIntVector <- function(x, n, orientation=F){
  if(orientation){
    d <- c(paste(names(x)[n], x[[n]]$pos, sep="|+|"), paste(names(x)[n], x[[n]]$neg, sep="|-|"))
  }else{
    d <- paste(names(x)[n], x[[n]], sep="->")
  }
  return(d)
}

convertList <- function(x, n, vect){
  
  b <- vect[unlist(x[[n]]$pos)]
  d <- grep("more", b)
  if(length(d) > 0){
    b <- b[-1 * d]
  }
  b <- b[!is.na(b)]
  
  e <- vect[unlist(x[[n]]$neg)]
  d <- grep("more", e)
  if(length(d) > 0){
    e <- e[-1 * d]
  }
  e <- e[!is.na(e)]
  
  f <- list(pos=unname(b), neg=unname(e))
  return(f)
  }

adjListToAdjMat <- function(x){
  allGenes <- union(x$V1, x$V2)
  M <- matrix(0, length(allGenes), length(allGenes))
  rownames(M) <- colnames(M) <- allGenes
  index1 <- sapply(x$V1, function(y){which(rownames(M) %in% y)})
  index2 <- sapply(x$V2, function(y){which(rownames(M) %in% y)})
  M[cbind(index1, index2)] <- x$V3
  return(M)
  }

prepareGRNforEB <- function(x){
  df <- c()
  for(i in 1:length(x)){
    a <- cbind(rep(names(x)[i], length(x[[i]]$pos)), x[[i]]$pos, rep("+", length(x[[i]]$pos)))
    b <- cbind(rep(names(x)[i], length(x[[i]]$neg)), x[[i]]$neg, rep("-", length(x[[i]]$neg)))
    df <- rbind(df, a, b)
  }
  colnames(df) <- c("FROM", "TO", "TYPE")
  return(df)
}

topologyAnalysis2 <- function(h, netReg, degOutput=NA, outFile="prova.rda", treatments,
                              controls, myAlpha=0.05, convert=T, singleSample=F, ncores=10, 
                              ssDegs){
  
  hsa.grn <- prepareGRNforEB(netReg)
  hsa.gs <- lapply(netReg, unlist)
  library(EnrichmentBrowser)
  if(singleSample){
    library(foreach)
    library(doSNOW)
    C <- ncores
    cl <- makeCluster(C, type = "SOCK")
    registerDoSNOW(cl)
    start.time <- Sys.time()
    pb <- txtProgressBar(min = 0, max = length(treatments), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    nbea.res <- foreach(.errorhandling = "pass", i= 1:length(treatments), .packages = c("EnrichmentBrowser", "SummarizedExperiment"), .options.snow = opts)%dopar%{
      
      original <- h
      degOutput <- ssDegs[[treatments[i]]]
      colnames(degOutput) <- c("log2FoldChange", "padj")
      common <- intersect(rownames(degOutput), rownames(original))
      newZ <- original[common, c(treatments[i], controls[!controls %in% treatments[i]])]
      if(convert){
        load("D://mcangiano/Glasgow/RNAseq_orthografts/startData/convertedGenesFromRawCountsFromMygene.rda")
        a <- ensToComplete[!is.na(ensToComplete)]
        a <- a[-1*length(grep("more", a))]
        a <- a[!duplicated(a)]
        degOutput <- degOutput[rownames(degOutput) %in% names(a), ]
        rownames(degOutput) <- a[rownames(degOutput)]
        newZ <- newZ[rownames(newZ) %in% names(a), ]
        rownames(newZ) <- unname(a[rownames(newZ)])
        common <- intersect(rownames(newZ), rownames(degOutput))
      }
      colData <- DataFrame(GROUP= ifelse(colnames(newZ) %in% treatments[i], 1, 0))
      rowData <- DataFrame(FC = degOutput[common, "log2FoldChange"], ADJ.PVAL= degOutput[common, "padj"])
      zSum <- SummarizedExperiment(assays=newZ[common, ], rowData=rowData, colData=colData, metadata=list(dataType="rseq"))
      tmp <-  tryCatch(expr = nbea(method="ggea", se = zSum, gs=hsa.gs, grn=hsa.grn, alpha = myAlpha, perm = 1000, padj.method = "fdr"), error=function(e){return(e)})
      return(tmp)
    }
    close(pb)
    stopCluster(cl)
    endTime <- Sys.time()
    names(nbea.res) <- treatments
  }else{
    if(convert){
      message("converting names")
      rownames(degOutput) <- sapply(rownames(degOutput), function(x){unlist(strsplit(x, "\\."))[1]})
      load("D://mcangiano/Glasgow/RNAseq_orthografts/startData/convertedGenesFromRawCountsFromMygene.rda")
      a <- ensToComplete[!is.na(ensToComplete)]
      a <- a[-1*length(grep("more", a))]
      a <- a[!duplicated(a)]
      degOutput <- degOutput[rownames(degOutput) %in% names(a), ]
      rownames(degOutput) <- a[rownames(degOutput)]
    }
    common <- intersect(rownames(degOutput), rownames(h))
    message(length(common))
    h <- h[common, c(treatments, controls)]
    colData <- DataFrame(GROUP= ifelse(colnames(h) %in% treatments, 1, 0))
    rowData <- DataFrame(FC = degOutput[common, "log2FoldChange"], ADJ.PVAL= degOutput[common, "padj"])
    zSum <- SummarizedExperiment(assays=h, rowData=rowData, colData=colData, metadata=list(dataType="rseq"))
    message("performing enrichment")
    nbea.res <- tryCatch(expr = nbea(method="ggea", se = zSum, gs=hsa.gs, grn=hsa.grn, alpha = myAlpha, padj.method = "fdr"), 
                         error=function(e){return(e)})
    message("done")
  }
  save(nbea.res, file=outFile)
  return(nbea.res)
}

makeSingleSampleDegs2 <- function(normCounts, normals, all=T){
  
  normCounts <- as.matrix(normCounts)
  if(all){
    notNormals <- colnames(normCounts)
  }else{
    notNormals <- colnames(normCounts)[!colnames(normCounts) %in% normals] 
  }
  sdNormals <- apply(normCounts[, normals], 1, sd, na.rm=T)
  sdNormals[sdNormals == 0] <- min(sdNormals[sdNormals != 0])
  meanNormals <- apply(normCounts[, normals], 1, mean, na.rm=T)
  meanNormals[meanNormals == 0] <- min(meanNormals[meanNormals != 0])
  fcMat <- log2(normCounts[, notNormals] / meanNormals)
  fcMat[is.infinite(fcMat)] <- 0
  pvMat <- (normCounts[, notNormals] - meanNormals) / sdNormals
  
  pvMat <- apply(pvMat, 2, FUN = function(x){1 - pnorm(q = abs(x), mean = 0, sd = 1)})
  pvMat <- apply(pvMat, 2, p.adjust, "fdr")
  allRes <- list()
  for(i in 1:ncol(fcMat)){
    df <- cbind(fcMat[, i], pvMat[, i])
    colnames(df) <- c("log2FoldChange", "qvalue")
    allRes[[i]] <- df
  }
  names(allRes) <- notNormals
  return(allRes)
} # edited on 29-10-20

extractNormScore <- function(x){
  library(plyr)
  myMat <- do.call("rbind.fill.matrix", lapply(x, function(y){z <- y$res.tbl$NORM.SCORE; names(z) <- y$res.tbl$GENE.SET; return(t(as.data.frame(z)))}))
  rownames(myMat) <- names(x)
  return(myMat)
} 

pvalueMat <- function(x, threshold=0.05, adjusted=T){
  library(plyr)
  if(adjusted){
    myMat <- do.call("rbind.fill.matrix", lapply(x, function(y){z <- y$res.tbl$ADJ.PVAL; names(z) <- y$res.tbl$GENE.SET; return(t(as.data.frame(z)))}))
  }else{
    myMat <- do.call("rbind.fill.matrix", lapply(x, function(y){z <- y$res.tbl$PVAL; names(z) <- y$res.tbl$GENE.SET; return(t(as.data.frame(z)))}))
  }
  rownames(myMat) <- names(x)
  myMat[myMat >= threshold & !is.na(myMat)] <- 100
  myMat[myMat != 100 & !is.na(myMat)] <- 1
  myMat[myMat == 100 | is.na(myMat)] <- 0
  attr(myMat, "adjusted pvalue") <- adjusted
  return(myMat)
}

convert.z.score<-function(z, one.sided=NULL) {
  if(is.null(one.sided)) {
    pval = pnorm(-abs(z));
    pval = 2 * pval
  } else if(one.sided=="-") {
    pval = pnorm(z);
  } else {
    pval = pnorm(-z);                                                                                 
  }
  return(pval);
}  

perPASanalysis <- function(z, graph, group=  list(ctrl=infoTable2$Sample_ID[which(infoTable2$type == "Solid_Tissue_Normal")],
                                                  treat=infoTable2$Sample_ID[which(infoTable2$type == "Primary_Tumor")]),  ncores, outFile, tfList){
  
  message("praparing data")
  startData <- lcy.tORz.score(data = z, group = group, pvalue=T)
  tfVerteces <- V(graph)[V(graph)$name %in% tfList]
  library(igraph)
  singlePaths <- list()
  pb <- txtProgressBar(min = 0, max = length(tfVerteces), style=3)
  for(i in 1:length(tfVerteces)){
    myVertex <- tfVerteces[i]
    targets <- head_of(graph = fullGraph, E(graph)[from(myVertex)])
    if(length(targets) > 0){
      singlePaths[[i]] <- induced.subgraph(graph = fullGraph, vids = c(myVertex, targets))
    }else{
      singlePaths[[i]] <- NA
    }
    setTxtProgressBar(pb, i)
  }
  names(singlePaths) <- names(tfVerteces)
  singlePaths <- singlePaths[unlist(lapply(singlePaths, function(x){!is.na(x[[1]])}))]
  close(pb)
  
  # 
  # myFunction <- function(i, tfVerteces, graph){
  #   library(igraph)
  #   myVertex <- tfVerteces[i]
  #   targets <- igraph::head_of(graph = graph, E(graph)[from(myVertex)])
  #   if(length(targets) > 0){
  #     singlePath <- igraph::induced.subgraph(graph = graph, vids = c(myVertex, targets))
  #     results <- lcy.perpas(g=singlePath, data = startData, group = group)
  #   }else{
  #     results <- rep(NA, ncol(z))
  #   }
  #   return(results)
  # }
  
  # perpas.res <- parSapply(cl = cl, X = 1:2, FUN = myFunction, tfVerteces=tfVerteces, graph=fullGraph)
  # stopCluster(cl)
  
  message("starting analysis")
  library(foreach)
  library(doSNOW)
  C <- ncores
  cl <- makeCluster(C, type = "SOCK")
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 0, max = length(singlePaths), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  perpas.res <- foreach(.errorhandling = "pass", i= 1:length(singlePaths), .packages = c("PerPAS"), .options.snow = opts)%dopar%{
  # library(igraph)
  # library(PerPAS)
  # perpas.res <- c()
  # for(i in 1:length(tfVerteces)){
  results <- lcy.perpas(g=singlePaths[[i]], data = startData, group = group)
  # perpas.res <- rbind(perpas.res, results)
  # setTxtProgressBar(pb, i)
  return(results)
  }
  close(pb)
  stopCluster(cl)
  message("done")
  names(perpas.res) <- names(singlePaths)
  finalMatrix <- do.call("rbind", perpas.res)
  save(perpas.res, file=outFile)
  return(perpas.res)
}

ggeaSynergies <- function(h, netReg, degOutput=res, outFile="prova.rda", treatments= infoTable$GSid[which(infoTable$CRPC == T)],
                          controls = infoTable$GSid[which(infoTable$CRPC == F)],
                          myAlpha=0.05, convert=T, singleSample=F, ncores=30, minIntersection=15){
  allCombs <- combn(x = names(netReg), m = 2)
  library(doMC)
  cl <- snow::makeCluster(ncores, type="SOCK")
  clusterExport(varlist = c("netReg", "allCombs"), cl = cl)
  allIntersections <- parSapply(cl = cl, X = 1:ncol(allCombs), FUN = function(x){length(intersect(netReg[[allCombs[1, x]]]$pos, 
                      netReg[[allCombs[2, x]]]$pos)) + length(intersect(netReg[[allCombs[1, x]]]$neg, 
                                                                        netReg[[allCombs[2, x]]]$neg))})
  stopCluster(cl)
  myCombs <- allCombs[, allIntersections >= minIntersection]
  allRes <- list()
  pb <- txtProgressBar(min = 1, max = ncol(myCombs), style = 3)
  for(i in 1:ncol(myCombs)){
    miniReg <- list(netReg[[myCombs[1, i]]], netReg[[myCombs[2, i]]], prepareIntersection(netReg[myCombs[1, i]], netReg[myCombs[2, i]]))
    names(miniReg) <- c(myCombs[1, i], myCombs[2, i], paste(myCombs[1, i], myCombs[2, i], sep="_"))
    allRes[[i]] <- tryCatch(expr = {topologyAnalysis(h, miniReg, degOutput = res)}, error=function(e){NA})
    setTxtProgressBar(pb, i)
    }
  close(pb)
  names(allRes) <- apply(myCombs, 2, paste, collapse="_")
  allRes <- allRes[!is.na(allRes)]
  # there is no synergy -  to be completed
  }

extractOra <- function(x, minOverlap=2, threshold=0.05){
  if(length(which(unlist(lapply(x, is.null)))) > 0){
    x <- x[ -1 * which(unlist(lapply(x, is.null)))]
  }
  allRes <- c()
  myNames <- c()
  for(i in 1:length(x)){
    res <- x[[i]]@result[x[[i]]@result$Count >= minOverlap & x[[i]]@result$qvalue < threshold, ][1, ]
    if(nrow(res) > 0){
      myNames <- c(myNames, rep(names(x)[i], nrow(res)))
    }else{
      myNames <- c(myNames, names(x)[i])
    }
    allRes <- rbind(allRes, res)
  }
  allRes <- as.data.frame(allRes)
  rownames(allRes) <- myNames
  # allGenesets <- lapply(x, function(y){y@geneSets})
  # allGenesets <- Reduce(union, allGenesets)
  # allRes$originalNumbers <- lapply(x[[1]]@geneSets[allRes$ID], length)
  return(allRes)
}

##########################################################################################
