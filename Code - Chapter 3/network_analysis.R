
#########################################################################################

setwd("E://mcangiano/Glasgow/RNAseq_orthografts/")
source("functionsAndLibraries.R")

######################################################################################### annotation

library(xlsx)
infoTable <- read.xlsx2("/mnt/data/mcangiano/7060_Glasgow/RNAseq/MSR8_SampleSheet_170306.xls", sheetIndex = 1, as.data.frame = T, startRow = 21, 
                        header = T, stringsAsFactors=F, colIndex = c(1:5))
infoTable$GSid <- colnames(rawCounts)
infoTable$cellType <- sapply(X = infoTable$Sample_ID, FUN = function(x){unlist(strsplit(x, "-"))[2]})
infoTable$CRPC <- rep(c(rep(F,3), rep(T, 3)), 3)
rownames(infoTable) <- infoTable$GSid
save(infoTable, file = "startData/infoTable.rda")

######################################################################################## normCounts

fl <- list.files(path = "../fpkm/", full.names = T)
fl2 <- list.files(path = "../fpkm/")
library(plyr)
normCounts <- c()
pb <- txtProgressBar(min = 0, max = length(fl), style = 3)
for(i in 1:length(fl)){
     a <- read.table(fl[i], stringsAsFactors = F, header = T)[, c(4,5,7,10)]
     a$full_name <- paste(a$gene_id, a$gene_short_name, a$locus, sep="_")
     b <- tapply(a$FPKM, a$full_name, sum)
     if(i == 1){
       normCounts <- cbind(normCounts, b)
     }else{
     normCounts <- cbind(normCounts, b[rownames(normCounts)])
     }
setTxtProgressBar(pb, i)     
}
close(pb)
colnames(normCounts) <- sapply(fl2, function(x){unlist(strsplit(x, ".fpkm"))[1]})
save(normCounts, file = "startData/normCounts.rda")

#########################################################################################  RTN 

load("/opt/shared/easystore2/MarioCangiano/generalData/tf_list_20181807.rda")
newTFs <- findFullNames(tfList, rownames(normCounts))
save(newTFs, file = "startData/tfListFrancesca.rda")

load("startData/normCounts.rda")
load("startData/tfListFrancesca.rda")
load("startData/infoTable.rda")
newTFs <- newTFs[!is.na(newTFs)] 
newTFs <- newTFs[!duplicated(newTFs)]
filtered <- normCounts[apply(normCounts, 1, max) > 0, ]
standardised <- t(apply(X = normCounts, MARGIN = 1, FUN = rangeNormalisation, range=c(0,1))) # should not be needed
standardised <- standardised[apply(standardised, 1, max) > 0, ]

rtni <- tni.constructor(expData=standardised, regulatoryElements=newTFs) # it doesn't work with 'filtered'
options(cluster=makeCluster(20, "SOCK"))
rtni <- tni.permutation(rtni, verbose = T, nPermutations = 900, parChunks = 30)
rtni <- tni.bootstrap(rtni) 
save(rtni, "RTN/standardisedCounts_pearson/standardisedNormCountsBootstrap.rda")
rtni <- tni.dpi.filter(rtni)
colnames(factorDf) <- infoTable$GSid
rtni@gexp <- rbind(rtni@gexp, 
                   ifelse(infoTable[colnames(rtni@gexp), "CRPC"] == T, yes = 1, no = 0))
rownames(rtni@gexp)[nrow(rtni@gexp)] <- "CRPC_status"
rtni <- tni.conditional(minRegulonSize = 3, medianEffect = T, verbose = T, # changed parameters, saved file is with default
                        sampling = 15, object = rtni, modulators = "CRPC_status", tfs = rtni@regulatoryElements) 
save(rtni, file="RTN/standardisedCounts_pearson/standardisedNormCountsConditional.rda")

reg <- tni.get(rtni, what= "regulons")
usedTfs <- rtni@regulatoryElements
g <- tni.graph(rtni) 
RTNsingleNet <- t(tni.get(rtni, what= "tnet"))
standRTN <- t(apply(RTNsingleNet, MARGIN = 1, rangeNormalisation, range=c(0,1)))
netReg <- prepareFilteredRegulons(out = standRTN, z = standardised, tf_genes = newTFs)
length(intersect(netReg$`ENSG00000156150_ALX3_1:110602615-110613322`$pos, reg$`ENSG00000156150_ALX3_1:110602615-110613322`)) # 291
save(RTNsingleNet, reg, netReg, g, file = "RTN/standardisedCounts_pearson/RTNsingleReg.rda")
makeGmt(x = reg, outFile = "RTN/standardisedCounts_pearson/RTNsingleReg", type="RTN")

######################################################################## filtering by shadowing and modulations

load("E:/mcangiano/Glasgow/RNAseq_orthografts/RTN/standardisedCounts_pearson/standardisedNormCountsConditional.rda")

adj <- rtni@results$tn.dpi
prepareDfAracne(adj = adj, outFile = adjfile)
regul <- viper::aracne2regulon(afile = adjfile, verbose = F, eset = z)
message("msviper")
signature <- ssmarina::rowTtest(x = z[, treatments], y = z[, controls])
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
nullmodel <- viper::ttestNull(x = z[, treatments], y = z[, controls], per = 1000, repos = TRUE, verbose = FALSE, cores=ncores)
msviper <- viper::msviper(signature, regul, nullmodel, verbose = F, pleiotropy = T, cores=ncores)
msviper <- viper::ledge(msviper)
combinatorial <- viper::msviperCombinatorial(msviper, regulators = 25, verbose = FALSE, cores=ncores) 
shadowViper <- viper::shadow(msviper, regulators = 25, verbose = F, cores=ncores) 
save(shadowViper$regulons, file = "RTN/standardisedCounts_pearson/RTNsingleRegAracne.shadowFiltered.rda")

cdt <- tni.get(rtni, what = "cdtrev")[[1]] # empty
modulatorsEffect <- rtni@results$conditional$effect[names(rtni@results$conditional$effect) %in% names(netRegShadowFiltered)]
a <- lapply(modulatorsEffect, function(x){table(as.vector(as.matrix(x[, -1])))})
a2 <- sum(unlist(lapply(a, function(x){is.na(x["0"])}))) # 0
b <- summary(unlist(lapply(a, length))) # only 0
modulatorsCounts <- rtni@results$conditional$count[names(rtni@results$conditional$effect) %in% names(netRegShadowFiltered)]
a <- do.call(what = "rbind", modulatorsCounts)
rownames(a) <- 1:nrow(a)

####################################################################################### network topology

library(igraph)

load("RTN/standardisedCounts_pearson/RTNsingleRegAracne.shadowFiltered.rda")
netReg <- lapply(netRegShadowFiltered, function(x){unname(unlist(x))})
summary(unlist(lapply(netReg, length)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   11.00   20.00   24.91   34.25  121.00 
myDf <- prepareGRNforEB(netRegShadowFiltered)
bigGraph <- graph_from_data_frame(myDf)
myCs <- constraint(graph = bigGraph)
sum(myCs == 1) # 4924
myCl <- closeness(bigGraph, normalized = T)
sum(myCl == min(myCl)) # 13029
length(intersect(names(myCl)[myCl==min(myCl)], names(myCs)[myCs == max(myCs)])) # 4924
myBw <- betweenness(graph = bigGraph, normalized = T)
hist(log10(myBw))
sum(myBw == max(myBw))
which.max(myBw) # ENSG00000081692_JMJD4_1:227918125-227923112
myDg <- degree(graph = bigGraph, mode = "in")
tt <- table(unlist(netReg))
head(sort(tt, decreasing = T), 10)
# ENSG00000106245_BUD31_7:99006263-99017239  ENSG00000106397_PLOD3_7:100849257-100861701 ENSG00000184860_SDR42E1_16:82031220-82045093 
# 10                                           10                                           10 
# ENSG00000204379_XAGE1A_X:52255218-52260363   ENSG00000019505_SYT13_11:45261851-45307870    ENSG00000044446_PHKA2_X:18910417-19002716 
# 10                                            9                                            9 
# ENSG00000108219_TSPAN14_10:82213921-82292879  ENSG00000111231_GPN3_12:110890288-110907073    ENSG00000117305_HMGCL_1:24128374-24165110 
# 9                                            9                                            9 
# ENSG00000118816_CCNI_4:77968310-77997158 
# 9 
fullDetails <- cbind(myCs, myCl, myBw, myDg)
colnames(fullDetails) <- c("constraint", "closeness", "betweenness", "degree")
library(xlsx)
write.xlsx(fullDetails, file = "/mnt/data/mcangiano/otherProjects/Cathal/survivalAnalysis/table_paper/networkStatistics.xlsx", sheetName = "statistics")

allIntersections <- sapply(X = 1:length(netReg), function(x){sapply(1:length(netReg), function(y){
  length(intersect(netReg[[x]], netReg[[y]]))/length(unique(c(netReg[[x]], netReg[[y]])))})})
rownames(allIntersections) <- colnames(allIntersections) <- names(netReg)
diag(allIntersections) <- 0
sum(rowSums(allIntersections) == 0) # 19
sum(colSums(allIntersections) == 0) # 19
allIntersections2 <- allIntersections
allIntersections2[allIntersections <= 0.1] <- 0

netRegGraph <- graph_from_adjacency_matrix(adjmatrix = allIntersections2, mode = "undirected", weighted = T, diag = F)
subgraph <- induced_subgraph(graph = netRegGraph, vids = names(degree(netRegGraph))[degree(netRegGraph) > 0])

regressionsRank <- read.csv("E://mcangiano/survivalAnalysis/table_paper/coxInternalDatasets.csv", header = T, stringsAsFactors = F)
myRank <- rowMeans(regressionsRank[, 2:3])
names(myRank) <- regressionsRank$X
myRank <- myRank[!is.na(myRank)]
library(network)
library(plotrix)
myRank2 <- color.scale(x = myRank, color.spec="rgb", extremes = c("red", "yellow"))
names(myRank2) <- names(myRank)
verCol <- myRank2[V(subgraph)$name]
verCol[is.na(verCol)] <- "grey"

plot(subgraph, vertex.label= sapply(V(subgraph)$name, function(x){ifelse(x %in% names(myRank2), unlist(strsplit(x, "_"))[2], NA)}), 
     vertex.color= verCol, 
     edge.color= color.scale(x = E(subgraph)$weight, color.spec="rgb", extremes = c("royalblue4", "green"), alpha = 0.5), 
     edge.width=1, vertex.label.cex = 0.9, vertex.label.font=2, vertex.label.color= "black",  
     vertex.size=log2(unlist(lapply(netReg[V(subgraph)$name], length))), layout=layout_with_dh(subgraph)) # layout_with_dh

plot(myRank, col=myRank2, ylab = "P-VALUE", font.lab=2, axes = F, ylim=c(0,1))
for(i in 1:length(myRank)){
  if(i == 1){
    rect(0, -0.001, 70, myRank[i], col = myRank2[i])
  }else{
    rect(0, myRank[i-1], 70, myRank[i], col = myRank2[i])
  }
}
axis(2, at = c(0, 0.1, 1), las=2, font.lab=2)

############################################################################### Zoom-in on JMJD6

myGenes <- c("ENSG00000070495_JMJD6_17:74708918-74722866")
netRegAttempt <- netRegShadowFiltered[myGenes] # NULL results as expected
netRegAttempt2 <- c(netRegAttempt, netRegShadowFiltered[unlist(netRegAttempt)])
library(igraph)
myNodes <- unique(unname(unlist(netRegAttempt2)))
myGrn <- prepareGRNforEB(netRegAttempt2)
myGraph <- graph_from_data_frame(d = myGrn, directed = F)
V(myGraph)$degree <- round(log2(degree(myGraph)), digits = 0) + 1
verCol <- myRank2[V(myGraph)$name]
verCol[is.na(verCol)] <- "white"
load("E://mcangiano/Glasgow/RNAseq_orthografts/startData/tfListFrancesca.rda")
newTFs <- newTFs[!is.na(newTFs)] 
newTFs <- newTFs[!duplicated(newTFs)]
V(myGraph)$name <- unname(sapply(V(myGraph)$name, function(x){unlist(strsplit(x, "_"))[2]}))

plot(myGraph, vertex.label=ifelse(V(myGraph)$name %in% names(newTFs), V(myGraph)$name, NA), 
     vertex.size=V(myGraph)$degree, vertex.color=verCol,layout=layout_with_dh(myGraph), 
     edge.width=3, vertex.label.cex = 3, vertex.label.font=2, vertex.label.color="black",
     edge.color=ifelse(E(myGraph)$TYPE == "+", "tomato", "steelblue"))

##################################################################################### percentage of commonly ssDegs

load("E://mcangiano/survivalAnalysis/startData/startData/patientsDegsMat_2.rda")
load("E://mcangiano/survivalAnalysis/startData/startingMatrices.rda")

upEmc <- sapply(rownames(emcMat), function(x){rownames(foldChangeEmc)[foldChangeEmc[, x] > 0 & qvalueEmc[, x] < 0.05]})
upUta <- sapply(rownames(tampereMat), function(x){rownames(foldChangeUta)[foldChangeUta[, x] > 0 & qvalueUta[, x] < 0.05]})
downEmc <- sapply(rownames(emcMat), function(x){rownames(foldChangeEmc)[foldChangeEmc[, x] < 0 & qvalueEmc[, x] < 0.05]})
names(upEmc) <- names(downEmc) <- rownames(emcMat)
downUta <- sapply(rownames(tampereMat), function(x){rownames(foldChangeUta)[foldChangeUta[, x] < 0 & qvalueUta[, x] < 0.05]})
names(upUta) <- names(downUta) <- rownames(tampereMat)
upDist <- unlist(lapply(combn(c(upEmc, upUta), 2, simplify = FALSE), function(x) {
  length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))
upDist2 <- cbind(t(combn(length(c(upEmc, upUta)), 2)), upDist)
upMat <- matrix(data = median(upDist2[, 3]), nrow = length(c(upEmc, upUta)), ncol = length(c(upEmc, upUta)))
upMat[upDist2[, 1:2]] <- upMat[upDist2[, 2:1]] <- upDist2[, 3]
downDist <- unlist(lapply(combn(c(downEmc, downUta), 2, simplify = FALSE), function(x) {
  length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))
downDist2 <- cbind(t(combn(length(c(downEmc, downUta)), 2)), downDist)
downMat <- matrix(data = median(downDist2[, 3]), nrow = length(c(downEmc, downUta)), ncol = length(c(downEmc, downUta)))
downMat[downDist2[, 1:2]] <- downMat[downDist2[, 2:1]] <- downDist2[, 3]
rownames(upMat) <- rownames(downMat) <- colnames(upMat) <- colnames(downMat) <- c(rownames(rotterdamMat), rownames(tampereMat))
library(heatmap3)
# myAnnot <- cbind(DATASET=as.factor(c(rep("EMC", length(upEmc)), rep("UTA", length(upUta)))), 
#                  PROGRESSION=as.factor(c(c("No", "Yes")[metaData2[names(upEmc), "PSAProgAfterRP01YN"] + 1], 
#                                          sampleCoding[names(upUta), "Progression"])))
# rownames(myAnnot) <- rownames(upMat)
# myAnnot2 <- myAnnot[!is.na(myAnnot[, 2]), ]
heatmap3(upMat[rownames(myAnnot), rownames(myAnnot)], scale = "none", ColSideColors = myAnnot, labCol = NA, labRow = NA, main = "Upregulated Jaccard Index")
heatmap3(downMat[rownames(myAnnot), rownames(myAnnot)], scale = "none", ColSideColors = myAnnot, labCol = NA, labRow = NA, main = "Downregulated Jaccard Index")

summary(unlist(lapply(upEmc, length)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 827    1971    2439    2563    2793    7395 
summary(unlist(lapply(upUta, length)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1098    1897    2406    2633    3063    6419 
summary(unlist(lapply(downEmc, length)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    45.0   126.0   175.9   217.0   952.0 
summary(unlist(lapply(downUta, length)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 44.0   149.5   282.0   363.2   450.0  1173.0 
Reduce(f = intersect, x = upEmc) # 0 
length(Reduce(f = union, x = upEmc)) # 30035
Reduce(f = intersect, x = upUta) # 7among which PCAT7 and four lncoding
length(Reduce(f = union, x = upUta)) # 21165
Reduce(f = intersect, x = downEmc) # none
length(Reduce(f = union, x = downEmc)) # 2687
Reduce(f = intersect, x = downUta) # 0 
length(Reduce(f = union, x = downUta)) # 2960

foldChangeIcgc[pvalueIcgc > 0.05 & !is.na(pvalueIcgc)] <- NA
foldChangeUta[pvalueUta > 0.05 & !is.na(pvalueUta)] <- NA
foldChangeEmc[pvalueEmc > 0.05 & !is.na(pvalueEmc)] <- NA
myGenes <- c("ENSG00000105707", "ENSG00000232677", "ENSG00000167107", "ENSG00000178750", "ENSG00000183317", "ENSG00000269190", 
             "ENSG00000156284", "ENSG00000234949", "ENSG00000125247", "ENSG00000260228")
myGenesUta <- foldChangeUta[myGenes, ]
myGenesEmc <- foldChangeEmc[myGenes, ]
library(plyr)
matrix.combined <- rbind.fill.matrix(myGenesUta, myGenesEmc)
matrix.combined2 <- matrix.combined
matrix.combined2[seq(1, nrow(matrix.combined), 2), ] <- matrix.combined[1:10, ]
matrix.combined2[seq(2, nrow(matrix.combined), 2), ] <- matrix.combined[11:20, ]
boxplot(t(myGenesUta))
boxplot(t(myGenesEmc))
boxplot(t(matrix.combined2), col=rep(c(brewer.pal(n = 2, name = "Set1")[1:2]), 10), axes=F)
axis(side = 2, at = c(min(matrix.combined2, na.rm = T), -2, 0, 2, max(matrix.combined2, na.rm = T)), labels = c(-6.9, -2, 0, 2, 5), las=2)
axis(side = 1, at = seq(1.5, 20, 2), labels = sapply(ensToComplete[myGenes], function(x){unlist(strsplit(x, "_"))[2]}), las=2)
# abline(h=0)
legend("bottomright", legend = c("EMC", "UTA"), fill = brewer.pal(n = 2, name = "Set1")[1:2])

library(plyr)
fullMat2 <- rbind.fill.matrix(a, t(foldChangeIcgc), t(foldChangeEmc), t(foldChangeUta))
rownames(fullMat2) <- c(rownames(a), colnames(foldChangeIcgc), colnames(foldChangeEmc), colnames(foldChangeUta))
fullMat3 <- fullMat2[rownames(fullMat), ]
save(fullMat3, file="startData/patientsDegsFullMat.rda")

#################################################################################################
