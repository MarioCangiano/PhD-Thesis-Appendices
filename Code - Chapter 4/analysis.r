
source("functions.R")
setwd("D:/mcangiano/Glasgow/Proteomics_analysis/")

######################################################################## protein network

ortograftAnalysis <- read.table("Orthograft Proteomics Protein Groups processed/proteinGroups.txt", stringsAsFactors = F, sep = "\t", header = T)
rownames(ortograftAnalysis) <- paste(ortograftAnalysis$Gene.names, ortograftAnalysis$Protein.IDs, sep = "@")
ortograftMatrix <- log10(ortograftAnalysis[, grep("LFQ.intensity.L.", colnames(ortograftAnalysis))] + 1)
colnames(ortograftMatrix) <- infoTable[sapply(colnames(ortograftMatrix), function(x){unlist(strsplit(x, "\\."))[4]}), "GSid"]
ortograftRatioMatrix <- log10((ortograftL+1)/(ortograftH+1))

options(keep.source = TRUE, width = 70, stringsAsFactors=FALSE, digits=2)
library(ProCoNA)
colnames(ortograftMatrix) <- colnames(ortograftMatrixCorrected) <- colnames(ortograftRatio) <- gsub(pattern = "\\.", replacement = "_", colnames(ortograftMatrix))
groups <- sapply(colnames(ortograftMatrix), function(x){unlist(strsplit(x, "_"))[4]})
crpcStatus <- c(rep(T, 3), rep(F, 3), rep(F, 3), rep(T, 3), rep(F, 3), rep(T, 3))
allowWGCNAThreads()
ortograftMatrix[ortograftMatrix == 0] <- 0.1
beta <- pickSoftThreshold(t(ortograftMatrix), networkType="signed", RsquaredCut=0.8)
proteinNetwork <- buildProconaNetwork(networkName="ortograftMatrix", pepdat=t(ortograftMatrix), networkType="signed", pearson = F, toPermTestPermutes=1000)

ortograftModules <- cleanModules(getModulesList(proteinNetwork, onlyGenes = T))
ortograftModulesSplitted <- makeProperList(moduleAnalysis(modules = ortograftModules, exp = ortograftMatrix, alg = "biclust"))
resOrtograftModulesSplitted <- executeEnrichment(x = ortograftModulesSplitted, ncores = 20) 

########################################################################################## gene/protein correlations

load("startData/allStartingMatrices.rda")
load("../RNAseq_orthografts/startData/infoTable.rda")
load("F:/mcangiano/Glasgow/RNAseq_orthografts/startData/normCounts.rda")
newNames <- unname(sapply(rownames(normCounts), function(x){toupper(unlist(strsplit(x, "_")))[2]}))
normCounts <- normCounts[!duplicated(newNames), ]
rownames(normCounts) <- unname(sapply(rownames(normCounts), function(x){toupper(unlist(strsplit(x, "_")))[2]}))

library(doSNOW)
c <- snow::makeCluster(8, type = "SOCK")
registerDoSNOW(cl = c)
pb <- txtProgressBar(max = length(myList), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
resList <- foreach(i=2:length(myList), .multicombine = F, .options.snow = opts) %dopar% {
  a <- myList[[i]]
  numberOfGenes <- c()
  allCors <- matrix(0, 18, 18)
  newNames <- unname(sapply(rownames(a), function(x){toupper(unlist(strsplit(x, "@")))[1]}))
  common <- intersect(newNames, rownames(normCounts))
  for(j in 1:18){
    for(k in 1:18){
      allCors[j, k] <- cor(a[common, j], normCounts[common, k], use = "complete.obs")
      numberOfGenes <- c(numberOfGenes, sum(!is.na(a[common, j])))
    }
  }
  rownames(allCors) <- colnames(a)
  colnames(allCors) <- colnames(normCounts)
  attr(allCors, which = "numberOfGenes") <- numberOfGenes
  return(allCors)
}
stopCluster(c)
stopImplicitCluster()
names(resList) <- names(myList)[-1]
save(resList, file = "startData/corRes_expression.rda")

load("startData/corRes_expression.rda")
myMat <- resList$ortograftL
colnames(myMat) <- rownames(myMat) <- infoTable[rownames(myMat), "Sample_ID"]
heatmapGG(df = myMat)

########################################################################################## DEPs

library(ROTS)
load("startData/allStartingMatrices.rda")
load("../RNAseq_orthografts/startData/infoTable.rda")
infoTable[colnames(myList$ortograftRatioMatrix), "CRPC"] # it's making 1 for CRPC and 0 for PC

rots.out <- ROTS(data=ortograftRatioMatrix, groups = ifelse(1:ncol(ortograftRatioMatrix) %in% c(1:3, 7:9, 13:15), yes = 0, no = 1),
                 B = 2000, K = 1500, progress = T) 
depCRPC_LH <- rots.out
save(depCRPC_LH, file = "Deps/CRPC_LH.rda")

load("../RNAseq_orthografts/DEG/DEG_matrices_rawCounts_disambiguated.rda")
load("Deps/CRPC_LH.rda")
load("../RNAseq_orthografts/startData/convertedGenesFromRawCountsFromMygene.rda")
geneDegs <- CRPC[, c("padj", "log2FoldChange")]
newNames <- sapply(ensToComplete[rownames(CRPC)], function(x){unlist(strsplit(x, "_"))[2]})
geneDegs <- geneDegs[!duplicated(newNames) & !is.na(newNames), ]
rownames(geneDegs) <- newNames[!duplicated(newNames) & !is.na(newNames)]
geneDegs <- geneDegs[, c("log2FoldChange", "padj")]
colnames(geneDegs) <- c("log2FoldChange", "qvalue")
protDegs <- cbind(-depCRPC_LH$logfc, depCRPC_LH$FDR) # the grouping was incorrect
rownames(protDegs) <- sapply(rownames(protDegs), function(x){toupper(unlist(strsplit(x, "@")))[1]})
a <- grep(";", rownames(protDegs))
protDegs <- protDegs[-a, ]
colnames(protDegs) <- c("log2FoldChange", "qvalue")
protDegs <- as.data.frame(protDegs, stringsAsFactors=F)
geneDegs <- as.data.frame(geneDegs[!is.na(geneDegs$qvalue), ])
save(geneDegs, protDegs, file = "startData/cleanedCRPCdegdep2.rda")

########################################################################################### genes and prots differential expression comparison

load("startData/cleanedCRPCdegdep2.rda") # updated on 220319
common <- intersect(rownames(geneDegs), rownames(protDegs))
common2 <- intersect(rownames(geneDegs)[geneDegs$qvalue < 0.05], rownames(protDegs)[protDegs$qvalue < 0.05])
library(eulerr)
fit <- euler(c("GENES" = nrow(geneDegs), "PROTEINS" = nrow(protDegs),
               "GENES&PROTEINS" = length(common)))
plot(fit, quantities = T)
fit <- euler(c("DEGS" = sum(geneDegs$qvalue < 0.05), "DEPS" = sum(protDegs$qvalue < 0.05),
               "DEGS&DEPS" = length(common2)))
plot(fit, quantities = T)

commonTable <- cbind(geneDegs[common2, ], protDegs[common2, ]) # written on 020621
colnames(commonTable) <- paste(c(rep("transcript_", 2), rep("peptides_", 2)), colnames(commonTable), sep = "")
Sys.setenv(JAVA_HOME = "C:/Program Files (x86)/Java/jdk-16/")
library(xlsx)
write.xlsx(x = commonTable, file = "C:/Users/MarioC/Google Drive/Thesis material/Thesis/Chapter 4 - Proteomics and transcriptomics/Additional_data/common_degs_deps.xlsx", sheetName = "CRPC_vs_PC_common_significant_features", row.names = T, col.names = T)

myCols <- rep(NA, length(common))
myCols[which(geneDegs[common, "qvalue"] < 0.05)] <- "red"
myCols[which(protDegs[common, "qvalue"] < 0.05)] <- "green"
names(myCols) <- common
myCols[common2] <- "cyan" 
plot(x = geneDegs[common, 'log2FoldChange'], col=myCols,lwd= 2, pch=16, xlab = "Log2 fold change - gene expression",
     ylab = "log2 fold change - protein expression",
     y = protDegs[common, 'log2FoldChange'], main="9 CRPC vs 9 PC orthografts")
legend("topright", legend = c("DEG - 256", "DEP - 152", "DE at both layers - 16"), 
       fill=c("red", "green", "cyan"), cex=0.8)
abline(h = 0)
abline(v = 0)

setwd(dir = "D:/mcangiano/Glasgow/Proteomics_analysis/")
load("startData/cleanedCRPCdegdep2.rda")
significantGenes <- rownames(geneDegs)[geneDegs$qvalue < 0.05]
significantGenes <- unname(sapply(significantGenes, function(x){unlist(strsplit(x, "\\."))[1]}))
significantProt <- rownames(protDegs)[protDegs$qvalue < 0.05]
significantProt <- unname(sapply(significantProt, function(x){unlist(strsplit(x, "\\."))[1]}))

library(clusterProfiler)
myGmt <- read.gmt(gmtfile = "../../generalData/biologicalInfo/h.all.v7.2.symbols.gmt")
enrichGenes <- enricher(gene = significantGenes, TERM2GENE = myGmt) # HALLMARK_ANDROGEN_RESPONSE and HALLMARK_UV_RESPONSE_DN
enrichProt <- enricher(gene = significantProt, TERM2GENE = myGmt) # HALLMARK_OXIDATIVE_PHOSPHORYLATION and HALLMARK_MYC_TARGETS_V1
write.xlsx(x = enrichGenes@result, file = "C:/Users/MarioC/Google Drive/Thesis material/Thesis/Chapter 4 - Proteomics and transcriptomics/Additional_data/Overrepresentation_analysis.xlsx",
           sheetName = "genes_enrichment", row.names = T, col.names = T)
write.xlsx(x = enrichProt@result, file = "C:/Users/MarioC/Google Drive/Thesis material/Thesis/Chapter 4 - Proteomics and transcriptomics/Additional_data/Overrepresentation_analysis.xlsx",
           sheetName = "proteins_enrichment", row.names = T, col.names = T, append = T)

#################################  bootstrap analysis on UGLA CRPCdex - all regulons - updated on 220319

load("startData/cleanedCRPCdegdep2.rda")
load("integrativeRegulons/singleLayerUndirected/integrativeGraphs.rda")

allRes <- list()
pb <- txtProgressBar(min = 1, max = length(integrativeRegulons), style = 3)
for(i in 1:length(integrativeRegulons)){
  
  allRes[[i]] <- tryCatch(expr = {calculateZscore3(y = integrativeRegulons[[i]], nPerm = 1000, ncpu = 7,  geneDegs = geneDegs,
                                                   protDegs = protDegs, pvalueT = 0.05)},
                          error=function(e){e})
  setTxtProgressBar(pb, i)
}
close(pb)
names(allRes) <- names(integrativeRegulons)
save(allRes, file = "enrichments/CRPCdegdep2.rda")

load("enrichments/CRPCdegdep2.rda")
a <- as.data.frame(do.call("rbind", lapply(allRes, FUN = extractBootInfo3)))
plot(a$TfScore_p, a$pcScore_p)
abline(h = 0.05)
abline(v = 0.05)
a[which(a$TfScore_p < 0.05 & a$pcScore_p < 0.05), ]
a <- cbind(a, apply(X = a[, 4:6], MARGIN = 2, p.adjust, method="fdr")) # the qvalue is exactly 0 because there was no permutation providing an higher score
colnames(a)[7:9] <- gsub(colnames(a)[7:9], pattern = "_p", replacement = "_q") 
a[which(a$TfScore_q < 0.05 & a$pcScore_q < 0.05), ]
uglaRes <- a
save(uglaRes, file = "enrichments/allRegsCRPC_newScore/enrichmentResults/CRPCdegdep2_res.rda")

############################################################# MID1 plot

getFeatures <- function(x, y){
  myEdges <- as_ids(E(x)[E(x)$origin == y])
  myFeatures <- sapply(myEdges, function(x){unlist(strsplit(x, "\\|"))[2]})
  return(unique(myFeatures))
}

library(igraph)
load("integrativeRegulons/singleLayerUndirected/integrativeGraphs.rda")
a <- integrativeRegulons$`ENSG00000101871_MID1_X:10413349-10851773`
b <- graph.attributes(a)
onlyGenes <- induced_subgraph(a, vids=c("MID1", getFeatures(a, "gene")))
onlyProts <- induced_subgraph(a, vids=c(getFeatures(a, "protein")))

load("startData/cleanedCRPCdegdep2.rda")
significantGenes <- rownames(geneDegs)[geneDegs$qvalue < 0.05]
significantGenes <- unname(sapply(significantGenes, function(x){unlist(strsplit(x, "\\."))[1]}))
significantProt <- rownames(protDegs)[protDegs$qvalue < 0.05]
significantProt <- unname(sapply(significantProt, function(x){unlist(strsplit(x, "\\."))[1]}))

sigGenes <- geneDegs[geneDegs$qvalue < 0.05, "log2FoldChange"]
names(sigGenes) <- significantGenes
genesCol <- ifelse(sigGenes[V(onlyGenes)$name] > 0, "red", "blue")
genesCol[is.na(genesCol)] <- "gray"

alternativeLayout <- layout_with_sugiyama(onlyGenes, layers=V(onlyGenes)$origin)$layout
plot(onlyGenes, vertex.label.cex=c(2,rep(1, length(V(onlyGenes))-1)), vertex.label.family="Calibri", layout = layout.circle(onlyGenes), vertex.color=genesCol,
     vertex.size=round(abs(E(onlyGenes)$weight)*10), edge.color=c("blue", "red")[as.numeric(as.factor(E(onlyGenes)$type))], vertex.label.dist=1,
     vertex.label.color="black", vertex.label.font = 2, edge.size =3)

load("startData/CORUMlistOnly.rda")
mySets <- CORUMlist[gsub(pattern = "_[A-Z]+[0-9]+", replacement = "", b$posComplexes)]
numberOfSets <- c()
setsList <- list()
for(i in V(onlyProts)$name){
  relatedSets <- unlist(sapply(1:length(mySets), function(x){if(i %in% mySets[[x]]){return(names(mySets)[x])}}))
  z <- length(relatedSets)
  setsList[[i]] <- relatedSets
  numberOfSets <- c(numberOfSets, z)
}
names(numberOfSets) <- names(setsList) <- V(onlyProts)$name
protsCol <- c("blue", "orange", "cyan")[as.numeric(as.factor(numberOfSets))]
# protsCol <- c("yellow", "orange", "purple", "pink")[as.numeric(as.factor(numberOfSets))]

plot(onlyProts, vertex.label.family="Calibri", layout = layout.fruchterman.reingold(onlyProts), vertex.color=protsCol,
     vertex.size=round(abs(E(onlyProts)$weight)*10), edge.color=c("blue", "red")[as.numeric(as.factor(E(onlyProts)$type))], vertex.label.dist=1,
     vertex.label.color="black", vertex.label.font = 2, edge.size =3, vertex.shape= ifelse(names(V(onlyProts)) %in% significantProt, yes = "square", no = "circle"))
legend("bottomleft", legend = c(1:3), fill = c("blue", "orange", "cyan"), title = "Number of complexes")
legend("bottomleft", legend = c("Mitochondrial rib.", "Multiple", "Nop56p pre rRNA", "Cytoplasmic rib."), fill = c("yellow", "orange", "purple", "pink"))

plot(onlyProts, vertex.label.family="Calibri", layout = layout.fruchterman.reingold(onlyProts), vertex.color=protsCol,
     vertex.size=round(abs(E(onlyProts)$weight)*10), edge.color=c("blue", "red")[as.numeric(as.factor(E(onlyProts)$type))], vertex.label.dist=1,
     vertex.label.color="black", vertex.label.font = 2, edge.size =3)
legend("bottomleft", legend = c(1:3), fill = c("blue", "red", "green"))
legend("bottomleft", legend = c("Mitochondrial rib.", "Multiple", "Nop56p pre rRNA", "Cytoplasmic rib."), fill = c("yellow", "orange", "purple", "pink"))

common <- intersect(names(sigGenes), names(sigProts))
sigProts <- protDegs[protDegs$qvalue < 0.05, "log2FoldChange"]
names(sigProts) <- significantProts
links <- unique(unlist(sapply(b$posComplexes, function(x){unlist(strsplit(x, "_"))[length(unlist(strsplit(x, "_")))]})))
myGenes <- sigGenes[V(a)$name]
myProts <- sigProts[V(a)$name]
myGenes[is.na(myGenes)] <- 0
myProts[is.na(myProts)] <- 0
names(myGenes) <- names(myProts) <- V(a)$name
vertexCols <- unlist(sapply(V(a)$name, function(x){if(myGenes[x] > 0 & myProts[x] > 0)
{"purple"}else{if(myGenes[x] > 0){"firebrick"}else{if(myProts[x] > 0){"forestgreen"}else{"grey"}}}}))
alternativeLayout <- layout_with_sugiyama(a, layers=V(a)$origin)$layout # reingold.tilford, star are nice, 
plot(a, vertex.label.family="Calibri", layout = layout.reingold.tilford(a), vertex.color=vertexCols, vertex.label.cex=c(2,rep(0.5, length(V(a))-1)),
     vertex.size=5, edge.color=c("blue", "red")[as.numeric(as.factor(E(a)$type))], vertex.label.dist=0,
     vertex.label.color="black", vertex.label.font = 2, edge.size =3, vertex.label=NA)
legend("bottomright", legend = c("Gene upregulation", "Protein upregulation", "Gene and protein upregulation"), fill = c("firebrick", "forestgreen", "purple"))

myEnds <- ends(graph = onlyProts, es = E(onlyProts))

###########################################################################
