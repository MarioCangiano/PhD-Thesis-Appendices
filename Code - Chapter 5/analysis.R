
####################################################################################### 

setwd("D:/mcangiano/Belfast/human_analysis_2/")
source("C:/Users/MarioC/Google Drive/Thesis material/Thesis/Chapter 3 - Paper/Scripts/functionsAndLibraries.R")

######################################################################################## info table

sampleSheet <- read.table("sampleCoding2.txt", sep="\t", quote = NULL, header = F, stringsAsFactors = F)
colnames(sampleSheet) <- c("GSid", "file")
myConv <- c("Mycophenolate mofetil", "Abiraterone", "ARN-509", "Mycophenolate mofetil and ARN-509", "Mycophenolate mofetil and Abiraterone",
            "benzyl alcohol and safflower oil", "benzyl alcohol and safflower oil", "Untreated")
names(myConv) <- c("ARM1", "ARM2", "ARM3", "ARM4", "ARM5", "Vehicle", "Untr")
sampleSheet$shortName <- myConv[sapply(sapply(sampleSheet$file, FUN = function(x){unlist(strsplit(x, "\\/"))[8]}), function(y){unlist(strsplit(y, "_"))[1]})]
sampleSheet$shortName <- gsub(pattern = " ", replacement = "_", sampleSheet$shortName)

# fl <- list.files(path = "disambiguate/", pattern = "*_summary.txt$", recursive = T, full.names = T)
# allSummaries <- c()
# for(i in 1:length(fl)){
#   a <- read.table(fl[i], sep = "\t", stringsAsFactors = F, quote = NULL, header = T)
#   allSummaries <- rbind(allSummaries, a)
# }
# myNames <- sapply(allSummaries$sample, function(x){unlist(strsplit(x, "-"))[4]})
# allSummaries$percHuman <- allSummaries$unique.species.A.pairs/apply(allSummaries[, 2:4], 1, sum)
# allSummaries$percMouse <- allSummaries$unique.species.B.pairs/apply(allSummaries[, 2:4], 1, sum)
# allSummaries$percAmb <- allSummaries$ambiguous.pairs/apply(allSummaries[, 2:4], 1, sum)
# save(allSummaries, file = "disambiguate/allSummaries.rda")
# par(mfrow=c(2, 1), mar=c(4,4,2,1))
# barplot(t(as.matrix(allSummaries[order(myNames), 2:4])), names.arg = sort(myNames), las=2,
#         col = c("green", "red", "orange"), main = "Disambiguity results on QUB dataset")
# barplot(t(as.matrix(allSummaries[order(myNames), 5:7])), names.arg = sort(myNames), las=2,
#         col = c("green", "red", "orange"), legend.text = c("human", "mouse", "ambiguous"),
#         args.legend = list(x=28, y=0.8, cex=0.7))
# par(mfrow=c(1,1), mar=c(5,4,4,4))
# 

############################################ PCA with normals, single and combination of drugs

load("startData/infoTable.rda")
sampleSheet$group2 <- "single"
sampleSheet$group2[grep("mofetil_and_", sampleSheet$group)] <- "combination"
sampleSheet$group2[c(grep("alcohol_and_", sampleSheet$group), grep("Untreated", sampleSheet$group))] <- "controls"
load("startData/normCounts.rda")

myPCa <- prcomp(normCounts)
plot(myPCa)
plot(myPCa$rotation, pch=16, lwd=4, col=as.numeric(as.factor(sampleSheet[rownames(myPCa$rotation), "group2"])), 
     main="PCA plot - normalized counts")
legend("bottomright", legend = levels(as.factor(sampleSheet[rownames(myPCa$rotation), "group2"])), fill= 1:3)

plot(myPCa$rotation, pch=16, lwd=4, col=as.numeric(as.factor(sampleSheet[rownames(myPCa$rotation), "group"])), 
     main="PCA plot - normalized counts")
# legend("bottomright", legend = levels(as.factor(sampleSheet[rownames(myPCa$rotation), "group"])), fill= 1:6)

############################################ common degs with Euler plots 

load("startData/degListsNew.rda")

ARM1 <- outListGenes$ARM1degs
ARM2 <- outListGenes$ARM2degs
# ARM1vsARM2 <- outListGenes$MycAbiDegs
ARM5vsARM2 <- outListGenes$ARM5vsARM2Degs
common1 <- intersect(rownames(ARM1), rownames(ARM2))
# common1 <- common1[(ARM2[common1, 1] < 0 & ARM1vsARM2[common1, 1] < 0) |
#                    (ARM2[common1, 1] > 0 & ARM1vsARM2[common1, 1] > 0) ]
common2 <- intersect(rownames(ARM2), rownames(ARM5vsARM2))
common3 <- intersect(rownames(ARM1), rownames(ARM5vsARM2))
common4 <- intersect(common2, common3)

library(eulerr)
fit1 <- euler(c("ARM2" = 24, "ARM1" = 0, "ARM5vsARM2" = 57,
                "ARM2&ARM1" = 0, "ARM2&ARM5vsARM2" = 6, "ARM1&ARM5vsARM2" = 0,
                "ARM2&ARM1&ARM5vsARM2" = 1))
plot(fit1, quantities = T, main="Common differentially expressed genes")

############################################################################### regulons enrichment - 150821 

load("../../Glasgow/RNAseq_orthografts/RTN/standardisedCounts_pearson/RTNsingleRegAracne.shadowFiltered.rda")
a <- netRegShadowFiltered$`ENSG00000119335_SET_9:131445702-131458679`
newDf <- cbind(a$pos, "positive")
colnames(newDf) <- c("Target", "Regulation")
library(xlsx)
write.xlsx(x = newDf, row.names = F, 
           file = "C:/Users/MarioC/Google Drive/Thesis material/Thesis/Chapter 5 - Treatment induced changes in gene regulatory network/Table 3.xlsx")

load("startData/degListsNew.rda")
set_genes <- unname(sapply(a$pos, function(x){unlist(strsplit(x, "_"))[1]}))
intersect(rownames(outListGenes$ARM5vsARM2Degs), set_genes) # 0

load("startData/normCounts.rda")
load("startData/infoTable.rda")
mySamples <- rownames(sampleSheet)[sampleSheet$group == "Mycophenolate_mofetil_and_Abiraterone" | sampleSheet$group == "Abiraterone"]
myRes <- topologyAnalysis2(h = normCounts, netReg = netRegShadowFiltered, degOutput = outListGenes$ARM5vsARM2Degs,
                        treatments = rownames(sampleSheet)[sampleSheet$group == "Mycophenolate_mofetil_and_Abiraterone"], controls = rownames(sampleSheet)[sampleSheet$group == "Abiraterone"], 
                        convert = T)

############################################################################### targets enrichment

myGenes <- unname(sapply(a$pos, function(x){(unlist(strsplit(x, "_"))[2])}))
myEnrichment <- enricher(gene = myGenes, TERM2GENE = read.gmt("../../Glasgow/RNAseq_orthografts/startData/c5.bp.v6.2.symbols.gmt"))
library(enrichplot)
cnetplot(myEnrichment)

load("startData/targetsEnrichment.rda")
a <- allEnrichmentsGenes$`ENSG00000119335_SET_9:131445702-131458679`@result
b <- a[a$qvalue < 0.05, c("ID", "GeneRatio", "pvalue", "qvalue")]
write.xlsx(x = b, row.names = F, 
           file = "C:/Users/MarioC/Google Drive/Thesis material/Thesis/Chapter 5 - Treatment induced changes in gene regulatory network/Table 4.xlsx")

###############################################################################
