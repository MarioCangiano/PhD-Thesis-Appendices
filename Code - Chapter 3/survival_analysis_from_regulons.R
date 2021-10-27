
#####################################################################################################

setwd("D:/mcangiano/survivalAnalysis/")
source("C:/Users/MarioC/Google Drive/Thesis material/Thesis/Chapter 3 - Paper/Scripts/functionsAndLibraries.R")
load("D:/mcangiano/Glasgow/RNAseq_orthografts/RTN/standardisedCounts_pearson/RTNsingleRegAracne.shadowFiltered.rda")

##################################################################################################### Tampere

load("regulonsEnrichment/datasets_deseq2/tampereDeseq2.rda") 
at <- attributes(normCounts)
mySamples <- rownames(at$metadata)
mySamples <- mySamples[-length(mySamples)] # "7068-01-004-058" is an outlier for gene body coverage
target <- ifelse(at$metadata[mySamples, "Progression"] == "Yes", "red", "green")
names(target) <- mySamples
target <- target[!is.na(target)]
meta <- at$metadata[names(target), ]
meta <- meta[!is.na(meta$Progression.free.days.after.diagnosis), ]
meta$Progression <- ifelse(meta$Progression == "Yes", 1, 0)
continuousCov <- c("Age.at.diagnosis", "PSA.at.diagnosis", "Progression.free.days.after.diagnosis", "Gleason.score")
discreteCov <- c("T.class", "M.class")
meta[, continuousCov] <- apply(meta[, continuousCov], MARGIN = 2, as.numeric)
meta$T.class2 <- as.numeric(factor(meta$T.class, levels = c("T1c", "T2", "T2a", "T3"), ordered = T))
meta$M.class <- factor(meta$M.class)
coxph(Surv(Progression.free.days.after.diagnosis, Progression) ~ Gleason.score, data = meta) 
coxph(Surv(Progression.free.days.after.diagnosis, Progression) ~ meta$T.class2, data = meta) 

load("ggeaEnrichment4_netRegShadowFiltered/ggeaTampereMat2.rda")
myMat <- tampereRes$pval
myMat2 <- myMat[names(target), ]
tampereMat <- myMat2
tampereMeta <- meta[rownames(myMat2), c("Progression.free.days.after.diagnosis", "Progression", "Gleason.score", "T.class")]
tampereCounts <- normCounts[, rownames(myMat2)]

myReg2 <- myMat2[, "ENSG00000070495_JMJD6_17:74708918-74722866"]
JMJD6 <- ifelse(myReg2 == 1, "Active", "Inactive")
myMeta0 <- cbind(meta[names(JMJD6), ], JMJD6)
myMeta0$Progression.free.months.after.diagnosis <- myMeta0$Progression.free.days.after.diagnosis/30
coxph(Surv(Progression.free.months.after.diagnosis, Progression) ~ myReg2, data = myMeta0) 
coxph(Surv(Progression.free.months.after.diagnosis, Progression) ~ myReg2 + Gleason.score + T.class2, data = myMeta0) 
ggsurvplot(legend.labs = c("JMJD6 ACTIVE", "JMJD6 INACTIVE"), xlab="Time (months)", pval.method = T, data = myMeta0, pval = T, risk.table="nrisk_cumevents",
           survfit(Surv(Progression.free.months.after.diagnosis, Progression) ~ JMJD6, data = myMeta0), title="UTA")

p <- ggsurvplot(legend.labs = c("JMJD6 ACTIVE", "JMJD6 INACTIVE"), xlab="Time (months)", pval.method = T, data = myMeta0, pval = F, risk.table="nrisk_cumevents",
           survfit(Surv(Progression.free.months.after.diagnosis, Progression) ~ JMJD6, data = myMeta0), title="UTA")
p$plot <- p$plot +
    annotate("text", x = 5, y = 0.25, label = "Log-rank\np= 6e-04", cex=4, col="black", vjust=0, hjust = 0, fontface=3)
p


load("signatures_ensembl/Georgescu.rda")
myVect <- apply(normCounts[signature, ], 2, mean)
myClass4 <- as.numeric(myVect > quantile(x = myVect, probs = 0.67)) + 1
names(myClass4) <- colnames(normCounts)
myClass4 <- myClass4[rownames(myMat2)]
myMeta0 <- cbind(meta[names(myClass4), ], myClass4)
coxph(Surv(Progression.free.days.after.diagnosis, Progression) ~ myClass4 , data = myMeta0) 
coxph(Surv(Progression.free.days.after.diagnosis, Progression) ~ myClass4 + Gleason.score + T.class2, data = myMeta0)

load("../Tampere/Rna-seq/startData/bromoGSVA.rda") 
myVect <- tampereBromo["Bromo.10", rownames(myMat2)]
myMeta0 <- cbind(meta[names(myVect), ], myVect)
coxph(Surv(Progression.free.days.after.diagnosis, Progression) ~ myVect , data = myMeta0)
coxph(Surv(Progression.free.days.after.diagnosis, Progression) ~ myVect + Gleason.score + T.class2, data = myMeta0) 

load("signatures_ensembl/Yang.rda")
myVect <- colSums(normCounts[names(signature), ] * signature)
myClass4 <- as.numeric(myVect > quantile(x = myVect, probs = 0.5)) 
names(myClass4) <- colnames(normCounts)
myClass4 <- myClass4[rownames(myMat2)]
myMeta0 <- cbind(meta[names(myClass4), ], myClass4)
coxph(Surv(Progression.free.days.after.diagnosis, Progression) ~ myClass4, data = myMeta0) 

################################################## EMC

load("regulonsEnrichment/datasets_deseq2/rotterdam.rda")
at <- attributes(normCounts)
meta <- at$metadata
meta$PSAProgAfterRP01MnthSinceRP[meta$PSAProgAfterRP01YN == 0] <- meta$LastFuMnthSincePCaDiag[meta$PSAProgAfterRP01YN == 0]
target <- meta$PSAProgAfterRP01YN
names(target) <- rownames(meta)
meta[meta==9] <- NA
meta$TMPRSS2_ERG[meta$TMPRSS2_ERG == "NO – reclassified from YES to NO after second PCR (TMPRSS2-PADI4 positive)"] <- "NO"
meta$fullGleason <- as.factor(meta$pGleason1 + meta$pGleason2)
meta$pTStad_num <- as.numeric(ordered(meta$pTStad, levels=c("x", "2x", "2a", "2b", "2c", "3x", "3a", "3b", "3c", "4x", "4a", "4b")))
target <- target[target != 9]
target <- ifelse(target == 0, "green", "red")
target <- target[names(target) %in% rownames(meta)[which(meta$MorNPosYN == 0)]] 
target <- target[-1*c(1:2)] # reclassified samples
target <- target[rownames(meta)[meta$EndoTreatYN == 0]]
target <- target[!is.na(target)]
coxph(Surv(PSAProgAfterRP01MnthSinceRP, PSAProgAfterRP01YN) ~ fullGleason, data = meta[names(target), ]) 
coxph(Surv(PSAProgAfterRP01MnthSinceRP, PSAProgAfterRP01YN) ~ pTStad_num , data = meta[names(target), ]) 

load("ggeaEnrichment4_netRegShadowFiltered/ggeaEmcMat2.rda")
myMat <- ggeaEmcRes$pval
myMat2 <- myMat[names(target), ]
meta2 <- meta[meta$EndoTreatYN == 0, ]
common <- intersect(rownames(myMat2), rownames(meta2))
myMat3 <- myMat2[common, ]
emcMat <- myMat3
emcMeta <- meta[rownames(myMat3), c("PSAProgAfterRP01MnthSinceRP", "PSAProgAfterRP01YN", "fullGleason", "pTStad")]
emcCounts <- normCounts[, rownames(myMat3)]

myReg2 <- myMat3[, "ENSG00000070495_JMJD6_17:74708918-74722866"]
JMJD6 <- ifelse(myReg2 == 1, "ACTIVE", "INACTIVE")
myMeta0 <- cbind(meta[names(JMJD6), ], JMJD6)
coxph(Surv(PSAProgAfterRP01MnthSinceRP, PSAProgAfterRP01YN) ~ JMJD6 , data = myMeta0) 
coxph(Surv(PSAProgAfterRP01MnthSinceRP, PSAProgAfterRP01YN) ~ myReg2 + pTStad_num, data = myMeta0) 
ggsurvplot(xlab="Time (months)", pval.method = T, data = myMeta0, pval = T, risk.table="nrisk_cumevents",
           survfit(Surv(PSAProgAfterRP01MnthSinceRP, PSAProgAfterRP01YN) ~ JMJD6, data = myMeta0), title="EMC")

load("signatures_ensembl/Georgescu.rda")
myVect <- apply(normCounts[signature, ], 2, mean)
myClass4 <- as.numeric(myVect > quantile(x = myVect, probs = 0.67)) + 1
names(myClass4) <- colnames(normCounts)
myClass4 <- myClass4[rownames(myMat3)]
myMeta0 <- cbind(meta[names(myClass4), ], myClass4)
coxph(Surv(PSAProgAfterRP01MnthSinceRP, PSAProgAfterRP01YN) ~ myClass4 , data = myMeta0) 

load("../Rotterdam/7046-04/startData/bromoGSVA.rda")
myVect <- rotterdamBromo["Bromo.10", rownames(myMat3)]
myMeta0 <- cbind(meta[names(myVect), ], myVect)
coxph(Surv(PSAProgAfterRP01MnthSinceRP, PSAProgAfterRP01YN) ~ myVect , data = myMeta0) 

load("signatures_ensembl/Yang.rda")
myVect <- colSums(normCounts[names(signature), ] * signature)
myClass4 <- as.numeric(myVect > quantile(x = myVect, probs = 0.5)) 
names(myClass4) <- colnames(normCounts)
myClass4 <- myClass4[rownames(myMat3)]
myMeta0 <- cbind(meta[names(myClass4), ], myClass4)
coxph(Surv(PSAProgAfterRP01MnthSinceRP, PSAProgAfterRP01YN) ~ myClass4, data = myMeta0) 

###################################################################### ICGC

load("D://mcangiano/ICGC/startData/deseqNormCounts.rda")
load("D://mcangiano/ICGC//startData/infoTable.rda")
meta <- infoTable[colnames(deseqNorm), ]
target <- as.numeric(meta$BCR)
names(target) <- rownames(meta)
target <- target[!is.na(target)] 
target <- ifelse(target == 0, "green", "red")
meta$Gleason[meta$Gleason == "NA"] <- NA
meta <- meta[names(target), ]
meta[, c(8, 10, 11, 13)] <- apply(meta[, c(8, 10, 11, 13)], 2, as.numeric)
meta$Stage_num <- as.numeric(factor(meta$Stage, levels = c("pT2a", "pT2c", "pT3a", "pT3b", "pT4"), ordered = T))
meta$Gleason_num <- as.numeric(factor(meta$Gleason, levels = c("3+3", "3+4", "4+3", "4+4", "4+5", "5+4", "5+5"), ordered = T))
conVect <- c(7,9,7,6,9,10,8)
names(conVect) <- unique(meta$Gleason)[!is.na(unique(meta$Gleason))]
meta$Gleason_num2 <- conVect[meta$Gleason]
target <- target[meta$Primary_.Treatment == "RadP"]
meta <- meta[names(target), ]
coxph(Surv(Time.from.surgery.to.BCR.or.lastFU, BCR) ~ Stage_num, data = meta) 
coxph(Surv(Time.from.surgery.to.BCR.or.lastFU, BCR) ~ Gleason_num2, data = meta) 

load("ggeaEnrichment4_netRegShadowFiltered/ggeaIcgcMat2.rda") 
library(Matrix.utils)
myMat <- ggeaIcgcRes$pval
a <- sapply(rownames(myMat), function(x){unlist(strsplit(x, "\\."))[1]})
myMat2 <- as.matrix(aggregate.Matrix(x = myMat, groupings = factor(a), fun="sum"))[names(target), ]
myMat2[myMat2 > 1] <- 1
meta2 <- meta[meta$Primary_.Treatment == "RadP", ]
common <- intersect(rownames(myMat2), rownames(meta2))
myMat3 <- myMat2[common, ]
icgcMat <- myMat3
icgcMeta <- meta[rownames(myMat3), c("Time.from.surgery.to.BCR.or.lastFU", "BCR","Gleason", "Stage")]
icgcCounts <- deseqNorm[, rownames(myMat3)]

myReg2 <- myMat3[, "ENSG00000070495_JMJD6_17:74708918-74722866"]
JMJD6 <- ifelse(myReg2 == 1, "ACTIVE", "INACTIVE")
myMeta0 <- cbind(meta[names(JMJD6), ], JMJD6)
coxph(Surv(Time.from.surgery.to.BCR.or.lastFU, BCR) ~ JMJD6 , data = myMeta0) 
coxph(Surv(Time.from.surgery.to.BCR.or.lastFU, BCR) ~ myReg2 + Stage_num + Gleason_num2, data = myMeta0) 
ggsurvplot(xlab="Time (months)", pval.method = T, data = myMeta0, pval = T, surv.median.line = "v", risk.table="nrisk_cumevents",
           survfit(Surv(Time.from.surgery.to.BCR.or.lastFU, BCR) ~ JMJD6, data = myMeta0), title="ICGC")

icgcDeseq <- deseqNorm
load("signatures_ensembl/Georgescu.rda")
myVect <- apply(icgcDeseq[signature, ], 2, mean)
myClass4 <- as.numeric(myVect > quantile(x = myVect, probs = 0.67)) + 1
names(myClass4) <- colnames(icgcDeseq)
myClass4 <- myClass4[rownames(myMat3)]
myMeta0 <- cbind(meta[names(myClass4), ], myClass4)
coxph(Surv(Time.from.surgery.to.BCR.or.lastFU, BCR) ~ myClass4 , data = myMeta0) 
coxph(Surv(Time.from.surgery.to.BCR.or.lastFU, BCR) ~ myClass4 + Stage_num + Gleason_num2, data = myMeta0)

load("../ICGC/startData/rawCounts.rda")
a <- sapply(rownames(rawCounts), function(x){unlist(strsplit(x, "\\."))[1]})
rawCounts2 <- as.matrix(aggregate.Matrix(x = rawCounts, groupings = factor(a), fun="mean"))
b <- sapply(colnames(rawCounts2), function(x){unlist(strsplit(x, "\\."))[1]})
rawCounts3 <- t(as.matrix(aggregate.Matrix(x = t(rawCounts2), groupings = factor(b), fun="mean")))
load("../Glasgow/RNAseq_orthografts/startData/convertedGenesFromRawCountsFromMygene.rda")
rownames(rawCounts3) <- sapply(ensToComplete[rownames(rawCounts3)], function(x){unlist(strsplit(x, "_"))[2]})
icgcBromo <- makeBromoScore(fpkms = rawCounts3, sig = "../ICGC/startData/signatures_symbols.gmt")
myVect <- icgcBromo["Bromo.10", rownames(myMat3)]
myMeta0 <- cbind(meta[names(myVect), ], myVect)
coxph(Surv(Time.from.surgery.to.BCR.or.lastFU, BCR) ~ myVect , data = myMeta0) 

load("signatures_ensembl/Yang.rda")
myVect <- colSums(icgcDeseq[names(signature), ] * signature)
myClass4 <- as.numeric(myVect > quantile(x = myVect, probs = 0.5)) 
names(myClass4) <- colnames(icgcDeseq)
myClass4 <- myClass4[rownames(myMat3)]
myMeta0 <- cbind(meta[names(myClass4), ], myClass4)
coxph(Surv(Time.from.surgery.to.BCR.or.lastFU, BCR) ~ myClass4, data = myMeta0)
coxph(Surv(Time.from.surgery.to.BCR.or.lastFU, BCR) ~ myClass4 + Stage_num + Gleason_num2, data = myMeta0) 

#############################################################################################