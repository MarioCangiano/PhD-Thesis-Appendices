
#################################################################### Tampere

load("/mnt/data/mcangiano/Transpot/7060_Tampere/Rna-seq/startData/infoTable.rda")
load("/mnt/data/mcangiano/Transpot/7060_Tampere/Rna-seq/startData/rawCounts.rda")
load("/mnt/data/mcangiano/Transpot/7060_Tampere/Rna-seq/DEGs/deseq2_res.rda")

colnames(rawCounts) <- rownames(sampleCoding)
# sum(colnames(normCounts) == colnames(normCountsDeseq)) # 58
stData <- Surv(sampleCoding[colnames(normCountsDeseq), "Progression.free.days.after.diagnosis"], 
               ifelse(sampleCoding[colnames(normCountsDeseq), "Progression"] == "Yes", 1, 0))
rownames(stData) <- colnames(normCountsDeseq)
normalsCounts <- normCountsDeseq[, rownames(sampleCoding)[sampleCoding$class == "BPH"]]
normCounts <- normCountsDeseq[, rownames(sampleCoding)[sampleCoding$class != "BPH"]]
sampleCoding$numberOfCounts <- colSums(rawCounts)[rownames(sampleCoding)]
analysisInfo <- read.table("/mnt/data/mcangiano/Transpot/7060_Tampere/Rna-seq/mappingQC_multiqc_report_data/multiqc_rseqc_read_distribution.txt", 
                           stringsAsFactors = F, header = T, sep = "\t")
rownames(analysisInfo) <- sapply(analysisInfo$Sample, function(x){unlist(strsplit(x, "\\."))[1]})
analysisInfo <- analysisInfo[, grep("pct", colnames(analysisInfo))]
sampleCoding$readDistributionMean <- apply(analysisInfo, 1, mean)[rownames(sampleCoding)]
sampleCoding$readDistributionSd <- apply(analysisInfo, 1, sd)[rownames(sampleCoding)]
fastqInfo <- read.table("/mnt/data/mcangiano/Transpot/7060_Tampere/Rna-seq/fastQC_multiqc_report_data/multiqc_fastqc.txt", 
                        stringsAsFactors = F, header = T, sep = "\t")
rownames(fastqInfo) <- fastqInfo$Sample
sampleCoding$GC_content <- fastqInfo[rownames(sampleCoding), "X.GC"] 
sampleCoding$targetClass <- ifelse((sampleCoding$Progression == "No" & sampleCoding$class =="PC") | sampleCoding$class == "BPH", yes = 0, no = 1)
sampleCoding$progressedAfterADTin7years <- sampleCoding$targetClass
sampleCoding$progressedAfterADTin7years[sampleCoding$targetClass == 0 & sampleCoding$Progression.free.days.after.diagnosis > 2520] <- 2
table(sampleCoding$progressedAfterADTin7years, sampleCoding$class)[, 3] 
# 0  1  2 
# 7 16  5 

attr(x = normCounts, which = "stData") <- stData[colnames(normCounts), ]
attr(x = normCounts, which = "metadata") <- sampleCoding[colnames(normCounts), ]
attr(x = normCounts, which = "covariates") <- c(5,7,8,10:15,18,20:ncol(sampleCoding))
attr(x = normCounts, which = "targetProg") <- c(17,24)
attr(x = normCounts, which = "targetClass") <- 25
attr(x = normCounts, which = "name") <- "tampere_deseq2"
attr(x = normCounts, which = "normals") <- normalsCounts
attr(x = normCounts, which = "code") <- "ensembl"
attr(x = normCounts, which = "symbols") <- sapply(ensToComplete[rownames(normCounts)], function(x){unlist(strsplit(x, "_"))[2]})
attr(x = normCounts, which = "completeName") <- ensToComplete[rownames(normCounts)]
save(normCounts, file = "datasets/tampereDeseq2.rda")

######################################################### Rotterdam patients 

load("/mnt/data/mcangiano/Transpot/7060_Rotterdam/7046-04/startData/infoTable.rda")
load("/mnt/data/mcangiano/Transpot/7060_Rotterdam/7046-04/startData/rawCounts.rda")
load("/mnt/data/mcangiano/Transpot/7060_Rotterdam/7046-04/DEGs/deseq2_res.rda")

normalsCounts <- normCountsDeseq[, rownames(metaData2)[metaData2$caseControl == "normal"]]
normCounts <- normCountsDeseq[, rownames(metaData2)[metaData2$caseControl != "normal"]]
metaData2$numberOfCounts <- colSums(rawCounts)[rownames(metaData2)] 
analysisInfo <- read.table("/mnt/data/mcangiano/Transpot/7060_Rotterdam/7046-04/mappingQC_multiqc_report_data/multiqc_rseqc_read_distribution.txt", 
                           header = T, sep="\t", quote=NULL, stringsAsFactors = F)
rownames(analysisInfo) <- sapply(analysisInfo$Sample, function(x){unlist(strsplit(x, "\\."))[1]})
analysisInfo <- analysisInfo[, grep("pct", colnames(analysisInfo))]
metaData2$readDistributionMean <- apply(analysisInfo, 1, mean)[rownames(metaData2)]
metaData2$readDistributionSd <- apply(analysisInfo, 1, sd)[rownames(metaData2)]
# metaData2["7046-004-015", "readDistributionMean"] <- mean(as.numeric(read.table("/mnt/data/mcangiano/Transpot/7060_Rotterdam/7046-04/rseqc/7046-004-015.readDistribution.log", 
#                                                                             header = F, quote=NULL, stringsAsFactors = F, fill = T)[6:15, 4]))
# metaData2["7046-004-015", "readDistributionSd"] <- sd(as.numeric(read.table("/mnt/data/mcangiano/Transpot/7060_Rotterdam/7046-04/rseqc/7046-004-015.readDistribution.log", 
#                 header = F, quote=NULL, stringsAsFactors = F, fill = T)[6:15, 4]))
fastqInfo <- read.table("/mnt/data/mcangiano/Transpot/7060_Rotterdam/7046-04/fastQC_multiqc_report_data/multiqc_fastqc.txt", 
                        header = T, sep="\t", quote=NULL, stringsAsFactors = F)
rownames(fastqInfo) <- fastqInfo$Sample
fastqInfo <- fastqInfo[grep("R1", rownames(fastqInfo)), ]
rownames(fastqInfo) <- sapply(rownames(fastqInfo), function(x){unlist(strsplit(x, "_"))[1]})
metaData2$GC_content <- fastqInfo[rownames(metaData2), "X.GC"]
metaData2$metastasisAfterEndo <- ifelse(metaData2$EndoTreatYN == 1 & metaData2$MorNPosYN == 1 & metaData2$MorNPosDate >= metaData2$EndoTreatDate
                                        & metaData2$caseControl == "tumor", 1, 0)
a <- metaData2[metaData2$EndoTreatYN == 1 & !is.na(metaData2$EndoTreatYN), ]
boxplot(LastFuMnthSincePCaDiag ~ MorNPosYN, a, xlab=c("MorNPosYN of patients treated with endocrine therapy"), ylab="Last follow up time", 
        main="Rotterdam_metastasis_time")
stripchart(LastFuMnthSincePCaDiag ~ MorNPosYN, vertical = TRUE, data = a, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
metaData2$metastasisAfterEndo <- ifelse(metaData2$EndoTreatYN == 1 & metaData2$MorNPosYN == 1 & metaData2$MorNPosDate >= metaData2$EndoTreatDate
                                        & metaData2$caseControl == "tumor", 1, 0)
metaData2$metastasisAfterEndoWithin13years <- metaData2$metastasisAfterEndo
metaData2$metastasisAfterEndoWithin13years[metaData2$MorNPosYN == 0 & metaData2$EndoTreatYN == 1 & metaData2$LastFuMnthSincePCaDiag > 156] <- 2
table(metaData2$metastasisAfterEndoWithin13years)

stData <- Surv(metaData2[colnames(normCounts), "PSAProgAfterRP01MnthSinceRP"], 
               metaData2[colnames(normCounts), "PSAProgAfterRP01YN"])
rownames(stData) <- colnames(normCounts)

attr(x = normCounts, which = "stData") <- stData
attr(x = normCounts, which = "metadata") <- metaData2[colnames(normCounts), c(2,5,16:47, 49:56)]
attr(x = normCounts, which = "covariates") <- c(3:12, 38:41)
attr(x = normCounts, which = "targetProg") <- c(13,15)
attr(x = normCounts, which = "targetClass") <- c(42)
attr(x = normCounts, which = "name") <- "Rotterdam_deseq2"
attr(x = normCounts, which = "normals") <- normalsCounts
attr(x = normCounts, which = "code") <- "ensembl"
attr(x = normCounts, which = "symbols") <- sapply(ensToComplete[rownames(normCounts)], function(x){unlist(strsplit(x, "_"))[2]})
attr(x = normCounts, which = "completeName") <- ensToComplete[rownames(normCounts)]
save(normCounts, file = "datasets/RotterdamDeseq2.rda")

############################################################################# ICGC

rawCounts <- read.table("startData/mrna135.gc19.GENE.RAW_COUNT.featureCounts.tsv.gz", header = T, sep = "\t", stringsAsFactors = F)
colnames(rawCounts) <- gsub("_T0[0-9]$", "", colnames(rawCounts))
rownames(rawCounts) <- rawCounts$sample
rawCounts <- rawCounts[, -1]
save(rawCounts, file = "startData/rawCounts.rda")
library(xlsx)
infoTable <- read.xlsx("EOPC clinical data.xlsx", sheetIndex = 2, startRow = 2, header = T, stringsAsFactors=F)
rownames(infoTable) <- infoTable$PID
save(infoTable, file = "startData/infoTable.rda")

load("startData/rawCounts.rda")
load("startData/infoTable.rda")
missing <- setdiff(colnames(rawCounts), rownames(infoTable))
names(missing) <- sapply(missing, function(x){unlist(strsplit(x, "\\."))[1]})
infoTable2 <- infoTable
for(i in 1:length(missing)){
    if(names(missing)[i] %in% rownames(infoTable2)){
        infoTable2 <- rbind(infoTable2, rep(infoTable2[names(missing)[i], ], 1))
    }else{
        infoTable2 <- rbind(infoTable2, rep(NA, ncol(infoTable2)))
    }
    rownames(infoTable2)[nrow(infoTable2)] <- missing[i]
}
rawCounts2 <- rawCounts[, colnames(rawCounts) %in% rownames(infoTable2)]
infoTable2 <- infoTable2[colnames(rawCounts2), ]
infoTable2$Type[is.na(infoTable2$Type)] <- "NA"
pc <- prcomp(x = rawCounts2, scale. = T)
plot(pc$rotation[, 1], pc$rotation[, 2], col=as.numeric(as.factor(infoTable2[rownames(pc$rotation), "Type"])), pch=15, lwd=3)
plot(pc$rotation[, 1], pc$rotation[, 2], col=as.numeric(as.factor(infoTable2[rownames(pc$rotation), "Clonal_status"])), pch=15, lwd=3)

a <- sapply(colnames(rawCounts2), function(x){unlist(strsplit(x, "\\."))[1]})
library(Matrix.utils)
rawCounts3 <- t(as.matrix(aggregate.Matrix(x = t(rawCounts2), groupings = factor(a), FUN=sum)))
infoTable3 <- infoTable2[-grep("\\.", rownames(infoTable2)), ]
infoTable3 <- infoTable3[colnames(rawCounts3), ]
infoTable3$Type[is.na(infoTable3$Type)] <- c("normal", "NA")
pc <- prcomp(x = rawCounts3, scale. = T)
plot(pc$rotation[, 1], pc$rotation[, 2], col=as.numeric(as.factor(infoTable3[rownames(pc$rotation), "Type"])), pch=15, lwd=3)
plot(pc$rotation[, 1], pc$rotation[, 2], col=as.numeric(as.factor(infoTable3[rownames(pc$rotation), "Clonal_status"])), pch=15, lwd=3)
plot(pc$rotation[, 1], pc$rotation[, 2], col=as.numeric(as.factor(infoTable3[rownames(pc$rotation), "ETS_status"])), pch=15, lwd=3)
infoTable3$Gleason[infoTable3$Gleason == "4+4" | infoTable3$Gleason == "5+5"] <- NA
infoTable3[is.na(infoTable3)] <- "NA"
library(DESeq2)
dds0 <- DESeqDataSetFromMatrix(countData = rawCounts3,
                               colData = infoTable3, 
                               design= ~ Type + Gleason + Stage + Clonal_status)
dds0 <- DESeq(dds0, parallel = T)
deseqNorm <- counts(dds0, normalize=T) # identical norm matrix by changing formula with CRPC
plotPCA(vst(dds0), intgroup = "Type")
rownames(deseqNorm) <- sapply(rownames(deseqNorm), function(x){unlist(strsplit(x, "\\."))[1]})
save(deseqNorm, file = "startData/deseqNormCounts_aggregated.rda")

#############################################################################
