
############################################################################################################# 

setwd("D://mcangiano/survivalAnalysis/")

load("D://mcangiano/Glasgow/RNAseq_orthografts/RTN/standardisedCounts_pearson/RTNsingleRegAracne.shadowFiltered.rda")

source("functionsAndLibraries.R")

########################################################################################## Tampere

load("regulonsEnrichment/datasets_deseq2/tampere.rda")
at <- attributes(normCounts)
fullMat <- cbind(normCounts, at$normals)

ssDegs <- makeSingleSampleDegs2(normCounts = fullMat, normals = colnames(at$normals), all = F)

foldChangeUta <- do.call(what = "cbind", args = lapply(ssDegs, function(x){x[, 1]}))
qvalueUta <- do.call(what = "cbind", args = lapply(ssDegs, function(x){x[, 2]}))

ggeaTampere1b <- topologyAnalysis2(h = fullMat, netReg = netRegShadowFiltered, outFile = "ggeaEnrichment4_netRegShadowFiltered/tampereRes1b.rda", treatments = colnames(normCounts), 
                                  controls = colnames(at$normals), singleSample = T, ncores = 20, ssDegs = ssDegs)

############################################################################################### Rotterdam

load("regulonsEnrichment/datasets_deseq2/rotterdam.rda")
at <- attributes(normCounts)
fullMat <- cbind(normCounts, at$normals)

ssDegs <- makeSingleSampleDegs2(normCounts = fullMat, normals = colnames(at$normals), all = F)

foldChangeEmc <- do.call(what = "cbind", args = lapply(ssDegs, function(x){x[, 1]}))
qvalueEmc <- do.call(what = "cbind", args = lapply(ssDegs, function(x){x[, 2]}))

ggeaRotterdam2 <- topologyAnalysis2(h = fullMat, netReg = netRegShadowFiltered, outFile = "ggeaEnrichment4_netRegShadowFiltered/ggeaRotterdamRun2.rda", treatments = colnames(normCounts), 
                                  controls = colnames(at$normals), singleSample = T, ncores = 20, ssDegs = ssDegs)

############################################################################################### ICGC

load("regulonsEnrichment/datasets_deseq2/icgc.rda")
at <- attributes(normCounts)
fullMat <- cbind(normCounts, at$normals)

ssDegs <- makeSingleSampleDegs2(normCounts = fullMat, normals = colnames(at$normals), all = F)

foldChangeIcgc <- do.call(what = "cbind", args = lapply(ssDegs, function(x){x[, 1]}))
qvalueIcgc <- do.call(what = "cbind", args = lapply(ssDegs, function(x){x[, 2]}))

ggeaIcgc2 <- topologyAnalysis2(h = fullMat, netReg = netRegShadowFiltered, outFile = "ggeaEnrichment4_netRegShadowFiltered/ggeaIcgcRun2.rda", treatments = colnames(normCounts), 
                             controls = colnames(at$normals), singleSample = T, ncores = 30, ssDegs = ssDegs)

########################################################################################## save results

load("ggeaEnrichment4_netRegShadowFiltered/tampereRes1b.rda")
tampereRes <- list(norm=extractNormScore(nbea.res), pval=pvalueMat(nbea.res, adjusted = T))
save(tampereRes, file = "ggeaEnrichment4_netRegShadowFiltered/ggeaTampereMat2.rda")

load("ggeaEnrichment4_netRegShadowFiltered/ggeaRotterdamRun2.rda")
ggeaEmcRes <- list(norm=extractNormScore(nbea.res), pval=pvalueMat(nbea.res, adjusted = T))
save(ggeaEmcRes, file = "ggeaEnrichment4_netRegShadowFiltered/ggeaEmcMat2.rda")

load("ggeaEnrichment4_netRegShadowFiltered/ggeaIcgcRun2.rda")
ggeaIcgcRes <- list(norm=extractNormScore(nbea.res), pval=pvalueMat(nbea.res, adjusted = T))
save(ggeaIcgcRes, file = "ggeaEnrichment4_netRegShadowFiltered/ggeaIcgcMat2.rda")

##########################################################################################

