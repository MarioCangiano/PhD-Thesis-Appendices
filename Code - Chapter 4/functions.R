
library(boot)
library(igraph)

extractBootInfo3 <- function(obj){
    a <- sapply(X = 1:ncol(obj$t), FUN = function(x){ifelse(obj$t0[x] > 0, sum(obj$t[, x] >= obj$t0[x], na.rm = T)/nrow(obj$t), 
                                                            sum(obj$t[, x] <= obj$t0[x], na.rm = T)/nrow(obj$t))})
    a <- c(obj$t0,a)
    names(a) <- c("TfScore", "pcScore", "MeanScores", "TfScore_p", "pcScore_p", "MeanScores_p")
    return(a)
}

correctForAliases <- function(x){
    library(limma)
    aliases <- alias2SymbolTable(alias = x, species = "Hs")
    aliases[is.na(aliases)] <- x[is.na(aliases)]
    message(sum(!x %in% aliases), " aliases corrected!")
    return(aliases)
}

getModulesList <- function(x, onlyGenes=T){
    newList <- list()
    for(i in 1:length(table(x@dynamicColors))){
        newList[[i]] <- toupper(x@peptides[which(x@dynamicColors == i)])
        if(onlyGenes){
            newList[[i]] <- sapply(newList[[i]], function(x){unlist(strsplit(x, "_"))[1]})
            a <- grep(";", newList[[i]])
            if(length(a) > 0){
                newList[[i]] <- newList[[i]][-1 * a]
            }
        }
    }
    names(newList) <- x@colorOrder[1:length(newList)]
    return(newList)
} 

cleanModules <- function(modules){
    i <- 1
    while(i <= length(modules)){
        if(length(modules[[i]]) > 1){
            a <- grep(";", modules[[i]])
            if(length(a) > 0){
                modules[[i]] <- modules[[i]][-1 * a]
            }
        }else{
            modules <- modules[-i]
            i <- i - 1
        }
        i <- i + 1
    }
    return(modules)
}

moduleAnalysis <- function(modules, exp, alg=NULL, mustLinks, cantLinks){
    expClean <- exp[!duplicated(sapply(rownames(exp), function(x){unlist(strsplit(x, split = "@"))[1]})), ]
    rownames(expClean) <- sapply(rownames(expClean), function(x){unlist(strsplit(x, split = "@"))[1]})
    a <- grep(";", rownames(expClean))
    expClean <- expClean[, ]
    submodules <- list()
    modules <- cleanModules(modules)
    pb <- txtProgressBar(min = 0, max = length(modules), style = 3)
    for(i in 1:length(modules)){
        subset <- expClean[intersect(rownames(expClean), modules[[i]]), ]
        subset <- subset[apply(subset, 1, function(x){sum(is.na(x))}) != ncol(subset), ]
        if(is.null(alg)){
            alg = "biclust"
        }
        switch(alg, 
               biclust={
                   require(biclust)
                   res <- biclust(x = as.matrix(subset), method = "BCCC")
                   submodules[[i]] <- list(biclust=res, clusters=extractClusters(res, rownames(subset)), alg=alg)
               },
               conclust={
                   require(conclust)
                   mustLink2 <- mustLinks[mustLinks[, 1] %in% rownames(subset) & mustLinks[, 2] %in% rownames(subset), ]
                   mustLink2 <- apply(mustLink2, 2, function(x){sapply(x, function(y){which(rownames(subset) == y)})})
                   cantLink2 <- cantLinks[cantLinks[, 1] %in% rownames(subset) & cantLinks[, 2] %in% rownames(subset), ]
                   cantLink2 <- apply(cantLink2, 2, function(x){sapply(x, function(y){which(rownames(subset) == y)})})
                   if(length(mustLink2) == 2){
                       mustLink2 <- t(as.data.frame(mustLink2))
                   }
                   if(length(cantLink2) == 2){
                       cantLink2 <- t(as.data.frame(cantLink2))
                   }
                   if(!is.null(dim(mustLink2)[1]) & !is.null(dim(cantLink2)[1])){
                       cantLink2 <- cantLink2[!apply(cantLink2, 1, paste, collapse="_") %in% apply(mustLink2, 1, paste, collapse="_") &
                                                  !apply(cantLink2, 1, paste, collapse="_") %in% apply(mustLink2, 1, function(z){paste(z[2], z[1], sep="_")}) , ]
                       res <- tryCatch(expr = ccls(data = subset, mustLink = mustLink2, cantLink = cantLink2), error=function(e){return(NA)})
                       if(!is.na(res) & sum(table(res)) != lenght(table(res))){
                           submodules[[i]] <- list(conclust=res, clusters=extractClusters2(res, rownames(subset)), alg=alg)
                       }else{
                           submodules[[i]] <- list(conclust=NA, clusters=unname(modules[[i]]), alg=alg)
                       }
                   }else{
                       submodules[[i]] <- list(conclust=NA, clusters=unname(modules[[i]]), alg=alg)
                   }
               }
        )
        setTxtProgressBar(pb, i)
    }
    close(pb)
    names(submodules) <- names(modules)
    return(submodules)
}

makeProperList <- function(x){
    newList <- c()
    k <- 0
    for(i in 1:length(x)){
        if(class(x[[i]]$clusters) == "list"){
            for(j in 1:length(x[[i]]$clusters)){
                k <- k + 1
                newList[[k]] <- x[[i]]$clusters[[j]]
                names(newList)[k] <- paste(names(x)[i], j, sep = "_")
            }
        }else{
            k <- k + 1
            newList[[k]] <- x[[i]]$clusters
            names(newList)[k] <- names(x)[i]
        }
    }
    return(newList)
}

executeEnrichment <- function(x, ncores){
    
    library(cogena)
    library(clusterProfiler)
    library(foreach)
    library(doSNOW)
    C <- ncores
    cl <- makeCluster(C)
    registerDoSNOW(cl)
    pb <- txtProgressBar(min = 1, max = length(x), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    x <- cleanModules(x)
    allOras <- foreach(i= 1:length(x), .errorhandling = "pass", .packages = "clusterProfiler", .options.snow = opts)%dopar%{
        overRepAnalysis <- enricher(gene = x[[i]], minGSSize = 2, TERM2GENE = read.gmt(gmtfile = "startData/CORUM.gmt"), pvalueCutoff = 1)
        return(overRepAnalysis)
    }
    stopCluster(cl)
    names(allOras) <- names(x)
    res <- extractOra(allOras)
    res$originalNumbers <- unlist(lapply(CORUMlist[res$ID], length))
    res <- res[order(res$ID, res$qvalue, decreasing = F), ]
    res <- res[!duplicated(res$ID), ]
    score <- sum(res$Count/res$originalNumbers, na.rm=T) * -log10(median(res$qvalue, na.rm = T))
    return(list(enrichment=res, score=score))
}

heatmapGG <- function(df, data.only = FALSE, cors=T) {  
    require(ggplot2) # ggplot2
    require(reshape2) # melt
    if(cors==F){
        test <- melt(cor(df,use = "na.or.complete"))
        listing <- list()
        for(x in 1:nrow(test)){
            listing[[x]] <- round(
                cor.test(
                    df[[test[x,1]]],
                    df[[test[x,2]]])$p.value,3
            )
        }
        test$Test <- unlist(listing)
    }else{
        test <- melt(df)
    }
    
    p1 <- ggplot(test, aes(Var1,Var2,fill = value, label = round(value,2))) + 
        geom_tile() + 
        geom_text() +
        scale_fill_distiller(palette = "Spectral", trans = "reverse") + 
        labs(
            x = "", 
            y = "",
            fill = "Correlation") + 
        theme_grey() +
        theme(
            axis.text.x = element_text(angle=90, size = 8),
            axis.text.y = element_text(size = 12),
            plot.background = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_blank()
        )
    if(data.only) {
        return(test)
    }else{
        return(p1)
    }
}

findOuterLinks <- function(myList){ # to be checked
    outerLinks <- c()
    innerLinks <- c()
    b <- unname(unlist(myList))
    pb <- txtProgressBar(min = 0, max = length(myList), style = 3)
    for(i in 1:length(myList)){
        a <- unlist(myList[[i]])
        innerLinks <- rbind(innerLinks, t(combn(x = a, m = 2)))
        d <- unique(b[!b %in% a])
        outerLinks <- rbind(outerLinks, expand.grid(a, d))
        setTxtProgressBar(pb, i)
    }
    close(pb)
    outerLinks <- outerLinks[!apply(outerLinks, 1, paste, collapse="_") %in% apply(innerLinks, 1, paste, collapse="_"), ]
    outer
    return(outerLinks)
}

findInnerLinks <- function(myList, ncores=8){
    
    library(foreach)
    require(igraph)
    library(doParallel)  
    cl <- parallel::makeCluster(ncores, type = "FORK")
    registerDoParallel(cl) 
    pb <- txtProgressBar(min = 1, max = length(myList), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    allLinks <- foreach(i= 1:length(myList), .errorhandling = "pass", .options.parallel = opts)%dopar%{
        myProts <- myList[[i]]
        if(length(myProts) > 1){
            a <- t(combn(x = myProts, m = 2))
            a <- rbind(a, cbind(a[, 2], a[, 1]))
        }else{
            a <- c()
        }
        return(a)
    }
    stopCluster(cl)
    return(do.call("rbind", allLinks))
}

listIntersection <- function(x, df=T){
    nms <- combn(names(x) , 2 , FUN = paste0 , collapse = "&" , simplify = FALSE )
    ll <- combn( x , 2 , simplify = FALSE )
    # if(ncores > 1){
    #   out <- parL
    # }else{
    out <- lapply( ll , function(x) ifelse(length( intersect( x[[1]] , x[[2]] ) ) > 1, 
                                           length( intersect( x[[1]] , x[[2]] ) )/ length(union(x[[1]], x[[2]])), 0))
    # }
    a <- unlist(setNames( out , nms ))
    if(df){
        b <- sapply(1:length(a), function(x){c(unlist(strsplit(names(a)[x], "&")), a[x])}, simplify = F)
        a <- as.data.frame(do.call("rbind", b), stringsAsFactors = F)
    }
    return(a)
} # modified on 201220 with jaccard index instead of intersection (for more than 1 in common)

cleanAdj <- function(x){
    x <- as.data.frame(x, stringsAsFactors=F)
    x$V3 <- as.numeric(x$V3)
    x <- x[!x$V1 == x$V2, ]
    return(x)
}

integrateRegulons3 <- function(netReg, modules, adjGenes, adjProt, multi=F){
    
    geneRange <- c(min(adjGenes$weigth), max(adjGenes$weigth))
    protRange <- c(min(adjProt$weight), max(adjProt$weight))
    
    rangeNormalisation <- function(x, range){ 
        x <- x[!is.na(x)]
        if(length(x) > 0){
            normData <- (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
            diff <- (range[2] - range[1])
            normData2 <- (normData * diff) + range[1] 
            normData2[is.na(normData2)] <- range[1]
        }else{
            normData2 <- c()
        }
        return(normData2)
    } 
    
    cleanAdj <- function(x){
        x <- as.data.frame(x, stringsAsFactors=F)
        x$V3 <- as.numeric(x$V3)
        x <- x[!x$V1 == x$V2, ]
        return(x)
    }
    
    make3Col <- function(x1, adj, threshold=0, split="@"){
        x2 <- sapply(colnames(adj)[which(abs(adj[x1, ]) > threshold)], function(x){unlist(strsplit(x, split))[1]})
        if(length(x2) > 0){
            x3 <- adj[x1, which(abs(adj[x1, ]) > threshold)]
            x1 <- unlist(strsplit(unname(rownames(adj)[x1]), split))[1]
            if(length(grep(";", x1)) == 0){
                a <- unname(cbind(toupper(rep(x1, length(x2))), toupper(x2), toupper(x3)))
                b <- grep(";", a[, 2])
                if(length(b) > 0){
                    a <- a[-1*b, ]
                }
                d <- grep("^$", a[, 2])
                if(length(d) > 0){
                    a <- a[-1*d, ]
                }
            }else{
                a <- c()
            }
        }else{
            a <- c()
        }
        return(a)
    }
    
    getGene <- function(x){
        sapply(x, function(y){unlist(strsplit(y, "_"))[2]})
    }
    
    getProt <- function(y, n){
        sapply(y, function(z){unlist(strsplit(z, "_"))[n]})
    }
    
    makeChain <- function(myList, reg){
        if(!is.na(myList)){
            if(nrow(myList) == 1){
                myList$reg <- rep(reg, nrow(myList))
                myList <- myList[, c(3,4,2,5)]
            }else{
                myList$from <- c(myList$to[nrow(myList)], myList$to[1:(nrow(myList)-1)])
                myList$reg <- rep(reg, nrow(myList))
                myList <- myList[, c(3,1,2,4)]
            }
            colnames(myList) <- c("from", "to", "weight", "reg")
        }else{
            myList <- c()
        }
        return(myList)
    }
    
    if(multi){
        expAdj <-  ends(graph = netReg[[i]], es = E(netReg[[i]]))
        expAdj <- as.data.frame(expAdj, stringsAsFactors = F)
        expAdj$weight <- E(netReg[[i]])$weight
        expAdj[, c(1:2)] <- apply(expAdj[, c(1:2)], 2, function(x){sapply(x, function(z){unlist(strsplit(z, "_"))[2]})})
        expAdj$weight[expAdj$weight > 0] <- expAdj$weight[expAdj$weight > 0]/geneRange[2]
        expAdj$weight[expAdj$weight < 0] <- expAdj$weight[expAdj$weight < 0]/-geneRange[1]
        expAdj$type <- E(netReg[[i]])$col
        posGenes <- unique(expAdj$V2)
        negGenes <- c()
        posIdx <- sapply(1:length(modules), function(y){match(allGenes, modules[[y]], nomatch=0)})
        negIdx <- c()
    }else{
        posGenes <- unlist(sapply(unlist(netReg[[1]]$pos), function(x){unlist(strsplit(x, "_"))[2]}))
        negGenes <- unlist(sapply(unlist(netReg[[1]]$neg), function(x){unlist(strsplit(x, "_"))[2]}))
        tf <- unlist(strsplit(names(netReg)[1], "_"))[2]
        posWeights <- filter(.data = adjGenes, IDs %in% paste(tf, posGenes, sep = "_"))
        negWeights <- filter(.data = adjGenes, IDs %in% paste(tf, negGenes, sep = "_"))
        expAdj <- as.data.frame(unname(cbind(rbind(posWeights, negWeights), c(rep("pos", nrow(posWeights)), rep("neg", nrow(negWeights))))), stringsAsFactors=F)
        expAdj$from <- sapply(expAdj[, 1], function(x){unlist(strsplit(x, "_"))[1]})
        expAdj$to <- sapply(expAdj[, 1], function(x){unlist(strsplit(x, "_"))[2]})
        rownames(expAdj) <- expAdj[, 1]
        expAdj <- expAdj[, c(4,5,2,3)]
        posIdx <- sapply(1:length(modules), function(y){match(getGene(posWeights$IDs), modules[[y]], nomatch=0)})
        negIdx <- sapply(1:length(modules), function(y){match(getGene(negWeights$IDs), modules[[y]], nomatch=0)})
    }
    
    pairPos <- which(posIdx != 0, arr.ind = T)
    if(length(pairPos) == 0){
        pairPos <- matrix(, nrow = 0, ncol = 0)
    }
    if(nrow(pairPos) > 0){
        complexesPos <- list()
        for(j in 1:nrow(pairPos)){
            message(j)
            links <- posGenes[pairPos[j, 1]]
            myModules <- names(modules)[pairPos[j, 2]] 
            myProts <- unname(unlist(modules[pairPos[j, 2]]))
            protComb <- apply(t(combn(myProts, 2)), 1, paste, collapse="_")
            posProtAdj <- filter(.data = adjProt, IDs %in% protComb)
            posProtAdj$from <- getProt(posProtAdj$IDs, 1) 
            posProtAdj$to <- getProt(posProtAdj$IDs, 2) 
            if(nrow(posProtAdj) > 1){
                complexesPos[[j]] <- posProtAdj %>% group_by(to) %>% summarise(mean(weight))
            }else{
                complexesPos[[j]] <- posProtAdj
            }
            names(complexesPos)[j] <- paste(myModules, links, sep="_")
        }
    }else{
        complexesPos <- list(NA)
    }
    
    pairNeg <- which(negIdx != 0, arr.ind = T)
    if(length(pairNeg) == 0){
        pairNeg <- matrix(, nrow = 0, ncol = 0)
    }
    if(nrow(pairNeg) > 0){
        complexesNeg <- list()
        for(j in 1:nrow(pairNeg)){
            message(j)
            links <- negGenes[pairNeg[j, 1]]
            myModules <- names(modules)[pairNeg[j, 2]] 
            myProts <- unname(unlist(modules[pairNeg[j, 2]]))
            protComb <- apply(t(combn(myProts, 2)), 1, paste, collapse="_")
            NegProtAdj <- filter(.data = adjProt, IDs %in% protComb)
            NegProtAdj$from <- getProt(NegProtAdj$IDs, 1) 
            NegProtAdj$to <- getProt(NegProtAdj$IDs, 2) 
            if(nrow(NegProtAdj) > 1){
                complexesNeg[[j]] <- NegProtAdj %>% group_by(to) %>% summarise(mean(weight))
            }else{
                complexesNeg[[j]] <- NegProtAdj
            }
            names(complexesNeg)[j] <- paste(myModules, links, sep="_")
        }
    }else{
        complexesNeg <- list(NA)
    }
    
    fullProtAdj <- rbind(do.call("rbind", lapply(complexesPos, makeChain, reg="pos")),
                         do.call("rbind", lapply(complexesNeg, makeChain, reg="neg")))
    colnames(fullProtAdj) <- colnames(expAdj) <- c("source", "target", "weight", "type")
    expAdj$origin <- rep("gene", nrow(expAdj))
    fullProtAdj$origin <- rep("protein", nrow(fullProtAdj)) # here we should add the name of the complex
    fullAdj <- rbind(expAdj, fullProtAdj)
    fullAdj <- fullAdj[!duplicated(fullAdj), ]
    g <- graph_from_data_frame(d = fullAdj,  directed = F)
    g <- set.graph.attribute(graph = g, name = "geneRange", value = geneRange)
    g <- set.graph.attribute(graph = g, name = "proteinRange", value = protRange)
    g <- set.graph.attribute(graph = g, name = "tf", value = tf)
    g <- set.graph.attribute(graph = g, name = "posComplexes", value = names(complexesPos))
    g <- set.graph.attribute(graph = g, name = "negComplexes", value = names(complexesNeg))
    # plot(g, layout=layout_with_sugiyama(g, layers=V(g)$origin)$layout, 
    #      vertex.size=5, edge.color=as.numeric(as.factor(E(g)$type)))
    return(g)
}

enrichSubGraph3 <- function(data, indices, x, protDegs, qvalueT=0.1){ # updated on 210319
    
    library(pracma)
    geneDegs <- data
    # colnames(geneDegs) <- c("log2FoldChange", "qvalue")
    geneDegs <- as.data.frame(geneDegs)
    protDegs <- as.data.frame(protDegs)
    if(!is.null(indices)){
        geneDegs$log2FoldChange <- geneDegs$log2FoldChange[indices]
        indices2 <- indices[indices %in% 1:nrow(protDegs)]
        protDegs$log2FoldChange[sort(indices2)] <- protDegs$log2FoldChange[indices2]
    }
    posGenesMax <- max(geneDegs[geneDegs$log2FoldChange > 0, "log2FoldChange"], na.rm = T) # is this really useful?
    negGenesMax <- -min(geneDegs[geneDegs$log2FoldChange < 0, "log2FoldChange"], na.rm = T)
    posProtMax <- max(protDegs[protDegs$log2FoldChange > 0, "log2FoldChange"], na.rm = T) # maybe we should filter only on the significant ones
    negProtMax <- -min(protDegs[protDegs$log2FoldChange < 0, "log2FoldChange"], na.rm = T)
    geneDegs <- geneDegs[geneDegs$qvalue <= qvalueT & !is.na(geneDegs$qvalue), ]
    protDegs <- protDegs[protDegs$qvalue <= qvalueT & !is.na(protDegs$qvalue), ]
    
    posGenes <- geneDegs[geneDegs$log2FoldChange > 0, ] # think about normalising before separation
    posGenes$weight2 <- posGenes$log2FoldChange/posGenesMax
    negGenes <- geneDegs[geneDegs$log2FoldChange < 0, ]
    negGenes$weight2 <- negGenes$log2FoldChange/negGenesMax
    degGenes2 <- rbind(posGenes, negGenes)
    degGenes2$weight3 <- 1 - degGenes2$qvalue
    
    # consider if we want to introduce transcription factor itself or not - default is not
    TfEdges <- as.data.frame(ends(graph = x, es = E(x)[E(x)$origin != "protein" ]), stringsAsFactors=F)
    TfEdges$weight1 <- E(x)$weight[E(x)[E(x)$origin != "protein"]]
    TfEdges <- TfEdges[!is.na(TfEdges$weight1), ]
    rownames(TfEdges) <- TfEdges$V2 # modified on 240121
    degGenesInTfEdges <- sum(rownames(degGenes2) %in% TfEdges$V2)
    if(degGenesInTfEdges > 0){
        myDegGenes <- rownames(degGenes2)[rownames(degGenes2) %in% TfEdges$V2]
        tfScoreRatio <- degGenesInTfEdges/nrow(TfEdges) # there could be a mistake here
        TfScore <- sum(nthroot(TfEdges[myDegGenes, "weight1"] * degGenes2[myDegGenes, "weight2"] * degGenes2[myDegGenes, "weight3"], 3), na.rm=T) * tfScoreRatio
    }else{
        TfScore <- 0
    }
    
    posProt <- protDegs[protDegs$log2FoldChange > 0, ] # think about using positive and negative complexes info
    posProt$weight2 <- posProt$log2FoldChange/posProtMax
    posProt$weight3 <- 1 - posProt$qvalue
    negProt <- protDegs[protDegs$log2FoldChange < 0, ]
    negProt$weight2 <- negProt$log2FoldChange/negProtMax
    negProt$weight3 <- 1 - negProt$qvalue
    
    protEdges <- as.data.frame(ends(graph = x, es = E(x)[E(x)$origin == "protein"]), stringsAsFactors=F) # this should be split by each complex
    protNodes <- unique(c(protEdges$V1, protEdges$V2))
    numberOfUniqueProts <- length(protNodes)
    protEdges$weight <- ifelse(E(x)$type[E(x)$origin == "protein"] == "pos", E(x)$weight[E(x)$origin == "protein"],
                               -E(x)$weight[E(x)$origin == "protein"])
    
    protEdgesPos <- protEdges[protEdges[, 1] %in% rownames(posProt) | protEdges[, 2] %in% rownames(posProt), ]
    numberOfUniquePosProts <- length(protNodes[protNodes %in% rownames(posProt)])
    if(numberOfUniquePosProts > 0){
        posProtNodes <- protNodes[protNodes %in% rownames(posProt)]
        protEdgesScorePos <- mean(nthroot(rep(mean(protEdgesPos$weight, na.rm=T), numberOfUniquePosProts) * posProt[posProtNodes, "weight2"] * 
                                              posProt[posProtNodes, "weight3"], 3), na.rm=T)
    }else{
        protEdgesScorePos <- 0
    }
    protEdgesNeg <- protEdges[protEdges[, 1] %in% rownames(negProt) | protEdges[, 2] %in% rownames(negProt), ]
    numberOfUniqueNegProts <- length(protNodes[protNodes %in% rownames(negProt)])
    if(numberOfUniqueNegProts > 0){
        negProtNodes <- protNodes[protNodes %in% rownames(negProt)]
        protEdgesScoreNeg <- mean(nthroot(rep(mean(protEdgesNeg$weight, na.rm=T), numberOfUniqueNegProts) * negProt[negProtNodes, "weight2"] * 
                                              negProt[negProtNodes, "weight3"], 3), na.rm=T)
    }else{
        protEdgesScoreNeg <- 0
    }
    posProtRatio <- numberOfUniquePosProts/numberOfUniqueProts
    negProtRatio <- numberOfUniqueNegProts/numberOfUniqueProts
    pcScore <- sum(c(protEdgesScorePos*posProtRatio, protEdgesScoreNeg*negProtRatio))
    pcScore <- ifelse(numberOfUniquePosProts == 0 & numberOfUniqueNegProts == 0, yes = 0, no = pcScore)
    
    # tfProtEdges <- TfEdges[TfEdges$V2 %in% protNodes, ]
    # ProtScores <- mean(tfProtEdges$weight1 * pcScore)
    # ProtScores <- ifelse(ProtScores > 0, sqrt(ProtScores), -sqrt(abs(ProtScores)))
    allScores <- mean(c(TfScore, pcScore), na.rm = T)
    finalScores <- c(TfScore, pcScore, allScores)
    return(finalScores)
} # last modified on 240221

calculateZscore3 <- function(y, nPerm=1000, geneDegs, protDegs, pvalueT=0.1, ncpu=20){
    geneDegs <- as.data.frame(geneDegs)
    protDegs <- as.data.frame(protDegs)
    library(boot)
    permScores <- boot(data = geneDegs, statistic = enrichSubGraph3, R = nPerm, x=y, protDegs=protDegs, parallel = "multicore", ncpus = ncpu)
    return(permScores)
}