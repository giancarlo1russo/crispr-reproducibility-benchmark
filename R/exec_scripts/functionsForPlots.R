#############################################################
scoreVsRankCorrPlot <- function(gg, testLibs=FALSE) {
gg <- lapply(names(gg), function(p) {
  pp <- gg[[p]]
  pp[["Genes"]] <- sapply(rownames(pp), function(y) unlist(strsplit(y, "\\."))[1])
  pp$data <- p
  return(pp)
})
ggAll <- do.call("rbind",gg)
ggAll <- ggAll[ggAll$Genes != "All",]
metrics <- unique(ggAll$data)
genes <- unique(ggAll$Genes)

comparisonLists <- list(Metric1 = lapply(metrics, function(z) paste(c(FALSE,TRUE),z, sep = ".")),
                        GenesLib = lapply(genes, function(z) paste(c(FALSE,TRUE),z, sep = ".")),
                        Metric2=lapply(c(FALSE,TRUE), function(z) paste(z,metrics, sep = ".")),
                        Genes = lapply(data.frame(combn(genes,2)), function(x) c(x)))
if (!testLibs) {
tempPlot1 <- ggplot(ggAll, aes(data, Correlation)) + geom_boxplot() + 
  stat_compare_means(label = "p.signif", vjust = 0.25, comparisons = list(c("rank","zscores"))) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(vars(Genes),vars(attributeType)) 
tempPlot2 <- ggplot(ggAll, aes(Genes, Correlation)) + geom_boxplot() + 
  stat_compare_means(comparisons = comparisonLists$Genes,label = "p.signif", vjust = 0.25) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(vars(attributeType),vars(data)) 
}else {

  ggAll$sameLib <-as.factor(sapply(ggAll$libUsed, addLibraryInfo))
  ggAll <- ggAll[!is.na(ggAll$Correlation),]
  tempPlot1 <-  ggplot(ggAll, aes(x=interaction(sameLib,data), Correlation)) + 
    geom_boxplot(aes(fill=sameLib)) + 
    stat_compare_means(comparisons = c(comparisonLists$Metric1, comparisonLists$Metric2), label = "p.signif", vjust = 0.25) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(vars(Genes),vars(attributeType), scales = "free") 
  tempPlot2 <- ggplot(ggAll, aes(x=Genes, Correlation)) +
    geom_boxplot() + 
    stat_compare_means(label = "p.signif", vjust = 0.25, comparisons = comparisonLists$Genes)  +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    facet_grid(vars(attributeType),vars(data), scales="free") 
}
return(list(tempPlot1,tempPlot2))
}

#############################################################
###### REAL VS RANDOM RANKS
realVsRandomRankCorrPlot <- function(fullListOfScoresRankOverlapsRealFile, 
                                     fullListOfScoresRankOverlapsRandomFile,testLibs=TRUE,
                                     sameLibReal=FALSE, sameLibRandom=FALSE) {
  finalDF_List <- lapply(c(zscores="zscores",rank="rank"), function(metric) {
    realRanks <- fullListOfScoresRankOverlapsRealFile[[metric]]
    randomRank <- fullListOfScoresRankOverlapsRandomFile[[metric]]
    if (testLibs) {
    if (sameLibReal) {
      realRanks <- realRanks[selectBasedOnlibs(realRanks),]
    } else {
      realRanks <- realRanks[selectBasedOnlibs(realRanks, sameLib = FALSE),]
  }
    if (sameLibRandom) {
      randomRank <- randomRank[selectBasedOnlibs(randomRank),]
    } else {
      randomRank <- randomRank[selectBasedOnlibs(randomRank, sameLib = FALSE),]
    }
    }
    realRanks[["Genes"]] <- sapply(rownames(realRanks), function(y) unlist(strsplit(y, "\\."))[1])
       realRankList <- split(realRanks, realRanks$attributeType)
    randomRankList <- split(randomRank, randomRank$AttributeType)
    listOfPAttributes <- lapply(names(realRankList), function(p) {
      realRankTemp <- realRankList[[p]]
      randomVals <- randomRankList[[p]]
      randomValsPermList <- split(randomVals, randomVals$PermutationNr)
      realRankTempGenesList <- split(realRankTemp, realRankTemp$Genes)
      listOfPAttributesTempList <- lapply(names(realRankTempGenesList), function(pp) {
      realVals <- realRankTempGenesList[[pp]][["Correlation"]]
      listOfP <- sapply(randomValsPermList, function(l) {
        wilcox.test(l[["Correlation"]], realVals, alternative="less")$p.value
      })
      return(data.frame(Attribute=p, Pvalues=listOfP, Genes=pp))
      })
      listOfPAttributesTemp <- do.call("rbind",listOfPAttributesTempList)
      return(listOfPAttributesTemp)
    })
    listOfPAttributes_DF <- do.call("rbind",listOfPAttributes)
    listOfPAttributes_DF$pAdj <- p.adjust(listOfPAttributes_DF$Pvalues)
    listOfPAttributes_DF$negLogP <- -log10(listOfPAttributes_DF$pAdj)
    listOfPAttributes_DF$Metric <- metric
    return(listOfPAttributes_DF)
  })
  finalDF <- do.call("rbind",finalDF_List)
   tempPlot <-  ggplot(finalDF, aes(Attribute, negLogP)) + geom_jitter() + 
     geom_abline(slope=0, intercept = -log10(0.05), color="red") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
     facet_grid(rows=vars(Metric), cols=vars(Genes))
   return(tempPlot)
}

#### REAL SAME ATTR. VS. DIFF. ATTR. 
#############################################################
sameVsDiffAttrCorrPlot <- function(fullListOfScoresRankOverlapsRealFile) {
  finalDF_List <- lapply(c(zscores="zscores",rank="rank"), function(pp) {
    realRanks <- fullListOfScoresRankOverlapsRealFile[[pp]]
    sameRanks <- realRanks[selectBasedOnlibs(realRanks, sameLib = TRUE),]
    sameRanks[["Genes"]] <- sapply(rownames(sameRanks), function(y) unlist(strsplit(y, "\\."))[1])
    sameRanksList <- split(sameRanks, sameRanks[["attributeType"]])
    diffRanks <- realRanks[selectBasedOnlibs(realRanks, sameLib = FALSE),]
    diffRanks[["Genes"]] <- sapply(rownames(diffRanks), function(y) unlist(strsplit(y, "\\."))[1])
    diffRanksList <- split(diffRanks, diffRanks[["attributeType"]])
      listOfPAttributes <- lapply(names(sameRanksList), function(p) {
        sameRanksListTemp <- sameRanksList[[p]]
        sameRanksGeneTempList <- split(sameRanksListTemp, sameRanksListTemp$Genes)
        diffRanksListTemp <- diffRanksList[[p]]
        diffRanksGeneTempList <- split(diffRanksListTemp, diffRanksListTemp$Genes)
    geneTypeBasedList  <- lapply(names(sameRanksGeneTempList), function(geneType){
      sameVals <- sameRanksGeneTempList[[geneType]][["Correlation"]]
      diffVals <- diffRanksGeneTempList[[geneType]][["Correlation"]]
      pVal <-   wilcox.test(diffVals,sameVals,  alternative="less")$p.value
      return(data.frame(Attribute=p, Pvalues=pVal, Genes=geneType))
    })
    listOfPAttributesByGeneTemp <- do.call("rbind",geneTypeBasedList)
    return(listOfPAttributesByGeneTemp)
      })
    listOfPAttributes_DF <- do.call("rbind",listOfPAttributes)
    listOfPAttributes_DF$pAdj <- p.adjust(listOfPAttributes_DF$Pvalues)
    listOfPAttributes_DF$negLogP <- -log10(listOfPAttributes_DF$pAdj)
    listOfPAttributes_DF$Metric <- pp
    return(listOfPAttributes_DF)
  })
  finalDF <- do.call("rbind",finalDF_List)
  tempPlot <-  ggplot(finalDF, aes(Attribute, negLogP)) + geom_jitter() + 
    geom_abline(slope=0, intercept = -log10(0.05), color="red") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(vars(Metric),vars(Genes),scales="free")
  return(tempPlot)
}



#### CHECK  CORRELATION LEVELS BETWEEN DIFFERENT ATTRIBUTE TYPES 
#############################################################
interAttrCorrPlot <- function(fullListOfScoresRankOverlapsRealFile) {
  pp <- lapply(fullListOfScoresRankOverlapsRealFile, function(x) {
    x[["Genes"]] <- sapply(rownames(x), function(y) unlist(strsplit(y, "\\."))[1]) 
    x <- x[x[["Genes"]]=="All",]
    return(x)
  })
mergedScoreRanksCorrelations  <- rbind(cbind(pp$zscores, ScoreType="Score"),
                                       cbind(pp$rank, ScoreType="Rank"))
attributeList <- unique(mergedScoreRanksCorrelations$attributeType)
ScoreRanksCorrelationPlotByAttrAll <- ggplot(mergedScoreRanksCorrelations, aes(attributeType,Correlation)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", comparisons = as.list(data.frame(combn(attributeList,2)))) +
  facet_wrap(vars(ScoreType))
return(ScoreRanksCorrelationPlotByAttrAll)
}


#### INTRA-ATTRIBUTE TO CHECK WHICH FEATURES ARE MORE REPRODUCIBLE
#############################################################
intraAttrCorrPlot <- function(fullListOfScoresRankOverlapsRealFile) {
  pp <- lapply(fullListOfScoresRankOverlapsRealFile, function(x) {
    x[["Genes"]] <- sapply(rownames(x), function(y) unlist(strsplit(y, "\\."))[1]) 
    x <- x[x[["Genes"]]=="All",]
    return(x)
  })
  mergedScoreRanksCorrelations  <- rbind(cbind(pp$zscores, ScoreType="Score"),
                                         cbind(pp$rank, ScoreType="Rank"))
  mergedScoreRanksCorrelationsSplit <- split(mergedScoreRanksCorrelations, 
                                             mergedScoreRanksCorrelations[["attributeType"]])
  listOfPlots <- lapply(mergedScoreRanksCorrelationsSplit, function(y) {
    pp2 <- DF_factorOrder(y, factorToOrder = "attributeName",orderingMetric = "Correlation",
                          aggregate=TRUE, aggrStat="mean")
    tempPlot <- ggplot(pp2, aes(attributeName,Correlation)) + geom_boxplot() +
      facet_wrap(vars(ScoreType), nrow = 2) + coord_flip() 
    return(tempPlot)
  })
  return(listOfPlots)
}

#### TOP PERCENT HYPERTEST
#############################################################
top_percent_comparison <- function(top_perc_overlap_DF,isReal=TRUE, testLibs=FALSE){
top_perc_overlap_DF$Sign <- ifelse(top_perc_overlap_DF$adjP < 0.05, "YES","NO")
top_perc_overlap_DF$Overlap <- ifelse(isReal, "Real","Random")
if (!testLibs) {
  df_to_return <- top_perc_overlap_DF[,c("negLog10P","Percentile","Sign","Overlap")]
} else {
  top_perc_overlap_DF$sameLib <-sapply(top_perc_overlap_DF$Libraries, function(x) {
    y <- length(unique(unlist(strsplit(x,","))))
    z<- ifelse(y==1, TRUE, FALSE)
    return(z)})
  df_to_return <- top_perc_overlap_DF[,c("negLog10P","Percentile","Sign","sameLib","Overlap")]
}
df_to_return$Attribute <- unique(top_perc_overlap_DF$AttributeType)
return(df_to_return)
}
