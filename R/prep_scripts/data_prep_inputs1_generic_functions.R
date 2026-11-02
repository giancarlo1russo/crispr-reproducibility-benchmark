
##########################################################
generate_overlaps_prep <- function(allInfoFilePath,attribute_list,attribute_levels) { 
  allInfoFile <- readRDS(allInfoFilePath)

  
  if (any(names(allInfoFile)=="SCORE")) {
    names(allInfoFile)[names(allInfoFile)=="SCORE"] <- "BCSCORE" 
  }

  allInfoFileSplitList <- lapply(attribute_list, function(attr) {
    print(attr)
    split_data <- split(allInfoFile, allInfoFile[[attr]])
    split_data <- split_data[names(split_data)!=""]
    attribute_level <-
      attribute_levels[
        attribute_levels$Attribute == attr,
      ]$Level
    if (attribute_level < 4) {
      print(attribute_level)
    higher_level_name <-
      attribute_levels[
        attribute_levels$Level == (attribute_level + 1),
      ]$Attribute
    list_of_hln <- unique(allInfoFile[[higher_level_name]])
    list_of_un_attr <-  sapply(list_of_hln, function(x) {
      list_of_un_attr_temp <- length(unique(allInfoFile[allInfoFile[[higher_level_name]]==x,][[attr]]))
      if (list_of_un_attr_temp > 1) {
        return(x)
      } else {
        return(NULL)
      }
    })
    list_of_un_attr[sapply(list_of_un_attr,is.null)] <- NULL
    list_of_un_attr <- unname(unlist(list_of_un_attr))
    return(list(split_data=split_data,hln_to_check=list_of_un_attr))
  } else {
    return(list(split_data=split_data,hln_to_check=NULL))
  }
  })
  
  attributeSummary <- lapply(attribute_list, function(attribute) {
    gg <- allInfoFileSplitList[[attribute]][["split_data"]]
    sapply(gg, function(y) length(unique(y[["EWHA_ID"]])))
  })
  
  attributesWithAtLeastTwoScreens <- lapply(attributeSummary, function(x) {
    pp <- names(x[x>=2])
    pp <- pp[pp!=""]
    names(pp) <- pp
    return(pp)
  })
  
  attributesWithAtLeastThreeScreens <- lapply(attributeSummary, function(x) {
    pp <- names(x[x>=3])
    pp <- pp[pp!=""]
    names(pp) <- pp
    return(pp)
  })
  
  return(list(with_two=attributesWithAtLeastTwoScreens,
              with_three=attributesWithAtLeastThreeScreens,
              all_list=allInfoFileSplitList))
}


##########################################################
read_top_bottom_genes <- function(database_matrix,attribute_level_matrix) {
  database_matrix <- read.delim(database_matrix)
  list_of_db <- database_matrix$Database
  list_of_top_genes <- database_matrix$Top_genes
  names(list_of_top_genes) <- list_of_db
  list_of_bottom_genes <- database_matrix$Bottom_genes
  names(list_of_bottom_genes) <- list_of_db
  return(list(bottom=as.list(list_of_bottom_genes), top=as.list(list_of_bottom_genes)))
}


#############################################

generate_gene_overlaps_real <- function(
    allInfoFileSplitList,
    attributeList,
    attributesWithAtLeastTwoScreens,
    percLevels,
    data_dir,database_name, md_file_path) {
  metadata_file <- readRDS(md_file_path)
  fullListOfHyperOverlapsReal <- lapply(attributeList, function(attribute) {
    
    print(attribute)
    
    attributeObject <- allInfoFileSplitList[[attribute]][["split_data"]]
    attributeSpecificNames <- attributesWithAtLeastTwoScreens[[attribute]]
    
    atDiffPercentiles <- lapply(percLevels, function(topPerc) {
      
      print(topPerc)
      
      tempHyperOverlap <- lapply(attributeSpecificNames, function(name) {
        
        print(name)
        
        pp <- deriveOverlapBetweenTop5PercGenes(
          attributeName = name,
          attributeObject = attributeObject,
          attributeType = attribute,
          topPerc = topPerc,
          includeLibInfo = TRUE, metadataFile = metadata_file
        )
        
        hyper_matrix(pp, libInfo = TRUE)
        
      })
      
      hyperOverlap <- do.call(rbind, tempHyperOverlap)
      
      hyperOverlap$adjP <- p.adjust(hyperOverlap$P, method = "BH")
      hyperOverlap$Percentile <- topPerc
      
      hyperOverlap
    })
    
    do.call(rbind, atDiffPercentiles)
  })
  output_path <- file.path(
    data_dir,
    paste0(
      database_name,
      "_fullListOfHyperOverlapsReal.rds"
    )
  )
  saveRDS(fullListOfHyperOverlapsReal, output_path)
  }

#############################################
generate_gene_overlaps_random <- function(
    allInfoFileSplitList,
    attributeList,
    percLevels,
    data_dir,
    database_name,
    nOuter = 20,
    nInner = 100,md_file_path
) {
  metadata_file <- readRDS(md_file_path)
  fullListOfHyperOverlapsRandom <- lapply(attributeList, function(attribute) {
    
    print(attribute)
    
    attributeObject <- allInfoFileSplitList[[attribute]][["split_data"]]
    
    atDiffPercentiles <- lapply(percLevels, function(topPerc) {
      print(topPerc)  
      hh <- lapply(seq_len(nOuter), function(i) {
        print(i)
        randomOverlapsList <- lapply(seq_len(nInner), function(j) {
          print(j)
          pp <- deriveOverlapBetweenTop5PercGenes(
            attributeObject = attributeObject,
            attributeType = attribute,
            isRandom = TRUE,
            topPerc = topPerc,
            includeLibInfo = TRUE, metadataFile = metadata_file
          )
          
          hyper_matrix(pp, libInfo = TRUE)
          
        })
        
        randomOverlaps <- do.call(rbind, randomOverlapsList)
        
        randomOverlaps$adjP <- p.adjust(
          randomOverlaps$P,
          method = "BH"
        )
        
        randomOverlaps$Percentile <- topPerc
        
        randomOverlaps
      })
      
      do.call(rbind, hh)
    })
    
    do.call(rbind, atDiffPercentiles)
  })
  
  output_file <- file.path(
    data_dir,
    paste0(database_name,
           "_fullListOfHyperOverlapsRandom.rds")
  )
  
  saveRDS(fullListOfHyperOverlapsRandom, output_file)
  
  invisible(fullListOfHyperOverlapsRandom)
}


#############################################
generate_score_rank_corr_real <- function(
    allInfoFileSplitList,
    attributeList,
    modes,
    attributesWithAtLeastTwoScreens,
    directions,
    NgenesTop,
    NgenesBottom,
    data_dir,
    database_name,
    metadataFile = NULL
) {
  metadata_file <- readRDS(md_file_path)
  fullListOfOverlapsScoresRank <- lapply(modes, function(mode) {
    
    fullListOfOverlapsScoresRankTopBottomAllList <-
      lapply(directions, function(directionToCheck) {
        
        fullListOfOverlapsScoresRankList <-
          lapply(attributeList, function(attribute) {
            
            print(attribute)
            
            attributeObject <- allInfoFileSplitList[[attribute]][["split_data"]]
            attributeSpecificNames <-
              attributesWithAtLeastTwoScreens[[attribute]]
            
            tempOverlapList <-
              lapply(attributeSpecificNames, function(attributeName) {
                
                print(attributeName)
                
                pp <- deriveCorrelBetweenZscoresOrRank(
                  attributeName = attributeName,
                  attributeObject = attributeObject,
                  attributeType = attribute,
                  mode = mode,
                  directionToCheck = directionToCheck,
                  includeLibInfo = TRUE,
                  NgenesTop = NgenesTop,
                  NgenesBottom = NgenesBottom,
                  md = metadata_file
                )
                
                pp
              })
            
            do.call(rbind, tempOverlapList)
          })
        
        do.call(rbind, fullListOfOverlapsScoresRankList)
      })
    
    do.call(
      rbind,
      fullListOfOverlapsScoresRankTopBottomAllList
    )
  })
  
  output_file <- file.path(
    data_dir,
    paste0(
      database_name,
      "_fullListOfScoresRankOverlapsRealFile.rds"
    )
  )
  
  saveRDS(
    fullListOfOverlapsScoresRank,
    output_file
  )
  }

#############################################
generate_score_rank_corr_random <- function(
    allInfoFileSplitList,
    attributeList,
    modes,
    directions,
    NgenesTop,
    NgenesBottom,
    data_dir,
    database_name,
    attribute_levels = NULL,
    nOuter = 20,
    nInner = 100,md_file_path
) {
  metadata_file <- readRDS(md_file_path)
  
  fullListOfRandomOverlapsScoresRank <- lapply(modes, function(mode) {
    
    fullListOfRandomOverlapsScoresRankTopBottomList <-
      lapply(directions, function(directionToCheck) {
        
        fullListOfRandomOverlapsScoresRankList <-
          lapply(attributeList, function(attribute) {
            
            print(attribute)
            attributeObject <- allInfoFileSplitList[[attribute]][["split_data"]]
            
            if (!is.null(attribute_levels)) {
              
              attribute_level <-
                attribute_levels[
                  attribute_levels$Attribute == attribute,
                ]$Level
              
              higher_level_name <-
                attribute_levels[
                  attribute_levels$Level == (attribute_level + 1),
                ]$Attribute
              
            } else {
              
              attribute_level <- NULL
              higher_level_name <- NULL
              
            }
            list_of_hln <- allInfoFileSplitList[[attribute]][["hln_to_check"]]
            hh <- lapply(seq_len(nOuter), function(yy) {
              
              randomOverlapsList <- lapply(seq_len(nInner), function(x) {
                
                deriveCorrelBetweenZscoresOrRank(
                  attributeObject = attributeObject,
                  attributeType = attribute,
                  attribute_level = attribute_level,
                  higher_level_name = higher_level_name,
                  mode = mode,
                  isRandom = TRUE,
                  includeLibInfo = TRUE,
                  directionToCheck = directionToCheck,
                  NgenesTop = NgenesTop,
                  NgenesBottom = NgenesBottom,md = metadata_file,
                  list_of_hln =list_of_hln
                )
                
              })
              
              tempOverlap <- do.call(
                rbind,
                randomOverlapsList
              )
              
              cbind(
                tempOverlap,
                AttributeType = attribute,
                PermutationNr = yy
              )
              
            })
            
            do.call(rbind, hh)
            
          })
        
        do.call(
          rbind,
          fullListOfRandomOverlapsScoresRankList
        )
        
      })
    
    do.call(
      rbind,
      fullListOfRandomOverlapsScoresRankTopBottomList
    )
    
  })
  
  output_file <- file.path(
    data_dir,
    paste0(
      database_name,
      "_fullListOfScoresRankOverlapsRandomFile.rds"
    )
  )
  
  saveRDS(
    fullListOfRandomOverlapsScoresRank,
    output_file
  )
}