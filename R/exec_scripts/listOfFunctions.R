# -----------------------------------------------------------------------------
# General utilities
# -----------------------------------------------------------------------------
#############################################################
DF_factorOrder <- function(DF, factorToOrder, orderingMetric,
                           aggregate = FALSE, aggrStat = NULL) {
  temp_DF <- DF
  
  if (aggregate) {
    if (is.null(aggrStat)) {
      stop("`aggrStat` must be provided when `aggregate = TRUE`.")
    }
    
    stat_fun <- match.fun(aggrStat)
    split_df <- split(DF, DF[[factorToOrder]])
    
    temp_DF <- do.call("rbind", lapply(split_df, function(df_part) {
      df_part[[aggrStat]] <- stat_fun(df_part[[orderingMetric]])
      df_part
    }))
    
    ordering_values <- temp_DF[[aggrStat]]
  } else {
    ordering_values <- temp_DF[[orderingMetric]]
  }
  
  ordered_levels <- unique(
    temp_DF[[factorToOrder]][order(ordering_values, decreasing = TRUE)]
  )
  
  temp_DF[[factorToOrder]] <- factor(
    temp_DF[[factorToOrder]],
    levels = ordered_levels,
    ordered = TRUE
  )
  
  temp_DF
}

#############################################################
reformatDataFramesForFinalAnalysis <- function(kk, mergeBy = NULL, keepAll = TRUE) {
  Reduce(function(x, y) {
    if (identical(mergeBy, "row.names")) {
      merged_df <- merge(x, y, by = 0, all = keepAll)
      rownames(merged_df) <- merged_df$Row.names
      merged_df <- subset(merged_df, select = -Row.names)
    } else {
      merged_df <- merge(x, y, by = mergeBy, all = keepAll)
    }
    
    merged_df
  }, kk)
}

#############################################################
coeffVar <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  if (is.na(mean_x) || mean_x == 0) {
    return(NA_real_)
  }
  
  sd(x, na.rm = TRUE) / mean_x
}

#############################################################
selectAtLeastTriplicates <- function(DF, varName) {
  counts <- table(DF[[varName]])
  names(counts[counts > 2])
}


# -----------------------------------------------------------------------------
# Library helpers
# -----------------------------------------------------------------------------
#############################################################
addLibraryInfo <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(NA_character_)
  }
  
  n_libs <- length(unique(unlist(strsplit(x, ","))))
  ifelse(n_libs == 1, "Same", "Different")
}

#############################################################
selectBasedOnlibs <- function(DF, sameLib = TRUE) {
  colName <- ifelse(any(grepl("Libraries", names(DF))), "Libraries", "libUsed")
  
  if (sameLib) {
    which(sapply(DF[[colName]], function(x) {
      length(unique(unlist(strsplit(x, ",")))) == 1
    }))
  } else {
    which(sapply(DF[[colName]], function(x) {
      length(unique(unlist(strsplit(x, ",")))) == 2
    }))
  }
}

#############################################################
librariesAndDiversityFunc <- function(attributeName, attribute, metadataFile) {
  metadata_subset <- metadataFile[!is.na(metadataFile[[attribute]]), ]
  libUsed <- metadata_subset[metadata_subset[[attribute]] == attributeName, "LIBRARY_ID"]
  
  list(
    libUsed = paste(libUsed, collapse = ","),
    numLibs = length(unique(libUsed))
  )
}

#############################################################
stepWiseGlm <- function(inputDF, dep, indep) {
  model_null <- lm(as.formula(paste(dep, "~ 1")), inputDF)
  model_full <- lm(as.formula(paste(dep, "~", indep)), inputDF)
  
  selected_model <- step(
    model_null,
    scope = formula(model_full),
    direction = "both",
    test = "Chisq",
    data = inputDF,
    trace = 0
  )
  
  selected_terms <- grep(
    "Intercept",
    names(selected_model$coefficients),
    invert = TRUE,
    value = TRUE
  )
  
  if (length(selected_terms) == 0) {
    final_model <- model_null
  } else {
    final_model <- lm(
      as.formula(paste(dep, "~", paste(selected_terms, collapse = "+"))),
      inputDF
    )
  }
  
  summary(final_model)$coefficients
}

#############################################################
create_library_corr_coeff_matrix <- function(DF, colnameLibs = "libUsed",
                                             colNameCorrCoeff = "Correlation") {
  libUsed <- DF[[colnameLibs]]
  
  if (colnameLibs == "libUsed") {
    list_of_libs <- unique(unlist(lapply(libUsed, function(x) unlist(strsplit(x, ",")))))
  } else {
    libUsed <- paste("C", libUsed, sep = "_")
    list_of_libs <- unique(libUsed)
  }
  
  design_df <- data.frame(matrix(0, ncol = length(list_of_libs) + 1, nrow = length(libUsed)))
  names(design_df) <- c(list_of_libs, colNameCorrCoeff)
  
  for (i in seq_along(libUsed)) {
    lib <- libUsed[i]
    
    if (colnameLibs == "libUsed") {
      rel_libs <- unlist(strsplit(lib, ","))
      design_df[i, c(rel_libs, colNameCorrCoeff)] <- c(rep(1, length(rel_libs)), DF[i, colNameCorrCoeff])
    } else {
      design_df[i, c(lib, colNameCorrCoeff)] <- c(1, DF[i, colNameCorrCoeff])
    }
  }
  
  indepVar <- paste(list_of_libs, collapse = "+")
  model_summary <- stepWiseGlm(design_df, dep = colNameCorrCoeff, indep = indepVar)
  model_summ_df <- data.frame(model_summary)[-1, , drop = FALSE]
  names(model_summ_df) <- c("Est", "Std_err", "t_val", "p_val")
  
  if (colnameLibs != "libUsed") {
    rownames(model_summ_df) <- gsub("C_", "", rownames(model_summ_df))
  }
  
  model_summ_df
}


# -----------------------------------------------------------------------------
# Overlap and correlation helpers
# -----------------------------------------------------------------------------
#############################################################
hyper_matrix <- function(geneListsAndAttributes, background = 10000, libInfo = FALSE) {
  if (is.null(geneListsAndAttributes)) {
    return(NULL)
  }
  
  gene.list <- geneListsAndAttributes[[1]]
  
  if (libInfo) {
    lib <- geneListsAndAttributes[[2]]
    combinationLibs <- expand.grid(lib, lib)
    screensUsed <- geneListsAndAttributes[[3]]
    attrType <- geneListsAndAttributes[[4]]
  } else {
    screensUsed <- geneListsAndAttributes[[2]]
    attrType <- geneListsAndAttributes[[3]]
  }
  
  combination <- expand.grid(names(gene.list), names(gene.list))
  combination$values <- NA
  combination <- reshape(combination, idvar = "Var1", timevar = "Var2", direction = "wide")
  rownames(combination) <- combination$Var1
  combination <- combination[, -1]
  colnames(combination) <- gsub("values.", "", colnames(combination), fixed = TRUE)
  
  attrNames <- paste(names(gene.list), collapse = ",")
  
  for (i in colnames(combination)) {
    for (j in rownames(combination)) {
      combination[j, i] <- length(intersect(gene.list[[j]], gene.list[[i]]))
    }
  }
  
  pval <- c()
  libComb <- c()
  k <- 1
  klib <- 1
  
  for (m in seq_along(gene.list)) {
    for (n in seq_along(gene.list)) {
      if (n > m) {
        pval[k] <- phyper(
          combination[m, n] - 1,
          length(gene.list[[m]]),
          background - length(gene.list[[m]]),
          length(gene.list[[n]]),
          lower.tail = FALSE
        )
        
        if (libInfo) {
          libComb[k] <- paste(unlist(combinationLibs[klib, ]), collapse = ",")
        }
        
        k <- k + 1
      }
      klib <- klib + 1
    }
  }
  
  if (libInfo) {
    return(data.frame(
      P = pval,
      negLog10P = -log10(pval),
      AttributeNames = attrNames,
      AttributeType = attrType,
      Libraries = libComb,
      Screens = screensUsed
    ))
  }
  
  data.frame(
    P = pval,
    negLog10P = -log10(pval),
    AttributeNames = attrNames,
    AttributeType = attrType,
    Screens = screensUsed
  )
}

#############################################################
generateBoxPlotOfHyGeomForOverGenes <- function(listOfListOfTop5PercGenes, libInfo = FALSE) {
  overlap_tables <- lapply(listOfListOfTop5PercGenes, hyper_matrix, libInfo = libInfo)
  overlap_df <- do.call("rbind", overlap_tables)
  
  overlap_df$adjP <- p.adjust(overlap_df$P, method = "BH")
  overlap_df$adjP[overlap_df$adjP == 0] <- 1e-300
  overlap_df$logAdjP <- -log10(overlap_df$adjP)
  overlap_df$AttributeNames <- with(overlap_df, reorder(AttributeNames, logAdjP, mean))
  
  plot_obj <- ggplot(overlap_df, aes(AttributeNames, logAdjP, fill = AttributeNames)) +
    geom_boxplot() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    geom_hline(yintercept = 1.3, color = "red")
  
  list(overlapHypGeomData = overlap_df, overlapHypGeomPlot = plot_obj)
}

#############################################################
generateBoxPlotOfCorrSameCell <- function(listOfCorrVector) {
  corr_df <- do.call("rbind", listOfCorrVector)
  corr_df <- na.omit(corr_df)
  corr_df$cellLine <- with(corr_df, reorder(cellLine, Correlation, mean))
  
  ggplot(corr_df, aes(cellLine, Correlation, fill = cellLine)) +
    geom_boxplot() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none"
    )
}

#############################################################
hypGeomCorrValueAndLibraryCorrelation <- function(pastedListOfHypGeomOverlapPvalues, mdFile) {
  pastedListOfHypGeomOverlapPvalues$Sign <- ifelse(
    pastedListOfHypGeomOverlapPvalues$adjP < 0.05,
    "YES",
    "NO"
  )
  
  listOfLibOccurrences <- lapply(c("YES", "NO"), function(sign_label) {
    tempLibs <- pastedListOfHypGeomOverlapPvalues[
      pastedListOfHypGeomOverlapPvalues$Sign == sign_label,
      "Libraries"
    ]
    as.data.frame(table(sapply(tempLibs, function(x) unlist(strsplit(x, ",")))))
  })
  
  merged_df <- merge(listOfLibOccurrences[[1]], listOfLibOccurrences[[2]], by = "Var1", all = TRUE)
  merged_df[is.na(merged_df)] <- 0
  rownames(merged_df) <- merged_df$Var1
  merged_df <- merged_df[, -1]
  names(merged_df) <- c("YES", "NO")
  
  merged_df
}

#############################################################
analyzeOutliers <- function(input_df = input_df, split_var = split_var,
                            variable_to_quant = variable_to_quant) {
  groups_to_compare <- names(table(input_df[[split_var]]))
  names(groups_to_compare) <- groups_to_compare
  
  final_outliers_attr <- lapply(groups_to_compare, function(group_name) {
    input_df$Aux <- factor(ifelse(input_df[[split_var]] == group_name, group_name, paste0("No_", group_name)))
    outliers <- identify_outliers(input_df, variable = variable_to_quant)
    
    outlier_counts <- data.frame(table(outliers[["Aux"]]))
    total_counts <- data.frame(table(input_df[["Aux"]]))
    merged_counts <- merge(outlier_counts, total_counts, by = "Var1")
    
    rownames(merged_counts) <- merged_counts$Var1
    names(merged_counts) <- c("Var1", "Outliers", "All")
    merged_counts <- subset(merged_counts, select = -Var1)
    
    no_group_row <- grep("No_", rownames(merged_counts))
    if (no_group_row == 1) {
      fisher.test(merged_counts, alternative = "l")$p.value
    } else {
      fisher.test(merged_counts, alternative = "g")$p.value
    }
  })
  
  temp_df <- data.frame(t(data.frame(final_outliers_attr)))
  names(temp_df) <- "Fisher_p"
  temp_df$Variable <- split_var
  temp_df$Significant <- temp_df$Fisher_p < 0.05
  temp_df
}

#############################################################
select_features_from_corr_data <- function(feature_type, feature_names, is_grep = FALSE,
                                           action = "remove", corr_df) {
  if (action == "remove") {
    temp_df <- corr_df
    
    if (is_grep) {
      for (feat_id in feature_names) {
        temp_df <- temp_df[grep(feat_id, temp_df[[feature_type]], invert = TRUE), ]
      }
    } else {
      temp_df <- temp_df[!temp_df[[feature_type]] %in% feature_names, ]
    }
    
    return(temp_df)
  }
  
  if (is_grep) {
    temp_df <- data.frame()
    for (feat_id in feature_names) {
      temp_df <- rbind(temp_df, corr_df[grep(feat_id, corr_df[[feature_type]]), ])
    }
    return(temp_df)
  }
  
  corr_df[corr_df[[feature_type]] %in% feature_names, ]
}


# -----------------------------------------------------------------------------
# Variance helpers
# -----------------------------------------------------------------------------
#############################################################
generateRankVarianceMatrix <- function(allInfoFile = allInfoFile, variableNames = variableNames,
                                       subsetTest = FALSE, datasetName = datasetName) {
  all_info_clean <- allInfoFile[grep("-", allInfoFile$HGNC, invert = TRUE), ]
  
  genesToTest <- if (subsetTest) {
    unique(all_info_clean$HGNC)[1:10]
  } else {
    unique(all_info_clean$HGNC)
  }
  
  gene_level_tables <- lapply(genesToTest, function(gene) {
    gene_df <- allInfoFile[allInfoFile$HGNC == gene, ]
    
    attribute_tables <- lapply(variableNames, function(variable) {
      valid_instances <- selectAtLeastTriplicates(gene_df, variable)
      filtered_df <- gene_df[gene_df[[variable]] %in% valid_instances, ]
      split_df <- split(filtered_df, filtered_df[[variable]])
      
      per_attribute_tables <- lapply(names(split_df), function(varName) {
        var_split_df <- split_df[[varName]]
        
        rank_values <- as.numeric(var_split_df[["RANK"]])
        score_values <- as.numeric(var_split_df[["BCSCORE"]])
        scaled_scores <- rescale(score_values, c(min(rank_values), max(rank_values)))
        
        rbind(
          data.frame(Value = coeffVar(rank_values), Type = "CoV", Metric = "Rank",
                     Attribute_name = variable, Attribute = varName, Dataset = datasetName),
          data.frame(Value = coeffVar(scaled_scores), Type = "CoV", Metric = "Score",
                     Attribute_name = variable, Attribute = varName, Dataset = datasetName),
          data.frame(Value = sd(rank_values), Type = "SD", Metric = "Rank",
                     Attribute_name = variable, Attribute = varName, Dataset = datasetName),
          data.frame(Value = sd(scaled_scores), Type = "SD", Metric = "Score",
                     Attribute_name = variable, Attribute = varName, Dataset = datasetName)
        )
      })
      
      do.call("rbind", per_attribute_tables)
    })
    
    attribute_df <- do.call("rbind", attribute_tables)
    if (!is.null(attribute_df)) {
      attribute_df$Gene <- gene
      attribute_df
    }
  })
  
  gene_level_tables[sapply(gene_level_tables, is.null)] <- NULL
  combined_df <- do.call("rbind", gene_level_tables)
  
  aggregate(
    Value ~ Attribute + Attribute_name + Metric + Type + Dataset,
    data = combined_df,
    FUN = mean
  )
}

#############################################################
analyze_coeff_of_variation <- function(dbName, rootDir) {
  metric_tables <- lapply(c("BCSCORE", "RANK"), function(metric) {
    overlap_tables <- lapply(c("Real", "Random"), function(overlap_type) {
      file_name <- paste(dbName, overlap_type, metric, "allDataNeeded.rds", sep = "_")
      file_path <- file.path(rootDir, "data", dbName, file_name)
      data_rds <- readRDS(file_path)
      
      do.call("rbind", lapply(names(data_rds), function(attribute_name) {
        df <- data.frame(data_rds[[attribute_name]])
        names(df) <- "Coef_of_var"
        df$Attribute <- attribute_name
        df$Metric <- metric
        df$Overlap <- overlap_type
        df
      }))
    })
    
    do.call("rbind", overlap_tables)
  })
  
  do.call("rbind", metric_tables)
}


# -----------------------------------------------------------------------------
# Genomic range helpers
# -----------------------------------------------------------------------------
#############################################################
overlapLibAndIsofExprcCoord <- function(x) {
  exonCoord[unlist(exonCoord$Parent) == paste("transcript", x, sep = ":"), c("seqid", "start", "end")]
}

#############################################################
extractExonCoord <- function(x, parent) {
  transcript_id <- paste("transcript", x, sep = ":")
  exonCoord[parent == transcript_id, c("seqid", "start", "end")]
}

#############################################################
range_overlap_mosaic_plot <- function(input_ranges) {
  ranges_overlap <- round(mean(sapply(random_libs_ranges, function(rlr) {
    length(findOverlaps(input_ranges, rlr))
  })), 0)
  
  ranges_size <- length(random_libs_ranges[[1]])
  bad_lib_overlap <- length(findOverlaps(input_ranges, bad_lib_ranges))
  bad_lib_size <- length(bad_lib_ranges)
  
  mosaicplot(
    matrix(c(ranges_overlap, ranges_size, bad_lib_overlap, bad_lib_size), nrow = 2),
    xlab = "Significant overlap",
    ylab = "Library"
  )
}
