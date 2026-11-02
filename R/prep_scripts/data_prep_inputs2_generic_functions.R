##########################################################
generate_all_info_list <- function(dbNames,data_dir) { 
  allInfoList <- list()
  for (dbName in dbNames) {
    all_info_file <- file.path(data_dir,paste0(dbName,"_compact_with_all_info_noNA.rds"))
    all_info_data <- readRDS(all_info_file)
    if (any(names(all_info_data)=="SCORE")) {
      names(all_info_data)[names(all_info_data)=="SCORE"] <- "BCSCORE" 
    }
    all_info_data$BCSCORE <- as.numeric(all_info_data$BCSCORE)
    all_info_data <- all_info_data[all_info_data$BCSCORE  >-20,]
    all_info_data$BCSCORE <- rescale(all_info_data$BCSCORE,to = c(0,1))
    allInfoList[[dbName]] <-  all_info_data
  }
  file_name <- file.path(data_dir,"allInfoFileList.rds")
  saveRDS(allInfoList,file_name)
}


##########################################################
generate_all_db_sd_df_db_list <- function(ds,data_dir,attribute_list) { 
    input_compact_with_all_info_noNA <- file.path(data_dir,
                                                  paste0(ds, "_compact_with_all_info_noNA.rds"))
    allInfoFile <- readRDS(input_compact_with_all_info_noNA)
    allDataNeeded <- generateRankVarianceMatrix(allInfoFile = allInfoFile, datasetName =  ds,
                                                variableNames=attribute_list,
                                                subsetTest=F)
    file_name <- file.path(data_dir, paste0(ds,"_sd_df.rds"))
    saveRDS(allDataNeeded, file_name)
}


##########################################################
generate_intraDB_list <- function(dbNames,data_dir,metrics) { 
  intraDBcorrList <- list()
  for (dbName in dbNames) {
    fullListOfScoresRankOverlapsRealFilePath <- file.path(data_dir,
                                                          paste0(dbName,"_fullListOfScoresRankOverlapsRealFile.rds"))
    fullListOfScoresRankOverlapsRealFile <- readRDS(fullListOfScoresRankOverlapsRealFilePath)
    temp_list <- fullListOfScoresRankOverlapsRealFile
    names(temp_list) <- metrics
    filtered_list <- lapply(metrics, function(metric) {
      t_df <- temp_list[[metric]]
      t_df[["Genes"]] <- sapply(rownames(t_df), function(y) unlist(strsplit(y, "\\."))[1])
      t_df[["sameLib"]] <- sapply(t_df[["libUsed"]], addLibraryInfo)
      return(t_df)
    })
    intraDBcorrList[[dbName]] <- filtered_list
  }
  file_name <- file.path(data_dir,"intraDBcorrList.rds")
  saveRDS(intraDBcorrList,file_name)
}



##########################################################

create_library_corr_coeff_matrix <- function(intraDB_list_file, dbNames, db_specific_name,
                                             metrics, colnameLibs = "libUsed",
                                             colNameCorrCoeff = "Correlation") {
  intraDB_list <- readRDS(intraDB_list_file)
all_db_corr_df_metric_list <- lapply(metrics, function(metric) {
  all_metric_list <-  lapply(dbNames, function(name) {
    pp <- cbind(intraDB_list[[name]][[metric]], Metric=metric, Database=name)
    return(pp)
  })
  pp2 <- do.call("rbind",all_metric_list)
  return(pp2)
}) 
corr_all_db_all_metrics <-  do.call("rbind",all_db_corr_df_metric_list)
corr_all_db_all_metrics_sc_with_dup <- corr_all_db_all_metrics
corr_all_db_all_metrics_sc_with_dup <- do.call("rbind",all_db_corr_df_metric_list)
corr_all_db_all_metrics_sc_with_dup$Metric <- ifelse(corr_all_db_all_metrics_sc_with_dup$Metric=="zcores", 
                                                     "Score","Rank")
corr_all_db_all_metrics_sc_with_dup <- corr_all_db_all_metrics_sc_with_dup[corr_all_db_all_metrics_sc_with_dup$Metric =="Score",]
names(corr_all_db_all_metrics_sc_with_dup)[names(corr_all_db_all_metrics_sc_with_dup)=="Database"]  <- "Dataset"
split_by_db <- split(corr_all_db_all_metrics_sc_with_dup,corr_all_db_all_metrics_sc_with_dup$Dataset)
DF <- split_by_db[[db_specific_name]]
DF <- DF[DF$Genes=="All",]
DF <- DF[DF$sameLib=="Different",]
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
    print(i)
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
  file_name <- file.path(data_dir,"prediction_of_library_effect.rds")
  saveRDS(model_summ_df,file_name)
}