##########################################################
prepare_cell_map_table <- function(depMap_dir) {
  depMap_model_file <- file.path(depMap_dir,"Model.csv")
  MD <- read.csv(depMap_model_file)
  df <- data.frame(CELL_LINE=gsub("NIH","",MD$StrippedCellLineName),TISSUE=MD$SampleCollectionSite,
  CANCER=MD$OncotreePrimaryDisease, SUBTYPE=MD$OncotreeSubtype)
  return(df)
}
 
 
##########################################################
prepare_biogrid_data <- function(biogrid_dir, cell_line_map,
                                 data_dir,biogrid_score_subset) { 
#### Generate compact info file
screen_files <- list.files(path = biogrid_dir, pattern = "*screen.tab.txt", 
                           full.names = T)
metadata_biogrid_file <- list.files(path = biogrid_dir, pattern = "*index.tab.txt",
                                    full.names = T)
metadata_biogrid <- read.delim(metadata_biogrid_file)
metadata_biogrid$CELL_LINE <- gsub("-","",metadata_biogrid$CELL_LINE)
metadata <- merge(cell_line_map,metadata_biogrid,by="CELL_LINE")
metadata <- metadata[metadata$LIBRARY_METHODOLOGY == "Knockout" & 
                       metadata$SCREEN_TYPE =="Negative Selection", ]

invisible(lapply(names(biogrid_score_subset),function(score_cat) {
  score_types <- biogrid_score_subset[[score_cat]]
  screen_list  <- metadata[metadata$SCORE.1_TYPE%in%score_types,"X.SCREEN_ID"]
  screenFilesToUse <- sapply(screen_list, function(s) {
    file.path(biogrid_dir,
              paste0("BIOGRID-ORCS-SCREEN_",s,"-2.0.18.screen.tab.txt"))
  })  
    names(screenFilesToUse) <- screen_list

### prepare and save metadata
metadata_to_save <- metadata[metadata$X.SCREEN_ID%in%screen_list,]
metadata_to_save$LIBRARY <- factor(metadata_to_save$LIBRARY)
metadata_to_save$LIBRARY_ID <- paste0("lib",as.integer(metadata_to_save$LIBRARY))
metadata_to_save <- metadata_to_save[,c("X.SCREEN_ID","CELL_LINE","TISSUE","CANCER","SUBTYPE",
                                        "LIBRARY_ID")]
names(metadata_to_save)[1] <- "EWHA_ID"
md_file_name <- file.path(data_dir,paste0(score_cat,"_metadata.rds"))
saveRDS(metadata_to_save, md_file_name)
###

screenDFList <- lapply(screen_list, function(x1){
  x <- screenFilesToUse[as.character(x1)]
  print(x)
  print(x1)
  y <- read.delim(x)
  y$RANK <- 1:nrow(y)
  screenName <- as.integer(unique(y["X.SCREEN_ID"]))
  MD <- metadata[metadata[["X.SCREEN_ID"]]==screenName,]
  DF <- data.frame(HGNC=y[["OFFICIAL_SYMBOL"]], BCSCORE=y[["SCORE.1"]], EWHA_ID=y[["X.SCREEN_ID"]],
                   RANK=y[["RANK"]],CELL_LINE=MD$CELL_LINE, LIBRARY=MD$LIBRARY,
                   TISSUE=MD$TISSUE, CANCER=MD$CANCER, SUBTYPE=MD$SUBTYPE, 
                   SCORE_TYPE=score_cat)
  return(DF)
})

biogrid_compact_with_all_info <- do.call("rbind",screenDFList)
biogrid_compact_with_all_info_noNA  <- na.omit(biogrid_compact_with_all_info)
biogrid_compact_with_all_info_noNA <- data.table(biogrid_compact_with_all_info_noNA)
file_name <- file.path(data_dir,paste0(score_cat,"_compact_with_all_info_noNA.rds"))
saveRDS(biogrid_compact_with_all_info_noNA, file_name)
  })
)
}


##########################################################
main_reformat_for_depmap <- function(y,metadata) {
  x <- unlist(y)
  colData <- x[!is.na(x) & !is.infinite(x) & names(x)!="screenName" & names(x)!="temp_id"]
  colData <- colData[order(as.numeric(colData))]
  
  screenName <- y[["screenName"]]
  MD <- metadata[metadata$ModelID==screenName,,drop=F]
  if (nrow(MD) > 0) {
    test1 <- data.frame(cbind(HGNC=names(colData), BCSCORE=colData, EWHA_ID=y[["temp_id"]],
                              CELL_LINE=gsub("NIH","",MD$StrippedCellLineName),TISSUE=MD$SampleCollectionSite,
                              CANCER=MD$OncotreePrimaryDisease, SUBTYPE=MD$OncotreeSubtype))
    test1$RANK <-  1:nrow(test1)
  }else {
    test1 <- NULL
  }
  return(test1)
}


##########################################################
prepare_depMap_data <- function(depMap_dir,data_dir,depMap_screens_to_use) { 
depMap_screen_file <-  file.path(depMap_dir,"ScreenGeneEffect.csv")
depMap_model_file <- file.path(depMap_dir,"Model.csv")
depMap_scores <- read.csv(depMap_screen_file, header = T, row.names = 1)
rownames(depMap_scores) <- sapply(rownames(depMap_scores), function(x) {
  p<- gsub("SC","ACH",x)
  return(p)
})
depMap_scores <- depMap_scores[depMap_screens_to_use,]
names(depMap_scores) <- sapply(names(depMap_scores), function(p) {
  unlist(strsplit(p, "\\."))[1]
})
depMap_scores$screenName <-  sapply(rownames(depMap_scores), function(p) {
  unlist(strsplit(p, "\\."))[1]
})
depMap_scores$temp_id <- rownames(depMap_scores)

metadata <- read.csv(depMap_model_file, header = T)
ss2 <-  apply(depMap_scores,1,main_reformat_for_depmap, metadata=metadata,simplify = F)
ss2[sapply(ss2, is.null)] <- NULL
depMap_compact_with_all_info <- do.call("rbind",ss2)
depMap_compact_with_all_info_noNA  <- na.omit(depMap_compact_with_all_info)
depMap_compact_with_all_info_noNA <- data.table(depMap_compact_with_all_info_noNA)
rds_file_name_all <- file.path(data_dir,"depMap_full_compact_with_all_info_noNA.rds")
saveRDS(depMap_compact_with_all_info_noNA,rds_file_name_all)
}


##########################################################
prepare_depMap_metadata <- function(depMap_dir,data_dir) { 
depMap_screen_file <-  file.path(depMap_dir,"ScreenGeneEffect.csv")
depMap_model_file <- file.path(depMap_dir,"Model.csv")
metadata <- read.csv(depMap_model_file, header = T)
screens <- read.table(pipe(paste0("cut -d \",\" -f 1 ", depMap_screen_file)))[,1]
screens <- gsub("SC","ACH",screens)
screen_form_list <- lapply(screens, function(x) {
  model_id <- unlist(strsplit(x,"\\."))[1]
  model_lib_temp <- unlist(strsplit(x,"\\."))[2]
  model_lib <- unlist(strsplit(model_lib_temp,"0"))[1]
  model_metadata <- metadata[metadata$ModelID==model_id,]
  if (nrow(model_metadata) >0) {
  temp_df <- data.frame(EWHA_ID=x,CELL_LINE=gsub("NIH","",model_metadata$StrippedCellLineName),
                                    TISSUE=model_metadata$SampleCollectionSite,
                                    CANCER=model_metadata$OncotreePrimaryDisease, 
                                    SUBTYPE=model_metadata$OncotreeSubtype,
                                    LIBRARY_ID=model_lib)
  return(temp_df)
}
})
screen_form_list[sapply(screen_form_list, is.null)]<- NULL
screen_form <- do.call("rbind",screen_form_list)
md_file_name <- file.path(data_dir,"depMap_full_metadata.rds")
saveRDS(screen_form,md_file_name)
}


