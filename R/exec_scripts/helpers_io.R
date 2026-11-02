#####################################################
load_analysis_inputs <- function(cfg) {
  setwd(cfg$root_dir)
  
  overlap_real <- lapply(cfg$db_names, function(db_name) {
    readRDS(paste0(db_name, "_fullListOfHyperOverlapsReal.rds"))
  })
  names(overlap_real) <- cfg$db_names
  
  overlap_random <- lapply(cfg$db_names, function(db_name) {
    readRDS(paste0(db_name, "_fullListOfHyperOverlapsRandom.rds"))
  })
  names(overlap_random) <- cfg$db_names
  metadata_file_names <- lapply(cfg$db_names, function(db_name) {
   paste0(db_name, "_metadata.rds")
   })
 
  
  list(
    db_lookup = setNames(names(cfg$db_names), cfg$db_names),
    all_info = readRDS("allInfoFileList.rds"),
    intra_db_corr = readRDS("intraDBcorrList.rds"),
    cov_summary = readRDS("all_db_sd_df_db_list.rds"),
    overlap_real = overlap_real,
    overlap_random = overlap_random,
    metadata_list = metadata_file_names, 
    library_metadata = read.delim("library_information.tsv"),
    prediction_of_library_effect = readRDS("prediction_of_library_effect.rds"),
    depMap_metadata = readRDS("depMap_full_metadata.rds"),
    heterogeneity =  read.csv("biological_databases/Heterogeneity.csv"),
    methyl_files = lapply(c(TSS="biological_databases/CCLE_RRBS_tss_CpG_clusters_20181022.txt",
                                            CGI="biological_databases/CCLE_RRBS_cgi_CpG_clusters_20181119.txt",
                                            TSS_1kb="biological_databases/CCLE_RRBS_TSS1kb_20181022.txt"),
                          read.delim, header=T),
    mutational_burden = read.csv("biological_databases/OmicsSomaticMutations.csv"),
    gdsc_DB_files = lapply(list.files(pattern = "GDSC*", recursive = T),read.csv),
    iso_data = readRDS("cell_line_vs_isoform_summary.rds"),
    attribute_levels = read.delim("attribute_level_matrix.tsv")
    )

  
}

#####################################################
combine_intra_db_correlations <- function(intra_db_corr, cfg) {
  combined <- lapply(cfg$metrics, function(metric_name) {
    do.call("rbind", lapply(names(intra_db_corr), function(db_name) {
      cbind(intra_db_corr[[db_name]][[metric_name]], Metric = metric_name, Database = db_name)
    }))
  })
  
  combined_df <- do.call("rbind", combined)
  combined_df <- combined_df[combined_df$attributeType != "HISTOLOGY", ]
  combined_df <- combined_df[!is.na(combined_df$Correlation), ]
  combined_df$Metric <- ifelse(combined_df$Metric == "zcores", "Score", "Rank")
  names(combined_df)[names(combined_df) == "Database"] <- "Dataset"
  combined_df
}

#####################################################
prepare_cov_summary <- function(cov_df,db_lookup) {
  names(cov_df)[names(cov_df) == "Attribute"] <- "attributeName"
  names(cov_df)[names(cov_df) == "Attribute_name"] <- "attributeType"
  cov_df2 <- cov_df[cov_df$attributeType != "HISTOLOGY", ]
  cov_df2 <- cov_df2[cov_df2$Type == "CoV", ]
  cov_df2$Dataset <- db_lookup[cov_df2$Dataset]
  return(cov_df2)
}

#####################################################
build_score_rank_summary <- function(all_info, cfg) {
  list_of_score_vs_rank <- lapply(names(cfg$db_names), function(db_name) {
    if (db_name != "BioGRID") {
    df <- all_info[[db_name]][, c("BCSCORE", "RANK")]
    df$BCSCORE <- as.numeric(df$BCSCORE)
    df <- df[!is.infinite(df$BCSCORE), ]
    df <- aggregate(BCSCORE ~ RANK  , data = df, FUN = mean)
    df$Dataset <- db_name
    df 
    } else {
      df <- all_info[[db_name]][, c("BCSCORE", "RANK","SCORE_TYPE")]
      df$BCSCORE <- as.numeric(df$BCSCORE)
      df <- df[!is.infinite(df$BCSCORE), ]
      df <- aggregate(BCSCORE ~ RANK + SCORE_TYPE , data = df, FUN = mean)
      df$Dataset <- db_name
      df 
    }
  })
  return(list_of_score_vs_rank)
}

#####################################################
prepare_overlap_summary <- function(overlap_real, overlap_random) {
  real_df <- do.call("rbind", lapply(overlap_real, function(db_entry) {
    do.call("rbind", lapply(db_entry, function(attr_df) {
      attr_df$Attribute <- "Same"
      attr_df
    }))
  }))
  
  random_df <- do.call("rbind", lapply(overlap_random, function(db_entry) {
    do.call("rbind", lapply(db_entry, function(attr_df) {
      attr_df$Attribute <- "Different"
      attr_df
    }))
  }))
  
  full_df <- rbind(real_df, random_df)
  full_df$Significant <- ifelse(full_df$adjP < 0.05, "YES", "NO")
  full_df$sameLib <- sapply(full_df$Libraries, addLibraryInfo)
  full_df$negLog10P[is.infinite(full_df$negLog10P)] <- 300
  full_df$Dataset <- sapply(rownames(full_df), function(x) strsplit(x, "\\.")[[1]][1])
  full_df$Percentile <- factor(full_df$Percentile)
  full_df
}

#####################################################
prepare_table_for_mosaic_plot <- function(db,attr,full_DF) {
    DF_01 <- full_DF[full_DF$Percentile=="0.1" & full_DF$Dataset == db & full_DF$AttributeType==attr,]
    DF_01[["Library+Attribute"]]  <- paste(DF_01$sameLib,DF_01$Attribute,  sep = "+")
    print(db)
    print(attr)
    ### They ae associated with the library in the real data 
    t_for_fisher_sa <- table(DF_01[DF_01$Attribute=="Same",c("Significant","sameLib")])
    t_for_fisher_la <- table(DF_01[((DF_01$Attribute=="Same" & DF_01$sameLib=="Different") | (DF_01$Attribute=="Different" & DF_01$sameLib=="Same")), 
                             c("Significant","Library+Attribute")])
    if (all(dim(t_for_fisher_sa)==2)) {
      pVal_sa <- fisher.test(t_for_fisher_sa)$p.value
      pVal_sa <- format(pVal_sa, scientific = TRUE, digits = 3)
      OR_sa <- round(fisher.test(t_for_fisher_sa)$estimate,2)
    } else {
      pVal_sa <- NULL
      OR_sa <- NULL
    }
    if (all(dim(t_for_fisher_la)==2)) {
      pVal_la <- fisher.test(t_for_fisher_la)$p.value
      pVal_la <- format(pVal_la, scientific = TRUE, digits = 3)
      OR_la <- round(fisher.test(t_for_fisher_la)$estimate,2)
    } else {
      pVal_la <- NULL
      OR_la <- NULL
    }
    return(list(t_for_fisher_la=t_for_fisher_la, t_for_fisher_sa=t_for_fisher_sa,
                pVal_sa=pVal_sa, pVal_la=pVal_la,OR_sa=OR_sa,OR_la=OR_la))
}     

#####################################################
correlation_vs_grnas <- function(inputs,input_df,prediction_of_library_effect) {
  df_temp <- input_df[input_df$Genes=="All",]
  df_temp <- df_temp[df_temp$sameLib=="Different",]
  prediction_of_library_effect$Sign <- ifelse(prediction_of_library_effect$p_val < 0.05,"Yes","No")
  prediction_of_library_effect$Direction <- ifelse(prediction_of_library_effect$t_val > 0 ,"Positive","Negative")
  bad_libs_temp <- rownames(prediction_of_library_effect[prediction_of_library_effect$Direction =="Negative" &
                                                           prediction_of_library_effect$Sign=="Yes",])
  good_libs_temp <- rownames(prediction_of_library_effect[prediction_of_library_effect$Direction =="Positive" &
                                                            prediction_of_library_effect$Sign=="Yes",])
  library_metadata <-  inputs$library_metadata
  library_metadata$Correlation  <- ifelse(library_metadata$LIBRARY_ID%in%bad_libs_temp, "Negative",
                                          ifelse(library_metadata$LIBRARY_ID%in%good_libs_temp, "Positive","Neutral"))
  
  library_metadata_filt <- library_metadata[library_metadata$NUM_gene >10000,]
  all_libs <- library_metadata_filt$LIBRARY_ID
  bad_libs <- intersect(bad_libs_temp, all_libs)
  good_libs <- intersect(good_libs_temp, all_libs)
  neutral_libs <- all_libs[!all_libs%in%c(bad_libs,good_libs)]
  
  ### check if # of gRNA is relevant 
  num_list <- lapply(c("NUM_gRNA","NUM_gRNA_PER_GENE"), function(num) {
    corr_list <- lapply(c("Positive","Negative","Neutral"), function(cc){
      tf <- data.frame(Correlation=library_metadata_filt[library_metadata_filt[["Correlation"]]==cc,"Correlation"],
                       Value=library_metadata_filt[library_metadata_filt[["Correlation"]]==cc,num],
                       Metric=num)
    })
    corr_df <- do.call("rbind",corr_list)
    return(corr_df)
  })
  num_corr <- do.call("rbind",num_list) 
  num_corr <- num_corr[num_corr$Correlation != "Neutral",]
  return(num_corr)
}

#####################################################
reformat_inputs <- function(inputs) {
  names(inputs$overlap_real) <- inputs$db_lookup[names(inputs$overlap_real)]
  names(inputs$overlap_random) <- inputs$db_lookup[names(inputs$overlap_random)]
  inputs$metadata_list$BioGRID <- NULL
  inputs$overlap_real$BioGRID <- NULL
  inputs$overlap_random$BioGRID <- NULL
  inputs$all_info$BioGRID <- NULL
  inputs$intra_db_corr$BioGRID <- NULL
  inputs$cov_summary <-  inputs$cov_summary %>%  mutate(across(where(is.character), factor))
  levels(inputs$cov_summary$Dataset) <- inputs$db_lookup[levels(inputs$cov_summary$Dataset)]
  inputs$db_lookup <- inputs$db_lookup[inputs$db_lookup!="BioGRID"]
  return(inputs)
}

#####################################################
prepare_db_summary_table <- function(db_md_file,cfg) {
  db_n <- gsub(pattern = "_metadata.rds",replacement = "",db_md_file)
  md_f <- readRDS(db_md_file)
  total_scr <- data.frame(nrow(md_f),row.names = "Total screens")
  names(total_scr) <- names(cfg$db_names)[cfg$db_names==db_n] 
  summary_df_list <- lapply(cfg[["biological_attributes"]], function(x) {
    x2 <- paste0(gsub(pattern = "_",replacement = " ",x),"s")
    rn <- c(paste0("Nr. of ",tolower(x2), " with >= 2 screens"),
            paste0("Nr. of ",tolower(x2), " with >= 3 screens"))
    two_attr <- sum(table(md_f[[x]]) >1)
    three_attr <- sum(table(md_f[[x]]) >2)
    temp_df <- data.frame(c(two_attr,three_attr),row.names = rn)
    names(temp_df) <- names(cfg$db_names)[cfg$db_names==db_n] 
    return(temp_df)
  })
  summary_df <- do.call("rbind",summary_df_list)
  rbind(total_scr,summary_df)
}

#####################################################
format_db_summary_table <- function(db_files,attributes) {
  gg <- lapply(db_files,prepare_db_summary_table ,cfg=cfg)
  gg2 <-  do.call("cbind",gg)
  kable(gg2, format = "latex", booktabs = TRUE,
        caption = "Summary of databases unsed in this study",
        label = "Tab_1") %>%
    kable_styling(latex_options = "scale_down") %>%
    cat(file = "Table_1.tex")
}
    


#####################################################
prepare_md_long <-function(metadata_file_names, metadata_files,attrs){
  long_df <- map_dfr(metadata_file_names, function(fn) {
    f <- metadata_files[[fn]]
    f1 <- readRDS(f)[,attrs] 
      dplyr::mutate(.data = f1, source = fn)
  }) %>%
    dplyr::distinct() %>%                                   # drop exact duplicate rows
    tidyr::pivot_longer(dplyr::all_of(attrs),
                        names_to = "column", values_to = "value") %>%
    dplyr::mutate(value  = stringr::str_trim(value),
                  column = factor(column, levels = attrs)) %>%
    dplyr::filter(!is.na(value))
  return(long_df)
}
                           
                           
#####################################################
summarize_attributes <- function(long_input) {
  
  # ---------------------------------------------------------------
  # 1) Number of unique values per column, per file
  # ---------------------------------------------------------------
  n_unique <- long_input %>%
    dplyr::group_by(source, column) %>%
    dplyr::summarise(n_unique = dplyr::n_distinct(value), .groups = "drop")
  
  p_unique <- ggplot(n_unique, aes(source, n_unique, fill = column)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_text(aes(label = n_unique), position = position_dodge(0.8),
              vjust =-0.3, size = 3) +
    scale_fill_brewer(palette = "Set2") +
    theme(legend.position = "right", legend.title = element_blank(), axis.title.x = element_blank()) +
    scale_y_break(breaks = c(185, 250))+
    scale_y_break(breaks = c(350, 970), scales = 0.5) + ylim(c(0,1270)) + ylab("")
  
  # ---------------------------------------------------------------
  # 2) Top 5 CANCER values per table, as stacked bars
  # ---------------------------------------------------------------
  top5 <- long_input %>%
    dplyr::filter(column == "CANCER") %>%
    dplyr::count(source, value, name = "n") %>%
    dplyr::group_by(source) %>%
    dplyr::slice_max(n, n = 5, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # One fixed colour per cancer type (named vector -> same colour in every table)
  cancer_levels <- top5 %>%
    dplyr::group_by(value) %>%
    dplyr::summarise(total = sum(n), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total)) %>%
    dplyr::pull(value)
  
  cancer_cols <- setNames(
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(cancer_levels)),
    cancer_levels
  )
  
  top5 <- top5 %>%
    dplyr::mutate(value = factor(value, levels = cancer_levels))
  
  p_top <- ggplot(top5, aes(x = source, y = n, fill = value)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = cancer_cols, name = "Cancer type") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme(legend.position = "right") + xlab("") + ylab("")
  return(list(summ_nr_attr_plot=p_unique,top5_attr_plot=p_top))
}


