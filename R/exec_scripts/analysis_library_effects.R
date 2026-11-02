
####### check same cell+same lib  vs same lib + diff cell 
check_ss_versus_ds <- function(db_n) {
  if (db_n=="iCSDB") {
    screen_info_df <- read.delim("../../information/screen_information.tsv")
    multi_screen_ind <- table(screen_info_df$SOURCE_ID) >1 
    multi_screens <- screen_info_df$SOURCE_ID[multi_screen_ind]
    pp <- lapply(multi_screens, function(s) {
      cl <- screen_info_df[screen_info_df$SOURCE_ID==s,"LIBRARY_ID"]
      cl_2 <- table(cl) > 1
      cl_2 <- factor(cl_2, levels = c("TRUE", "FALSE"))
      return(table(cl_2))
    })
    
  } else if (db_n=="biogrid") {
    screen_info_df <- read.delim("biogrid/BIOGRID-ORCS-SCREEN_INDEX-2.0.18.index.tab.txt")
    multi_screen_ind <- names(which(table(screen_info_df$AUTHOR) >1))
    multi_screens <- screen_info_df[screen_info_df$AUTHOR%in%multi_screen_ind,]
    attrs <- c("Same library" = "LIBRARY",
               "Same cell line" = "CELL_LINE")
    pp3 <-lapply(names(attrs), function(nm) {
      nn <- attrs[nm]
      pp <- lapply(multi_screen_ind, function(s) { 
        cl <- screen_info_df[screen_info_df$AUTHOR==s,nn]
        cl_2 <- table(cl) > 1
        cl_2 <- factor(cl_2, levels = c("TRUE", "FALSE"))
        return(table(cl_2))
      })
      pp2 <- data.frame(do.call("rbind",pp))
      return(data.frame("Dl"=sum(pp2[["FALSE."]]),"Sl"=sum(pp2[["TRUE."]]), row.names = nm))
    })
    yy <- data.frame(do.call("rbind",pp3))
    names(yy) <- c("Different lab", "Same lab")
    yy <- as.table(as.matrix(yy))
    plot_name <- "Fig_3c.pdf"
    pval_rf1 <- sprintf("%.2e", fisher.test(yy)$p.value)
    plot_title <- substitute(
      "Fisher p = " ~ base %*% 10^exponent, 
      list(
        base = strsplit(format(pval_rf1, scientific = TRUE), "e")[[1]][1],
        exponent = strsplit(format(pval_rf1, scientific = TRUE), "e")[[1]][2]
      )
    )
    
    save_plot_pdf(cfg, plot_name,mosaicplot(yy,main = plot_title,
                                            col = c("firebrick", "steelblue"),cex.axis = 0.6),6,5)
  }
}

####### lollipop for the two-step library regression 
two_step_lib_lolli <- function(df) {
  df$p_val[df$p_val == 0] <- .Machine$double.xmin   # p = 0 -> largest dot
  df$neg_log10_p <- -log10(df$p_val)
  
  l_plot <- ggplot(df, aes(x = reorder(LIBRARY_NAME, t_val), y = t_val)) +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_segment(aes(xend = LIBRARY_NAME, yend = 0), color = "grey40", linewidth = 0.8) +
    geom_point(aes(size = neg_log10_p, fill = t_val > 0),
               shape = 21, color = "black") +
    scale_size_continuous(range = c(3, 12), name = expression(-log[10](p))) +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                      guide = "none") +
    coord_flip() +
    labs(x = NULL, y = "t-value") +
    theme_minimal(base_size = 14)
  
  save_plot_pdf(cfg, "Fig_5a.pdf",l_plot,8, 6)
}

############################################################################
prepare_loolipop_plot <- function(df) {
  df <- pivot_longer(df, cols = c(Frac_diff, Mean_frac), 
                     names_to = "metric", values_to = "value")
  df$plot_value <- ifelse(df$metric == "Frac_diff", 
                          -df$value, 
                          df$value)
  df <- df %>%  mutate(across(where(is.character), factor))
  df <- df %>%
    mutate(
      Top_percentile = factor(Top_percentile,
                              levels = sort(unique(Top_percentile))),
      metric         = factor(metric)   # ensure discrete shape mapping
    )
  
  datasets  <- sort(unique(df$Dataset))
  n_ds      <- length(datasets)
  
  group_spacing <- 3
  dodge_width   <- 1.2
  
  offsets <- seq(-dodge_width/2, dodge_width/2, length.out = n_ds)
  names(offsets) <- datasets
  
  df <- df %>%
    mutate(
      x_num    = as.numeric(Top_percentile) * group_spacing,
      x_dodged = x_num + offsets[Dataset]
    )
  
  group_positions <- sort(unique(df$x_num))
  boundaries <- (head(group_positions, -1) + tail(group_positions, -1)) / 2
  
  f_plot <- ggplot(df, aes(x = x_dodged, y = plot_value, color = Dataset)) +
    geom_vline(xintercept = boundaries, color = "grey85", linewidth = 0.25) +
    geom_segment(aes(x = x_dodged, xend = x_dodged, y = 0, yend = plot_value),
                 linewidth = 0.9) +
    geom_point(aes(shape = metric), size = 2) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.25) +
    facet_wrap(~ Attribute) +
    scale_x_continuous(
      breaks = group_positions,
      labels = levels(df$Top_percentile)
    ) +
    scale_y_continuous(labels = function(v) abs(v)) +   # left side reads as positive too
    coord_flip() +                                       # makes the divergence horizontal
    labs(
      x = "Top percentile",
      color = "Dataset",
      shape = "Metric",
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title.x = element_blank(),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8,hjust=0.5),
      legend.position = "bottom",
      legend.title.position = "top",
      legend.spacing.x = unit(2.5, "cm"),
      legend.background = element_rect(
        color = "black",       # Border line color
        linewidth = 0.5,       # Border line thickness
        fill = "grey98"        # Background color inside the box (use "transparent" if you want it clear)
      )
    ) +
    guides(
      color = guide_legend(
        ncol = 2,       # Force exactly 2 columns
        byrow = TRUE    # Fill row-by-row (Group 1 & 2 on top, 3 & 4 on bottom)
      ),
      shape = guide_legend(
        ncol = 1,       # Force exactly 2 columns
        byrow = TRUE    # Fill row-by-row (Group 1 & 2 on top, 3 & 4 on bottom)
      )
    )
  return(f_plot)
}


############################################################################
prepare_summary_of_overlap_fisher_fun <- function(cfg,full_overlap_df){
summary_of_overlap_fisher_sa <- data.frame()
summary_of_overlap_fisher_la <- data.frame()

n <- 0
n1 <- 0
db_to_check  <- names(cfg$db_names)[!names(cfg$db_names)%in%c("BioGRID")]
for (db in db_to_check) {
  for (attr in cfg$biological_attributes) {
    t_p_OR_list <-  prepare_table_for_mosaic_plot(db,attr,full_overlap_df)
    print(t_p_OR_list$t_for_fisher_la)
    if (!is.null(t_p_OR_list$pVal_la)) {
      if (attr=="CELL_LINE" & db=="iCSDB")  {
        plot_name_la <- "Fig_S5.pdf"
      } else if (attr=="CELL_LINE" & db=="BioGRID_exp")  {
        plot_name_la <- "Fig_S4.pdf"
      }else {
        plot_name_la <- paste("Mosaic_LA",db,attr,"pdf",sep = ".")
      }
      save_plot_pdf(cfg, plot_name_la,mosaicplot(t_p_OR_list$t_for_fisher_la,
                                                 main = paste0("Fisher p = ",t_p_OR_list$pVal_la,"; OR = ",t_p_OR_list$OR_la),
                                                 xlab = "Significant overlap",
                                                 ylab = "Library + Attribute",col = c("firebrick", "steelblue")),6,6)
      n <- n+1
      summary_of_overlap_fisher_la[n,"Attribute"] <- attr
      summary_of_overlap_fisher_la[n,"DB"] <- db
      summary_of_overlap_fisher_la[n,"p_value"] <- as.numeric(t_p_OR_list$pVal_la)
      summary_of_overlap_fisher_la[n,"OR"] <- as.numeric(t_p_OR_list$OR_la)
    }
    if (!is.null(t_p_OR_list$pVal_sa)) {
      n1 <- n1+1
      summary_of_overlap_fisher_sa[n1,"Attribute"] <- attr
      summary_of_overlap_fisher_sa[n1,"DB"] <- db
      summary_of_overlap_fisher_sa[n1,"p-value"] <- as.numeric(t_p_OR_list$pVal_sa)
      summary_of_overlap_fisher_sa[n1,"Odds ratio"] <- as.numeric(t_p_OR_list$OR_sa)
    }
  }
}
return(list(summary_of_overlap_fisher_sa=summary_of_overlap_fisher_sa,
            summary_of_overlap_fisher_la=summary_of_overlap_fisher_la))
}


############################################################################
overlap_ratios_fun <- function(db,full_overlap_df) {
  print(db)
  dataset_df <- full_overlap_df[full_overlap_df$Dataset == db, ]
  percentile_split <- do.call("rbind", lapply(levels(dataset_df$Percentile), percentile_split_fun,
                                              dataset_df=dataset_df,db=db))
  percentile_split$Dataset <- db
  return(percentile_split)
}

############################################################################
percentile_split_fun <-  function(percentile_value,dataset_df,db) {
  print(percentile_value)
  percentile_df <- dataset_df[dataset_df$Percentile == percentile_value, ]
  split_by_attribute <- split(percentile_df, percentile_df$AttributeType)
  temp_df <- do.call("rbind", lapply(split_by_attribute,split_by_attribute_fun,
                                     db_name=db,percentile_value=percentile_value))
  return(temp_df)
}

############################################################################
split_by_attribute_fun <-  function(attribute_df,percentile_value,db_name) {
  split_by_overlap <- split(attribute_df, attribute_df$Attribute)
  
  sign_fractions <- sapply(split_by_overlap, sign_fractions_fun)
  diff_ov <- c(diff=abs(unname(diff(sign_fractions))),
               m=mean(sign_fractions))
  final_df <- data.frame(Top_percentile=percentile_value,
                         Dataset=db_name,
                         Frac_diff=diff_ov["diff"],
                         Mean_frac=diff_ov["m"],
                         Attribute=unique(attribute_df$AttributeType))
  return(final_df)
}

############################################################################
sign_fractions_fun <- function(x)  {
  freq_df <- data.frame(table(x$Significant) / nrow(x))
  freq_df$AttributeType <- unique(x$AttributeType)
  names(freq_df) <- c("Significant", "Fraction", "AttributeType")
  is_yes  <- nrow(freq_df[freq_df$Significant == "YES", ]) == 1
  if (is_yes) {
    freq_df <- freq_df[freq_df$Significant == "YES", ] 
  }else {
    freq_df <-  data.frame(Significant="YES",Fraction=0,
                           AttributeType=unique(freq_df$AttributeType))
  }
  return(freq_df$Fraction)
}

############################################################################
prepare_for_plot <- function(ol_summ,cfg)  { 
  df <- as_tibble(ol_summ)
  max_finite_OR <- max(df[["OR"]][is.finite(df[["OR"]])])
  df <- df %>%
    mutate(
      OR_plot = ifelse(is.infinite(OR), 2 * max_finite_OR,OR),
      sig     = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05")
    )
  
  df <- df %>%
    mutate(
      Attribute = factor(Attribute, levels = c("CELL_LINE","SUBTYPE","CANCER","TISSUE")),
      DB        = factor(DB, levels = names(cfg$db_names))
    )
  return(df)
}

############################################################################
build_outlier_summary <- function(score_only_df) {
  split_by_db <- split(score_only_df, score_only_df$Dataset)
  
  outlier_tables <- lapply(names(split_by_db), function(dataset_name) {
    input_df <- split_by_db[[dataset_name]]
    input_df <- input_df[input_df$attributeType != "HISTOLOGY", ]
    
    outlier_df <- do.call("rbind", lapply(c("attributeType", "sameLib", "Genes"), function(split_var) {
      print(split_var)
      df <- analyzeOutliers(input_df = input_df, split_var = split_var, variable_to_quant = "Correlation")
      df$Dataset <- dataset_name
      df$Variable_name <- gsub("\\.", "", rownames(df))
      df[df$Variable_name != "TRUE", ]
    }))
    
    outlier_df$Significant <- outlier_df$Fisher_p < 0.05
    outlier_df <- reshape(
      outlier_df[, c("Dataset", "Variable_name", "Significant")],
      direction = "wide",
      timevar = "Variable_name",
      idvar = "Dataset"
    )
    names(outlier_df) <- sub("^Significant\\.", "", names(outlier_df))
    rownames(outlier_df) <- NULL
    
 #   if (dataset_name == "depMap_full") {
#      outlier_df$CELL_LINE <- NA
#    }
    
    outlier_df
  })
  
  do.call("rbind", outlier_tables)
}

############################################################################
build_random_vs_real_library_df <- function(score_only_df, cfg,all_comb=F) {
  random_df <- do.call("rbind", lapply(names(cfg$db_names), function(db) {
    db_name <- cfg$db_names[db]
    print(db_name)
    if (db=="iCSDB") {
      df <- readRDS( paste0(db_name, "_fullListOfScoresRankOverlapsRandomFile.2.rds"))
      df <- df[["zscores"]] 
    } else {

    df <- readRDS( paste0(db_name, "_fullListOfScoresRankOverlapsRandomFile.rds"))
    df <- df[["zcores"]]
    }
    df$sameLib <- sapply(df$libUsed, addLibraryInfo)
    df$Dataset <- db
    df$Attribute <- "Different"
    df[["Library+Attribute"]] <- paste(df$sameLib, df$Attribute, sep = "+")
    df[, c("Correlation", "AttributeType", "Library+Attribute", "Dataset")]
  }))
  if (all_comb==F) {
    real_df <- score_only_df[score_only_df$sameLib == "Different", c("Correlation", "attributeType", "sameLib", "Dataset")]
  } else {
  real_df <- score_only_df[, c("Correlation", "attributeType", "sameLib", "Dataset")]
  }
  names(real_df) <- c("Correlation", "AttributeType", "sameLib", "Dataset")
  real_df$Attribute <- "Same"
  real_df[["Library+Attribute"]] <- paste(real_df$sameLib, real_df$Attribute, sep = "+")
  real_df <- real_df[, c("Correlation", "AttributeType", "Library+Attribute", "Dataset")]
  
  out <- rbind(real_df, random_df)
  out[out$AttributeType != "HISTOLOGY", ]
}


############################################################################
run_library_effect_analysis <- function(repro, inputs, cfg) {
  score_only_df <- repro$score_only
  outlier_summary <- build_outlier_summary(score_only_df)
  
  category_map <- data.frame(
    set = names(subset(outlier_summary, select = -Dataset)),
    Variable = c(rep("Biological attribute", 4), rep("gRNA library used", 2), rep("Gene subset", 4))
  )
  stripe_cols <- c("Biological attribute" = "steelblue", "gRNA library used" = "firebrick", "Gene subset" = "darkgreen")
  outlier_upset_plot <- ComplexUpset::upset(
    outlier_summary,
    names(subset(outlier_summary, select = -Dataset)),
    annotations = list(
      "Dataset" = ggplot(mapping = aes(fill = Dataset)) +
        geom_bar(stat = "count", position = "fill") +
        theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              legend.position = "right",
              legend.key.size = unit(0.4, "cm"), 
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 9))
    ),
    set_sizes = FALSE,
    base_annotations = list(),
    stripes = upset_stripes(mapping = aes(color = Variable), colors = stripe_cols, data = category_map),
    height_ratio = 3
  )
  save_plot_pdf(cfg, "Fig_3a.pdf", outlier_upset_plot,6,8) 
  
  random_vs_real <- build_random_vs_real_library_df(score_only_df, cfg,all_comb = T)
  random_vs_real <- random_vs_real[random_vs_real[["Library+Attribute"]]!="Different+Different",]
  random_vs_real_plot <- ggplot(random_vs_real,
                                aes(Dataset, Correlation, fill = .data[["Library+Attribute"]])) +
    geom_boxplot() +
    stat_compare_means(label = "p.signif", vjust = 0.05) +
    facet_wrap(vars(AttributeType)) +
    theme(axis.text.x = element_text(angle = 90),
          legend.key.size = unit(0.4, "cm"), 
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 13,hjust = 0.5),
          legend.position = "top",
          legend.title.position = "top"
          ) +xlab("") + ylim(c(-1.0,1.1))
  save_plot_pdf(cfg, "Fig_3b.pdf",random_vs_real_plot,8,8)
  
  invisible(check_ss_versus_ds("biogrid"))
  
  #### report "good" and "bad" libraries 
  library_summary <- merge(inputs$library_metadata, inputs$prediction_of_library_effect, by.y=0, by.x="LIBRARY_ID",all=F)
  library_summary <- library_summary[library_summary$NUM_gene > 10000,c("LIBRARY_NAME","AUTHOR","NUM_gRNA",
                                                                        "NUM_gene","NUM_gRNA_PER_GENE","Est","Std_err",
                                                                        "t_val","p_val")]
  invisible(two_step_lib_lolli(library_summary))
  library_summary$p_val <- scientific( library_summary$p_val, digits = 3) 
  kable(library_summary, format = "latex", booktabs = TRUE,
        caption = "Results from the two-step regresssion analysis of library occurrence vs correlation scores.",
        label = "Tab_1") %>%
    kable_styling(latex_options = "scale_down") %>%
    cat(file = file.path(cfg$plot_dir,"Table_1.tex"))
  
  ###############
  full_overlap_df <- prepare_overlap_summary(inputs$overlap_real, inputs$overlap_random)
  overlap_ratios <- do.call("rbind", lapply(names(cfg$db_names)[1:4], overlap_ratios_fun,
                                            full_overlap_df=full_overlap_df))
  overlap_lollipop_plot <- prepare_loolipop_plot(overlap_ratios)

  save_plot_pdf(cfg, "Fig_4a.pdf",overlap_lollipop_plot, 8,11)
  
  summary_of_overlap_fisher <- prepare_summary_of_overlap_fisher_fun(cfg = cfg,
                                                                    full_overlap_df=full_overlap_df)
  
  ol_summ_plot_df <- prepare_for_plot(summary_of_overlap_fisher$summary_of_overlap_fisher_la,
                                      cfg = cfg)
  ol_summ_plot <- ggplot(ol_summ_plot_df[ol_summ_plot_df$DB!="DepMap",], 
                         aes(x = DB, y = Attribute, size = OR_plot, fill = sig)) +
    geom_point(shape = 21, color = "black", stroke = 0.8) +
    scale_fill_manual(values = c("p < 0.05" = "red", "p >= 0.05" = "blue"),
                      name = "Significance") +
    scale_size_continuous(name = "Odds ratio", range = c(2, 13),) +
    labs(x = "Database", y = "Attribute") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust=1),
      axis.text.y = element_text(hjust = 1),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size=11,hjust = 0.5),
      legend.position = "bottom",
      legend.title.position = "top",
      legend.spacing.x = unit(3, "cm"),
      legend.background = element_rect(
        color = "black",       # Border line color
        linewidth = 0.5,       # Border line thickness
        fill = "grey98"        # Background color inside the box (use "transparent" if you want it clear)
      )
    )  +   
    guides(
      size = guide_legend(
        ncol = 2,       # Force exactly 2 columns
        byrow = TRUE    # Fill row-by-row (Group 1 & 2 on top, 3 & 4 on bottom)
      ),
      fill = guide_legend(
        ncol = 1,       
        byrow = TRUE,
        override.aes = list(
          size = 4,
          shape = 21)
      )
    ) + xlab("") + ylab("")
  save_plot_pdf(cfg, "Fig_4b.pdf",ol_summ_plot,8, 9)

  same_vd_diff_lib_no_plus_plot <- ggplot(score_only_df, aes(Dataset, Correlation, fill = sameLib)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90)) +
    stat_compare_means(label = "p.signif") + labs(fill="gRNA library used") +
    xlab("")
  save_plot_pdf(cfg, "Fig_S3.pdf",same_vd_diff_lib_no_plus_plot,6, 6)
  
  corr_vs_grnas <- correlation_vs_grnas(inputs,repro_results$score_only,
                                        inputs$prediction_of_library_effect)
  corr_vs_grnas_plot <- ggplot(corr_vs_grnas, aes(Metric,Value,fill=Correlation)) + geom_boxplot() +
    stat_compare_means(label="p.signif", vjust=-0.5)   + 
    facet_wrap(vars(Metric), scales="free") +
    theme(axis.text.x = element_blank()) + xlab("") + scale_fill_manual(values=c("firebrick", "steelblue")) +
    scale_y_continuous(labels = scales::label_comma())
  save_plot_pdf(cfg, "Fig_5b.pdf",corr_vs_grnas_plot,6, 6)
  
  
  list(summary_of_overlap_fisher_la=summary_of_overlap_fisher$summary_of_overlap_fisher_la,
              full_overlap_df = full_overlap_df)
}

