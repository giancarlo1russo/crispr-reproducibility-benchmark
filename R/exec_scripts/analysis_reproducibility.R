###################################
run_reproducibility_analysis <- function(inputs, cfg, score_vs_rank=FALSE) {
  
  if (!score_vs_rank) {
    corr_df <- combine_intra_db_correlations(inputs$intra_db_corr, cfg)
    cov_df <- prepare_cov_summary(cov_df = inputs$cov_summary,db_lookup=inputs$db_lookup)
    corr_all <- corr_df[corr_df$Genes == "All", ]
    corr_non_all <- corr_df[corr_df$Genes != "All", ]
    score_only <- corr_df[corr_df$Metric == "Score", ]
  
  corr_all_plot_2 <- ggplot(corr_all, aes(Dataset, Correlation, fill = Metric)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90)) +
    stat_compare_means(label = "p.signif") +
    facet_wrap2(vars(attributeType),scales="free_y",remove_labels = "x") + xlab("") 
  save_plot_pdf(cfg, "Fig_S1.pdf",corr_all_plot_2,8, 8)
  
  
  cov_df_plot_2 <-ggplot(cov_df, aes(Dataset, Value, fill = Metric)) +
    geom_boxplot() +
    stat_compare_means(label = "p.signif") +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap2(vars(attributeType),scales="free_y",remove_labels = "x") + 
  ylab("Coefficient of variation") +xlab("")
  save_plot_pdf(cfg, "Fig_S2.pdf",cov_df_plot_2,8, 8)

  
  corr_all_plot <- ggplot(corr_all, aes(Dataset, Correlation, fill = Metric)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90)) +
    stat_compare_means(label = "p.signif") + xlab("")
  save_plot_pdf(cfg, "Fig_2a.pdf",corr_all_plot,6, 6)
  
  cov_df_plot <- ggplot(cov_df, aes(Dataset, Value, fill = Metric)) +
    geom_boxplot() +
    stat_compare_means(label = "p.signif") +
    theme(axis.text.x = element_text(angle = 90)) + 
    ylab("Coefficient of variation")+xlab("")
  save_plot_pdf(cfg, "Fig_2b.pdf",cov_df_plot,6, 6)
  
  top_genes_plot <- ggplot(corr_non_all[corr_non_all$Genes=="Top",], 
                           aes(Dataset, Correlation, fill = Metric)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90)) +
    stat_compare_means(label = "p.signif") 
  save_plot_pdf(cfg, "Fig_2c.pdf",top_genes_plot,6, 6)
  
  
  compare_df <- rbind(
    cov_df[cov_df$Metric == "Score", c("Type", "Value", "Dataset", "attributeName", "attributeType")],
    transform(
      setNames(
        corr_all[corr_all$Metric == "Score", c("Correlation", "Dataset", "attributeName", "attributeType")],
        c("Value", "Dataset", "attributeName", "attributeType")
      ),
      Type = "Correlation"
    )
  )
  compare_ds_for_plot <- compare_df[compare_df$attributeType=="CELL_LINE",]
  dataset_comp_plot <- ggplot(compare_ds_for_plot, 
                              aes(Type, Value, fill = Dataset)) +
    geom_boxplot() +
    stat_compare_means(label = "p.signif") +  ylab("") +xlab("") +
    scale_fill_brewer(palette = "Set2")
  save_plot_pdf(cfg, "Fig_3d.pdf",dataset_comp_plot, 6,6)

  invisible(create_summary_plots(inputs, cfg))
  return(list(corr_df = corr_df, cov_df = cov_df, score_only = score_only))
  
} else {
  score_rank_df <- build_score_rank_summary(inputs$all_info, cfg)
  names(score_rank_df) <- names(inputs$all_info)
  bg_temp <- score_rank_df$BioGRID[,names(score_rank_df$DepMap)]
  score_rank_df_no_biogrid <- rbind(score_rank_df$DepMap,
                                    score_rank_df$iCSDB,bg_temp)
  score_rank_df_no_biogrid_plot <- ggplot(score_rank_df_no_biogrid, aes(RANK, BCSCORE)) +
    geom_point() +
    facet_wrap(vars(Dataset), scales = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
    scale_x_continuous(labels = scales::label_comma())
  save_plot_pdf(cfg, "Fig_1a.pdf",score_rank_df_no_biogrid_plot,12, 8) 
  
  score_rank_df_plot <- ggplot(score_rank_df$BioGRID, aes(RANK, BCSCORE)) +
    geom_point() +
    facet_wrap(vars(SCORE_TYPE), scales = "free",nrow=2) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_x_continuous(labels = scales::label_comma())
  save_plot_pdf(cfg, "Fig_1b.pdf",score_rank_df_plot,8, 8)

}
}

############# create summary plots
create_summary_plots <- function(inputs, cfg) {
long_input <- prepare_md_long(metadata_file_names = names(inputs$metadata_list), 
                              metadata_files = inputs$metadata_list,
                              attrs = cfg$biological_attributes)
summary_attr_plot_list <- summarize_attributes(long_input)
summary_attr_plot  <- (print(summary_attr_plot_list$summ_nr_attr_plot) | print(summary_attr_plot_list$top5_attr_plot))+
  plot_annotation(tag_levels = "A")
save_plot_pdf(cfg, "Fig_S10.pdf",summary_attr_plot,16, 8)
}