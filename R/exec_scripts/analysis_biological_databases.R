######################################################
mut_burd_aggr_list_fun <- function(attr,depMap_only,mut_burd_with_MD,var_type) {
  print(attr)
  formula_freq <- as.formula(paste0("Freq ~",attr))
  depMap_only_sub <- depMap_only[depMap_only$attributeType==attr,]
  dm_aggr_corr <-  aggregate(Correlation ~ attributeName,data=depMap_only_sub, FUN=mean)
  if (attr=="TISSUE"){
    dm_aggr_corr[["attributeName"]] <- toTitleCase(dm_aggr_corr[["attributeName"]])
    dm_aggr_corr[["attributeName"]] <- gsub("_"," ",dm_aggr_corr[["attributeName"]])
  }
  names(dm_aggr_corr)[1] <- attr
  dm_aggr_freq <- aggregate(formula_freq,data=mut_burd_with_MD, FUN=sum)
  dm_burd_corr_merged <- merge(dm_aggr_freq,dm_aggr_corr, by=attr, all=F)
  dm_burd_corr_merged$Variant_type <- var_type
  names(dm_burd_corr_merged)[1] <- "AttributeName"
  dm_burd_corr_merged$AttributeType <- attr
  return(dm_burd_corr_merged)
}

######################################################
temp_mut_burd_list_fun <- function(var_type,mutational_burden,cfg,depMap_only,
                                   depMap_metadata_2) {
  print(var_type)
  temp_mut_burd <- data.frame(table(mutational_burden[mutational_burden$VariantType==var_type,"ModelID"]))
  names(temp_mut_burd) <- c("EWHA_ID","Freq")
  mut_burd_with_MD <- merge(temp_mut_burd, depMap_metadata_2, by="EWHA_ID", all=F)
  mut_burd_aggr_list <- lapply(cfg$biological_attributes, mut_burd_aggr_list_fun,
                               depMap_only=depMap_only,mut_burd_with_MD=mut_burd_with_MD,
                               var_type=var_type)
  mut_burd_aggr <- do.call("rbind", mut_burd_aggr_list)
  return(mut_burd_aggr)
}


######################################################
cor_fun <- function(df,V1,V2,is_log=F) {
  if (is_log) {
  c1 <-  cor(log(df[[V1]]),df[[V2]])
  }else {
  c1 <-  cor(df[[V1]],df[[V2]])
  }
  c2 <- paste0("R2 = ",sprintf("%.2f", c1))
  return(c2)
}
######################################################
annotate_R2 <- function(df,V1,V2,x_pos,y_pos,is_log=F,var=var) {
to_plot_df_cor <- ddply(.data=df,as.quoted(var), cor_fun,
                        V1=V1,V2=V2,is_log=F)
to_plot_df_cor$x_pos <- x_pos
to_plot_df_cor$y_pos <- y_pos
return(to_plot_df_cor)
}
######################################################
standardize_cancer_names <- function(df) {
  df$attributeName <- gsub("_", " ", df$attributeName)
  df <- na.omit(df)
  df[grep("Renal", df$attributeName), "attributeName"] <- "Renal Clear Cell Carcinoma"
  df[grep("Lung Cancer", df$attributeName), "attributeName"] <- "Lung Adenocarcinoma"
  df[grep("astric", df$attributeName), "attributeName"] <- "Stomach Adenocarcinoma"
  df[grep("Leukemia",df$attributeName),]$attributeName <- "Leukemia"
  df[grep("Lymphoma",df$attributeName),]$attributeName <- "Lymphoma"
  df
}

######################################################
get_tissue_for_cancer <- function(cancer_name, all_info, cfg){ 
  for (db_name in names(cfg$db_names)){
    pp <- all_info[[db_name]]
    y <- unique(unlist(pp[pp[["CANCER"]]==cancer_name,"TISSUE"]))
    if (length(y)>0) {
      if (!is.na(y[1])) {
        y2 <- toupper(y[1])
        y2 <- gsub("_"," ",y2)
        break   
      }
    } else {
      y2 <- NA
    }
  }
  return(y2)
}

######################################################
run_cancer_analysis <- function(repro, inputs, cfg) {
    setwd(cfg$biological_db_dir)
    score_only_df <- repro$score_only
    
    cancer_df <- score_only_df[
      score_only_df$attributeType == "CANCER" & score_only_df$sameLib == "Different" & 
        score_only_df$Genes=="All",
    ]
    cancer_df <- standardize_cancer_names(cancer_df)
    cancer_df <- aggregate(Correlation ~ attributeName, data = cancer_df, FUN = mean)
    cancer_df$attributeName <- reorder(factor(cancer_df$attributeName), cancer_df$Correlation)
    cancer_df_for_plot <- na.omit(cancer_df)
    cancer_df_for_plot$Type <- "Cancer"
    kk <- na.omit(cancer_df)
    cancer_df_plot_supp <- ggplot(kk, aes(attributeName, Correlation)) +
      geom_bar(stat="identity") +
      coord_flip() +
      scale_x_discrete(labels = label_wrap(60)) +
      xlab("Mean correlation") + theme(axis.text.y = element_text(size=6))
    save_plot_pdf(cfg, "Fig_S6.pdf",cancer_df_plot_supp,12, 8)
   
    kk$Type <- "Cancer"
    cancer_df_plot_main <- ggplot(kk, aes(Type, Correlation)) +
      geom_boxplot() + xlab("")  +  geom_text_repel(
        data = filter(kk, Correlation <0.4), aes(label=attributeName),
        size = 3,hjust = 0,nudge_x = 0.2,
        box.padding = 0.5,        # Adds breathing room around the text
        point.padding = 0.3,      # Prevents the text from covering the dot
        segment.color = "grey50"  # Draws a clean line from text to dot
      )  + coord_flip()
    save_plot_pdf(cfg, "Fig_6a.pdf",cancer_df_plot_main,9, 7)
    
    
    
    
  cancer_df$TISSUE <- sapply(cancer_df$attributeName, get_tissue_for_cancer, 
                             all_info = inputs$all_info, cfg = cfg)
  tissue_df <- aggregate(Correlation ~ TISSUE, data = cancer_df, FUN = mean)
  tissue_df$TISSUE <- reorder(factor(tissue_df$TISSUE ), tissue_df$Correlation)
  tissue_corr_plot <- ggplot(tissue_df, aes(TISSUE, Correlation)) +
    geom_boxplot() +
    coord_flip() +
    scale_x_discrete(labels = label_wrap(60))
  save_plot_pdf(cfg, "not_used.pdf",tissue_corr_plot,6, 6)
  cancer_df <- cancer_df[,c("attributeName","Correlation")]
  names(cancer_df)[1] <- "Name"
  cancer_df$Attribute <- "Cancer"
  names(tissue_df)[1] <- "Name"
  tissue_df$Attribute <- "Tissue"
  
  return(list(Cancer=cancer_df,Tissue=tissue_df,Fig_6a=cancer_df_plot_main))
  # Keep the remaining burden / methylation / rearrangement / isoform blocks here.
  # The important refactor is that they now live in one isolated module rather than in the main script.
}

######################################################
mutation_db_burden_analysis <- function(repro,inputs,cfg,dbName="DepMap") {
  depMap_only <- repro$score_only
  depMap_only <- depMap_only[depMap_only$Dataset == dbName,]
  mutational_burden <- inputs$mutational_burden
  mutational_burden <- mutational_burden[mutational_burden$VariantType!= "substitution",]
  z <-  c(VariantType="VariantType")
  depMap_metadata_2 <- inputs$depMap_metadata
  depMap_metadata_2$EWHA_ID <- sapply(depMap_metadata_2$EWHA_ID, function(x) unlist(strsplit(x,"\\."))[1])
  depMap_metadata_2 <- depMap_metadata_2[!duplicated(depMap_metadata_2$EWHA_ID),]
  mutation_types <- names(table(mutational_burden$VariantType))
  temp_mut_burd_list <- lapply(mutation_types, temp_mut_burd_list_fun,
                               mutational_burden=mutational_burden,cfg=cfg,
                               depMap_only=depMap_only,depMap_metadata_2=depMap_metadata_2)
  temp_mut_burd <- do.call("rbind",temp_mut_burd_list)
  temp_mut_burd$AttributeType <- factor(temp_mut_burd$AttributeType, 
                                        levels=cfg$biological_attributes)
  to_plot_df <- temp_mut_burd[temp_mut_burd$AttributeType=="CANCER",]
  to_plot_df$Variant_type <- str_to_title(to_plot_df$Variant_type)
  to_plot_df_cor <- annotate_R2(to_plot_df,V1="Freq",V2="Correlation",
                                x_pos = rep(0.5,n=3),y_pos = c(10,10,11.5),
                                is_log=T,var="Variant_type")
  mut_burd_plot <- ggplot(to_plot_df, aes(Correlation,log(Freq))) + geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(data=to_plot_df_cor, aes(x = x_pos, y = y_pos,
                                       label=V1), size=5) +
    facet_wrap(vars(Variant_type), scales="free")  +
    ylab("Mutational burden")
    save_plot_pdf(cfg, "Fig_6b.pdf",mut_burd_plot,9,7)
    return(list(Fig_6b=mut_burd_plot))
}

######################################################
#### MUTATIONS: heterogeneity from dentro et al.

hetero_dentro <- function(het_df,cancer_df) {
het_df[["Percentage of subclones"]] <- apply(het_df[,grep("SC",names(het_df))],1,mean)
names(het_df)[names(het_df)=="SNVs"] <- "Percentage of mutated cells"
het_merged_df <- merge(cancer_df,het_df,by.x="Name",by.y=0,all=F)
het_merged_df <- het_merged_df[,c("Name","Correlation",
                                  "Percentage of subclones","Percentage of mutated cells")]
names(het_merged_df)[1] <- "Cancer"
het_merged_df_pl <- pivot_longer(het_merged_df,
                                 cols=c("Percentage of subclones","Percentage of mutated cells"),
                                 names_to="Type",values_to="Percentage")
het_merged_df_pl_cor <-annotate_R2(het_merged_df_pl,V1="Correlation",
                                   V2="Percentage", x_pos = rep(0.65,n=2),
                                   y_pos = c(0.47,0.383),var="Type")
het_merged_df_plot <- ggplot(het_merged_df_pl, 
       aes(Correlation,Percentage)) + geom_point() + 
  facet_wrap(vars(Type), scales="free") +  
  geom_smooth(method = "lm", se = FALSE)+
  geom_text(data=het_merged_df_pl_cor, aes(x = x_pos, y = y_pos,
                                           label=V1), size=5) +
  ylab("")

save_plot_pdf(cfg, "Fig_S7.pdf",het_merged_df_plot ,8,6)
return(list(SuppFig_4b=het_merged_df_plot))
}


######################## methyl burden
methyl_burden_prep <- function(n,meth_list,inputs,cfg) {
corr_df <-  combine_intra_db_correlations(inputs$intra_db_corr, cfg)
corr_df <- corr_df[corr_df$Genes=="All",]
methyl_ds <- meth_list[[n]]
methyl_ds_cl <- as.data.frame(lapply(methyl_ds[,4:ncol(methyl_ds)], as.numeric))
weights <- as.numeric(methyl_ds$avg_coverage)

### intersect with cell lines 
cell_line_names <- sapply(names(methyl_ds_cl), function(x) unlist(strsplit(x,"_"))[1])
names(methyl_ds) <- cell_line_names
methyl_df_burden <- data.frame(Cell_line=cell_line_names, 
                               Meth_level=apply(methyl_ds_cl, 
                                                2, function(x) weighted.mean(x,weights,na.rm=T)))
database <- "DepMap"
corr_depmap <- corr_df[corr_df$Dataset ==  database &
                         corr_df$attributeType=="CELL_LINE",
                       c("attributeName","Correlation")]
merged_df <- na.omit(merge(corr_depmap,methyl_df_burden,
                           by.y="Cell_line",by.x="attributeName",all=F))
merged_df$Meth_data <- n
return(merged_df)
}

########################
methyl_plot_fun <- function(input_df_list) {
  df <- do.call("rbind",input_df_list) 
  cor_data <- annotate_R2(df = df, V1="Correlation",
                        V2="Meth_level", x_pos = rep(0.55,3),
                        y_pos = c(0.55,0.45,0.35),var="Meth_data")

m_plot <- ggplot(df, aes(Correlation,Meth_level)) + geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(vars(Meth_data),scales="free") +
  geom_text(data=cor_data, aes(x = x_pos, y = y_pos,
                                           label=V1), size=5) +
  ylab("Methylation level")
  
save_plot_pdf(cfg, "Fig_S8.pdf",m_plot,10,6)
return(list(SuppFig_4c=m_plot))

}


#########################################################
### Exon dependence
######################################################### 
exon_dependance_plot_fun <- function(iso_data,input_df) {
input_df <- input_df[input_df$attributeType=="CELL_LINE" & 
                       input_df$Genes == "All",]
iso_merged <- merge(iso_data, input_df, by.x="CELL_LINE", by.y="attributeName",
                    all=F)
iso_merged_aggr <- aggregate(Correlation ~ CELL_LINE + Av_expressed_isof +Av_disp_expr,
                               data=iso_merged, FUN=mean)
names(iso_merged_aggr)[2:3] <- c("Mean isoform expression",
                               "Mean isoform dispersion")
iso_merged_aggr_pl <- pivot_longer(iso_merged_aggr,
                                   cols=c("Mean isoform expression","Mean isoform dispersion"),
                                 names_to="Type",values_to="Value")

cor_data <- annotate_R2(df = iso_merged_aggr_pl, V1="Value",
                        V2="Correlation", x_pos = c(1.53,0.7),
                        y_pos = rep(c(0.92),2),var="Type")

exon_plot <- ggplot(iso_merged_aggr_pl, aes(Value,Correlation)) + geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(vars(Type),scales="free") +
  geom_text(data=cor_data, aes(x = x_pos, y = y_pos,
                               label=V1), size=5) 
save_plot_pdf(cfg, "Fig_S9.pdf",exon_plot, 8,6)
return(list(SuppFig_4d=exon_plot))
}
