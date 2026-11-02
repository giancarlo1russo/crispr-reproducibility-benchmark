setwd("/Users/giancarlorusso/mainFolder_GMC_VU/Gediminas/crisprScreenDatabases/new_scripts/")
source("config.R")
source("listOfFunctions.R")
source("functionsForPlots.R")

source("helpers_io.R")
source("helpers_plotting.R")
source("analysis_reproducibility.R")
source("analysis_library_effects.R")
source("analysis_biological_databases.R")

setwd(cfg$root_dir)

inputs <- load_analysis_inputs(cfg)
names(inputs$all_info) <- sapply(names(inputs$all_info), 
                                 function(x) names(cfg$db_names)[which(x==cfg$db_names)])
names(inputs$intra_db_corr) <- names(inputs$all_info) 
invisible(run_reproducibility_analysis(inputs, cfg,score_vs_rank = T))

#### remove full BioGRID
cfg$db_names <- cfg$db_names[-match("BioGRID", names(cfg$db_names))]
inputs <- reformat_inputs(inputs)
####

#### analyse reproducibility 
repro_results <- run_reproducibility_analysis(inputs, cfg)

#### analyse library effect
library_results <- run_library_effect_analysis(repro_results, inputs, cfg)

#### integrate with other biological datadases 
cancer_corr_results <- run_cancer_analysis(repro_results, inputs, cfg)

mutational_burden_plot <- mutation_db_burden_analysis(repro_results,inputs,cfg,
                                                 dbName="DepMap")
dentro_data_plot <- hetero_dentro(inputs$heterogeneity,cancer_corr_results$Cancer)
methy_plot_df <- lapply(names(inputs$methyl_files),methyl_burden_prep,
                        meth_list=inputs$methyl_files,inputs=inputs,cfg=cfg)
mehtyl_plot <- methyl_plot_fun(methy_plot_df)
exon_dependance_plot<- exon_dependance_plot_fun(input_df = inputs$intra_db_corr$DepMap$zcores,
                               iso_data = inputs$iso_data)


