### Input files and dirs
root_dir <- "/Users/giancarlorusso/mainFolder_GMC_VU/Gediminas/crisprScreenDatabases/new_data_and_scripts"
data_dir <- file.path(root_dir,"data")
script_dir <- file.path(root_dir,"scripts")

### Source needed files 
source(file.path(script_dir,"data_prep_inputs2_generic_functions.R"))
source(file.path(script_dir,"load_libraries.R"))
source(file.path(script_dir,"listOfFunctions.R"))

### Define/read in generic inputs
attribute_list <- c(CELL_LINE="CELL_LINE",SUBTYPE="SUBTYPE",
                    CANCER="CANCER",TISSUE="TISSUE")

### Run main functions 
database_names <- c("depMap_full","biogrid_exp","biogrid_s","icsdb","biogrid")
modes <- c(zcores="zcores",rank="rank")
invisible(generate_all_info_list(dbNames = database_names,data_dir = data_dir))
invisible(generate_intraDB_list(dbNames = database_names,data_dir = data_dir,
metrics=modes))
sd_db_list <- invisible(lapply(database_names,generate_all_db_sd_df_db_list,
                                  data_dir = data_dir,attribute_list))
sd_db_files <- file.path(data_dir, paste0(database_names,"_sd_df.rds"))
sd_list_file <- lapply(sd_db_files, readRDS)
names(sd_list_file) <- database_names
all_db_sd_df_db_list <- do.call("rbind",sd_list_file)
saveRDS(all_db_sd_df_db_list,"all_db_sd_df_db_list.rds")
invisible(create_library_corr_coeff_matrix(intraDB_list_file=intraDB_list_file, dbNames=database_names,
                                         db_specific_name = "icsdb",
                                          metrics=modes, colnameLibs = "libUsed",
                                          colNameCorrCoeff = "Correlation"))
