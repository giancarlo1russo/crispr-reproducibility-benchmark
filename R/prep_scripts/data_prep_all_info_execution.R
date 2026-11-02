root_dir <- "/Users/giancarlorusso/mainFolder_GMC_VU/Gediminas/crisprScreenDatabases/new_data_and_scripts"
source(file.path(root_dir,"scripts/data_prep_functions.R"))
source(file.path(root_dir,"scripts/load_libraries.R"))
source(file.path(root_dir,"scripts/listOfFunctions.R"))

### input files and dirs
depMap_dir <- file.path(root_dir,"data/depMap")
biogrid_dir <- file.path(root_dir,"data/biogrid")
data_dir <- file.path(root_dir,"data")

### prepare metadata 
invisible(prepare_depMap_metadata(depMap_dir,data_dir))
depMap_metadata <- file.path(data_dir,"depMap_full_metadata.rds")
depMap_screens_to_use <- readRDS(depMap_metadata)[["EWHA_ID"]]

### prepare data: all info
biogrid_score_subset <- list(biogrid_exp=c("Bayes Factor","STARS Score"),
                             biogrid_s=c("CERES score","CRISPR Score (CS)","Z-score"))
cell_line_map <- prepare_cell_map_table(depMap_dir)
invisible(prepare_biogrid_data(biogrid_dir, cell_line_map,
                               data_dir,biogrid_score_subset))
invisible(prepare_depMap_data(depMap_dir,data_dir,depMap_screens_to_use))


