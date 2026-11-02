### Input files and dirs
root_dir <- "/Users/giancarlorusso/mainFolder_GMC_VU/Gediminas/crisprScreenDatabases/new_data_and_scripts"
data_dir <- file.path(root_dir,"data")
script_dir <- file.path(root_dir,"scripts")

### Source needed files 
source(file.path(script_dir,"data_prep_inputs1_generic_functions.R"))
source(file.path(script_dir,"load_libraries.R"))
source(file.path(script_dir,"listOfFunctions.R"))

### Define/read in generic inputs
attribute_list <- c(CELL_LINE="CELL_LINE",SUBTYPE="SUBTYPE",
                    CANCER="CANCER",TISSUE="TISSUE")
perc_levels <- c(point1=0.1,point5=0.5,one=1,five=5)
modes <- c(zcores="zcores",rank="rank")
directions <- c(All="All",Top="Top", Bottom="Bottom", Middle="Middle")
database_matrix <- file.path(data_dir,"database_matrix.tsv")
top_bottom_genes <- read_top_bottom_genes(database_matrix)
attribute_levels <-  read.delim(file.path(data_dir,"attribute_level_matrix.tsv"))

### Run main functions 
database_names <- c("depMap_full","biogrid_exp","biogrid_s","icsdb")
for (databaseName in database_names) {
  print(databaseName)
  allInfoFilePath <- file.path(data_dir,paste0(databaseName,"_compact_with_all_info_noNA.rds"))
  md_file_path <- file.path(data_dir,paste0(databaseName,"_metadata.rds"))
  prep_data <- generate_overlaps_prep(allInfoFilePath,attribute_list,attribute_levels)
  invisible(generate_gene_overlaps_real(allInfoFileSplitList = prep_data$all_list,
                                                attributeList = attribute_list,
                                               attributesWithAtLeastTwoScreens=prep_data$with_two,
                                               percLevels = perc_levels,data_dir = data_dir,
                                         database_name = databaseName, md_file_path=md_file_path
                                         ))
 invisible(generate_gene_overlaps_random(allInfoFileSplitList = prep_data$all_list,
                                              attributeList = attribute_list,
                                              percLevels = perc_levels,data_dir = data_dir,
                                          md_file_path=md_file_path,database_name=databaseName))
  invisible(generate_score_rank_corr_real(allInfoFileSplitList= prep_data$all_list,
                                                attributeList = attribute_list,
                                                modes = modes,data_dir = data_dir,
                                                attributesWithAtLeastTwoScreens=prep_data$with_two,
                                                directions = directions,database_name=databaseName,
                                                NgenesTop = top_bottom_genes$top[[databaseName]],
                                                NgenesBottom = top_bottom_genes$bottom[[databaseName]]))
   invisible(generate_score_rank_corr_random(allInfoFileSplitList= prep_data$all_list,
                                                   attributeList = attribute_list,
                                                   modes = modes,database_name=databaseName,
                                                   directions = directions,data_dir = data_dir,
                                                   NgenesTop = top_bottom_genes$top[[databaseName]],
                                                   NgenesBottom = top_bottom_genes$bottom[[databaseName]],
                                             attribute_levels = attribute_levels, md_file_path=md_file_path))
}
