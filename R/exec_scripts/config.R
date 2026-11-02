suppressPackageStartupMessages({
  library(ComplexUpset)
  library(GenomicRanges)
  library(cowplot)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(scales)
  library(stringr)
  library(tidyr)
  library(tools)
  library(plyr)
  library(patchwork)
  library(ggrepel)
  library(knitr)
  library(kableExtra)
  library(ggh4x)
  library(tidyverse)
  library(tidytext)
  library(ggbreak)
})

cfg <- list(
  root_dir= "/Users/giancarlorusso/mainFolder_GMC_VU/Gediminas/crisprScreenDatabases/new_data_and_scripts/data",
  plot_dir = "/Users/giancarlorusso/mainFolder_GMC_VU/Gediminas/crisprScreenDatabases/finalFigsAndTables_new",
  biological_db_dir = "/Users/giancarlorusso/mainFolder_GMC_VU/Gediminas/crisprScreenDatabases/new_data_and_scripts/data/biological_databases",
  db_names = c(DepMap="depMap_full", BioGRID_exp="biogrid_exp",
               BioGRID_s="biogrid_s",iCSDB="icsdb",BioGRID="biogrid"),
  literature_db_names = c("BioGRID_s","BioGRID_exp", "iCSDB"),
  metrics = c("zcores", "rank"),
  biological_attributes = c("CELL_LINE", "SUBTYPE", "CANCER", "TISSUE")
)
