## Benchmarking CRISPR Screen Reproducibility Across Literature Repositories and DepMap

This repository contains the R workflows and analysis pipeline for:
Manuscript: A Systematic Multi-Repository Assessment of Reproducibility in Cancer CRISPR Knockout 
Screens
Authors: G.Gumbis, G.Russo
Journal: [IMJS / Accepted]
Corresponding author: [giancarlo.russo@gmc.vu.lt]

### Structure of the repository
.
├── .gitignore
├── LICENSE.txt
├── README.md
├── data/
│   └── README.md
├── figures/
│   ├── Fig_1a.pdf
│   ├── Fig_1b.pdf
│   ├── Fig_2a.pdf
│   ├── Fig_2b.pdf
│   ├── Fig_2c.pdf
│   ├── Fig_3a.pdf
│   ├── Fig_3b.pdf
│   ├── Fig_3c.pdf
│   ├── Fig_3d.pdf
│   ├── Fig_4a.pdf
│   ├── Fig_4b.pdf
│   ├── Fig_5a.pdf
│   ├── Fig_5b.pdf
│   ├── Fig_6a.pdf
│   ├── Fig_6b.pdf
│   ├── Fig_S1.pdf
│   ├── Fig_S2.pdf
│   ├── Fig_S3.pdf
│   ├── Fig_S4.pdf
│   ├── Fig_S5.pdf
│   ├── Fig_S6.pdf
│   ├── Fig_S7.pdf
│   ├── Fig_S8.pdf
│   ├── Fig_S9.pdf
│   └── Fig_S10.pdf
└── R/
    ├── exec_scripts/
    │   ├── analysis_biological_databases.R
    │   ├── analysis_library_effects.R
    │   ├── analysis_reproducibility.R
    │   ├── config.R
    │   ├── functionsForPlots.R
    │   ├── helpers_io.R
    │   ├── helpers_plotting.R
    │   ├── listOfFunctions.R
    │   └── run_all.R
    └── prep_scripts/
        ├── data_prep_all_info_execution.R
        ├── data_prep_functions.R
        ├── data_prep_inputs1_execution.R
        ├── data_prep_inputs1_generic_functions.R
        ├── data_prep_inputs2_execution.R
        ├── data_prep_inputs2_generic_functions.R
        ├── listOfFunctions.R
        └── load_libraries.R

### Data Sources & Acquisition
This pipeline operates directly from primary raw datasets hosted by the 
respective portals. No manual preprocessing is required. Details on the sources of  data in 
data/README.md.

### R version: 4.6.1, tested on: macOS

### Required R libraries: 

library(base)
library(ComplexUpset)
library(corrplot)
library(cowplot)
library(data.table)
library(dplyr)
library(factoextra)
library(GenomicRanges)
library(ggbreak)
library(ggh4x)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gplots)
library(kableExtra)
library(knitr)
library(patchwork)
library(plyr)
library(RColorBrewer)
library(readr)
library(RImagePalette)
library(rstatix)
library(scales)
library(stats)
library(stringr)
library(tidytext)
library(tidyr)
library(tidyverse)
library(tools)
library(UpSetR)
library(utils)
library(vegan)


### Reproducing the Analysis & Figures

#### The preparatory scripts are run in the following order: 

1. R/prep_scripts/data_all_info_execution.R	
2. R/prep_scripts/data_prep_inputs1_execution.R
3. R/prep_scripts/data_prep_inputs2_execution.R

#### The main script is: 
1. R/exec_scripts/run_all.R 

### All figures are in folder ./figures
