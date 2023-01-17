---
title: "Install Packages"
author: "S.Y. CHIA"
date: "January 16 2023"
output: html_document
editor_options: 
chunk_output_type: console
---

# install the packages that will be needed to install other packages
  renv::restore(packages = c("devtools",
                             "BiocManager",
                             "remotes"))


# install some of the basic packages
renv::restore(packages = c( "ggplot2",
                            "dplyr",
                            "tibble",
                            "tidyr", 
                            "car", 
                            "agricolae", 
                            "ggpubr", 
                            "emmeans", 
                            "multcompView", 
                            "plotrix", 
                            "lme4", 
                            "tidyverse", 
                            "RColorBrewer", 
                            "colortools", 
                            "GGally", 
                            "lattice", 
                            "ggthemes", 
                            "rmarkdown",
                            "DHARMa",
                            "patchwork", 
                            "moments", 
                            "glmmTMB"
                            ))


# install all other packages
renv::restore()


# update all packages
renv::update()


#save renv state
renv::snapshot()


#check status
renv::status()

  
  
