# setup: load dependencies and install if necessary using pacman
#===============================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load("here",
               "dplyr", 
               "sjPlot",
               "psych",
               "tibble",
               "ggplot2",
               "GGally",
               "car", 
               "lmtest",
               "knitr", 
               "openxlsx", 
               "ggrain", 
               "ggdist",
               "elucidate", 
               "tidyr", 
               "corrplot", 
               "afex", 
               "ggpubr", 
               "see", 
               "reshape2", 
               "stats", 
               "multcomp", 
               "xtable", 
               "cocor",
               "kableExtra", 
               "modelsummary", 
               "broom", 
               "purrr",
               "conflicted",
               "patchwork",
               "magick",
               "ggplotify",
               "stringr",
               "performance")


# Avoid namespace conflicts 
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("ggsave", "ggplot2")
conflict_prefer("area", "patchwork")

