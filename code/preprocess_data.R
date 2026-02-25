#===============================================================================
# DATA PREPROCESSING AND LOADING
#===============================================================================
# This script loads and preprocesses data
# Author: Suhas Vijayakumar

# Useful abbreviations: 
#===============================================================================
# AF    = Arcuate Fasciculus 
# SLF   = Superior Longitudinal Fasciculus (SLF II/III)
# WM    = white matter
# GM    = gray matter
# L     = left hemipshere
# R     = right hemipshere
# AQ    = Asymmetry Quotient 
# PCA   = Principal Component Analysis (PC1, PC2 - first two components of PCA)
#===============================================================================
library(tidyverse)
library(here)

here::i_am("code/preprocess_data.R")

#===============================================================================
# function: to calculate asymmetry quotient (AQ)
#===============================================================================
asymmetry_quotient <- function(array_R_hemi, array_L_hemi) {
  # Convert inputs to numeric (in case they are factors/characters)
  R <- as.numeric(array_R_hemi)
  L <- as.numeric(array_L_hemi)
  
  # Case 1: both NA → result should be NA
  both_na <- is.na(R) & is.na(L)
  
  # Replace NA with 0
  R[is.na(R)] <- 0
  L[is.na(L)] <- 0
  
  # Compute AQ
  AQ <- 2 * (R - L) / (R + L)
  
  # Handle division by zero explicitly (R + L == 0)
  AQ[R + L == 0] <- NA
  
  # Restore both-NA cases
  AQ[both_na] <- NA
  
  return(AQ)
}

# Load preprocessed data
#===============================================================================
if (!file.exists(here("data", "preprocessed_data.csv"))) {
  cat("Preprocessed data not found — running preprocessing...\n")
  
  # load raw data
  #=============================================================================
  raw_data  <- read.csv(here("data","raw_data.csv"), header = TRUE, stringsAsFactors = FALSE)
  data      <- raw_data

  # PCA: flip PC1 component
  #=============================================================================
  data$PC1      <- -1*data$PC1 

  # syntactic simplicity: flip scores
  #===============================================================================
  data$syntactic_simplicity <- -1*data$syntactic_simplicity

  # Rename columns to avoid confusion
  #===============================================================================
  data <- data %>% 
    dplyr::rename(
      PC_thin = PC1,
      PC_long = PC2,
      syntactic_complexity = syntactic_simplicity
    )

  # calculate asymmetry quotients (AQ) and store in variable "data"
  # AQ < 0: left-lateralized, AQ > 0: right-lateralized
  #===============================================================================
  # for AF and SLF II/III based on average FA
  data$AQ_AF_averageFA  <- asymmetry_quotient(data$AF_R_averageFA, data$AF_L_averageFA)
  data$AQ_SLF_averageFA <- asymmetry_quotient(data$SLF_R_averageFA, data$SLF_L_averageFA)

  # for AF and SLF II/III based on WM tract volume
  data$AQ_AF_volumeWM   <- asymmetry_quotient(data$AF_R_volumeWM, data$AF_L_volumeWM)
  data$AQ_SLF_volumeWM  <- asymmetry_quotient(data$SLF_R_volumeWM, data$SLF_L_volumeWM)

  # for AF GM termination volume
  data$AQ_AF_PMv        <- asymmetry_quotient(data$AF_R_PMv, data$AF_L_PMv)
  data$AQ_AF_BA44       <- asymmetry_quotient(data$AF_R_BA44, data$AF_L_BA44)
  data$AQ_AF_BA45       <- asymmetry_quotient(data$AF_R_BA45, data$AF_L_BA45)
  data$AQ_AF_STG        <- asymmetry_quotient(data$AF_R_STG, data$AF_L_STG)
  data$AQ_AF_MTG        <- asymmetry_quotient(data$AF_R_MTG, data$AF_L_MTG)
  data$AQ_AF_ITG        <- asymmetry_quotient(data$AF_R_ITG, data$AF_L_ITG)
  data$AQ_AF_TP         <- asymmetry_quotient(data$AF_R_TP, data$AF_L_TP)

  # for SLF II/III GM termination volume
  data$AQ_SLF_PMv       <- asymmetry_quotient(data$SLF_R_PMv, data$SLF_L_PMv)
  data$AQ_SLF_BA44      <- asymmetry_quotient(data$SLF_R_BA44, data$SLF_L_BA44)
  data$AQ_SLF_BA45      <- asymmetry_quotient(data$SLF_R_BA45, data$SLF_L_BA45)
  data$AQ_SLF_SMG       <- asymmetry_quotient(data$SLF_R_SMG, data$SLF_L_SMG)
  data$AQ_SLF_AG        <- asymmetry_quotient(data$SLF_R_AG, data$SLF_L_AG)
  data$AQ_SLF_SPL       <- asymmetry_quotient(data$SLF_R_SPL, data$SLF_L_SPL)

  # Convert factors
  #===============================================================================
  data <- data %>%
    mutate(
      subject  = as.factor(subject),
      group    = factor(group, levels = c("con", "exp")),
      training = factor(training, levels = c("pre", "post"))
    )
  
  # Define variable groups for easy reference
  #===============================================================================
  var <- list()
  
  # white matter tract measures
  var$AQ_WM_measures          <- c("AQ_AF_averageFA", "AQ_AF_volumeWM", "AQ_SLF_averageFA", "AQ_SLF_volumeWM")
  var$WM_measures             <- c("AF_L_averageFA", "AF_R_averageFA", "AF_L_volumeWM", "AF_R_volumeWM", "SLF_L_averageFA", "SLF_R_averageFA", "SLF_L_volumeWM", "SLF_R_volumeWM")
  
  # AF: gray matter termination measures
  var$AQ_AF_GM_terminations   <- c("AQ_AF_PMv", "AQ_AF_BA44", "AQ_AF_BA45", "AQ_AF_STG", "AQ_AF_MTG", "AQ_AF_ITG", "AQ_AF_TP")
  var$AF_L_GM_terminations    <- c("AF_L_PMv", "AF_L_BA44", "AF_L_BA45", "AF_L_STG", "AF_L_MTG", "AF_L_ITG", "AF_L_TP")
  var$AF_R_GM_terminations    <- c("AF_R_PMv", "AF_R_BA44", "AF_R_BA45", "AF_R_STG", "AF_R_MTG", "AF_R_ITG", "AF_R_TP")
  
  # SLF II/III: gray matter termination measures
  var$AQ_SLF_GM_terminations  <- c("AQ_SLF_PMv", "AQ_SLF_BA44", "AQ_SLF_BA45", "AQ_SLF_SMG", "AQ_SLF_AG", "AQ_SLF_SPL")
  var$SLF_L_GM_terminations   <- c("SLF_L_PMv", "SLF_L_BA44", "SLF_L_BA45", "SLF_L_SMG", "SLF_L_AG", "SLF_L_SPL")
  var$SLF_R_GM_terminations   <- c("SLF_R_PMv", "SLF_R_BA44", "SLF_R_BA45", "SLF_R_SMG", "SLF_R_AG", "SLF_R_SPL")
  
  # Prior motor skill experience
  var$prior_experience <- c("gross_motor_experience", "fine_motor_experience")
  
  # Language measures (using renamed variables)
  var$language_measures <- c("syntactic_complexity", "AGL_d_grammaticality", "AGL_d_chunk_strength")
  
  # Toolmaking performance 
  var$toolmaking_performance <- c("toolmaking_performance")
  
  # PCA components (using renamed variables)
  var$pca_measures <- c("PC_thin", "PC_long")
  
  
  # Canonical display order for tables and figures (single source of truth)
  #===============================================================================
  # ROI label order per tract (for axis labels and table row order)
  roi_order <- list(
    AF   = c("PMv", "BA44", "BA45", "STG", "MTG", "ITG", "TP"),
    SLF  = c("PMv", "BA44", "BA45", "SPL", "AG", "SMG")
  )
  
  # WM tract order when reporting (AF, SLF)
  tract_order <- c("AF", "SLF")
  
  # Variable names in display order: AQ GM terminations 
  display_order <- list(
    AQ_WM      = var$AQ_WM_measures,
    AQ_AF_GM   = var$AQ_AF_GM_terminations,
    AQ_SLF_GM  = var$AQ_SLF_GM_terminations
  )
  
  # Define families of brain variables for analyses
  #===============================================================================
  brain_families <- list(
    AQ_WM   = var$AQ_WM_measures,
    AQ_AF   = var$AQ_AF_GM_terminations,
    AQ_SLF  = var$AQ_SLF_GM_terminations,
    WM      = var$WM_measures,
    AF_L    = var$AF_L_GM_terminations,
    AF_R    = var$AF_R_GM_terminations,
    SLF_L   = var$SLF_L_GM_terminations,
    SLF_R   = var$SLF_R_GM_terminations
  )
  
  # Define color sets frequently used in figures
  #===============================================================================
  color <- list()
  
  color$AF                <- "#D37676"
  color$pre$AF_L          <- "#D37676"
  color$pre$AF_R          <- "#F2C1C1"
  color$post$AF_L         <- "#D37676"
  color$post$AF_R         <- "#F2C1C1"
  
  color$SLF               <- "#08A8EB"
  color$pre$SLF_L         <- "#08A8EB"
  color$pre$SLF_R         <- "#A1E4F7"
  color$post$SLF_L        <- "#08A8EB"
  color$post$SLF_R        <- "#A1E4F7"
  
  
  color$exp               <- "#E79805"
  color$pre$exp           <- "#E79805"
  color$post$exp          <- "#FAD77B"
  
  color$con               <- "gray30"
  color$pre$con           <- "gray30"
  color$post$con          <- "gray70"
  
  color$correlation$low     <- "#039EE3"
  color$correlation$mid     <- "#FEFFFE" 
  color$correlation$high    <- "#FF0000"
  
  # Print a quick data summary to the console
  cat("Quick data summary:\n")
  cat(" - Total subjects:", n_distinct(data$subject), "\n")
  cat(" - Experimental group:", n_distinct(data$subject[data$group == "exp"]), "\n")
  cat(" - Control group:", n_distinct(data$subject[data$group == "con"]), "\n")
  cat(" - Pre-training observations:", sum(data$training == "pre"), "\n")
  cat(" - Post-training observations:", sum(data$training == "post"), "\n")
  
  
  # Save preprocessed data
  #===============================================================================
  write.csv(data, file = here("data", "preprocessed_data.csv"), row.names = FALSE)
  
  # Save preprocessed data and display-order config to data folder
  #===============================================================================
  save(data, var, brain_families, color, roi_order, tract_order, display_order,
       file = here("data/preprocessed_data.RData"))

} else {
  cat("Preprocessed data found — loading it.\n")
  load(here("data/preprocessed_data.RData"))
}

