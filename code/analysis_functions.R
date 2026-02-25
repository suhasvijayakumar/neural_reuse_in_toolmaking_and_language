#===============================================================================
# ANALYSIS FUNCTIONS
#===============================================================================
# This script contains all statistical analysis functions for the language-toolmaking study
# Author: Suhas Vijayakumar
# Date: 

# Constants
#===============================================================================
DEFAULT_FDR_THRESHOLD <- 0.25  # Benjamini-Hochberg FDR threshold
DEFAULT_ALPHA_NORM <- 0.05      # Alpha level for normality tests
DEFAULT_MU <- 0                 # Default mean for one-sample tests

# Load required libraries
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(stats)


# HELPERS:

#' Remove underscores and replace them with spaces
#'
#' @param df Data frame with column names, row names, and character values to process
#' @return Data frame with underscores replaced by spaces in column names, row names, and character columns
underscore_to_space <- function(df) {
  colnames(df) <- gsub("_", " ", colnames(df))
  rownames(df) <- gsub("_", " ", rownames(df))
  
  df_modified <- df %>%
    mutate(across(where(is.character), ~ gsub("_", " ", .)))
  
  return(df_modified)
}

#' Check for normality and run appropriate statistical test
#'
#' For each column in the data frame, performs Shapiro-Wilk normality test.
#' If data is normal (p >= alpha_norm), performs one-sample t-test.
#' Otherwise, performs Wilcoxon signed-rank test.
#'
#' @param df Data frame with numeric columns to test
#' @param mu Population mean for one-sample tests (default: 0)
#' @param alpha_norm Alpha level for normality test (default: 0.05)
#' @return Data frame with columns: Variable, n, Shapiro_W, Shapiro_p, Test, Statistic, p_value
check_normality_and_test <- function(df, mu = DEFAULT_MU, alpha_norm = DEFAULT_ALPHA_NORM) {
  cols <- colnames(df)
  
  map_dfr(cols, function(col) {
    x <- df[[col]]
    x <- x[!is.na(x)]
    n <- length(x)
    
    # guard: too few or constant values
    if (n < 3 || length(unique(x)) < 2) {
      return(tibble(
        Variable   = col,
        n          = n,
        Shapiro_W  = NA_real_,
        Shapiro_p  = NA_real_,
        Test       = NA_character_,
        Statistic  = NA_real_,
        p_value    = NA_real_
      ))
    }
    
    sh <- tryCatch(stats::shapiro.test(x), error = function(e) NULL)
    W   <- if (!is.null(sh)) unname(sh$statistic) else NA_real_
    pSh <- if (!is.null(sh)) sh$p.value            else NA_real_
    
    if (!is.na(pSh) && pSh >= alpha_norm) {
      tt  <- stats::t.test(x, mu = mu)
      tst <- "t-test"
      stat <- unname(tt$statistic)
      pval <- tt$p.value
    } else {
      wt  <- stats::wilcox.test(x, mu = mu, exact = FALSE)
      tst <- "Wilcoxon signed-rank"
      stat <- unname(wt$statistic)
      pval <- wt$p.value
    }
    
    tibble(
      Variable   = col,
      n          = n,
      Shapiro_W  = W,
      Shapiro_p  = pSh,
      Test       = tst,
      Statistic  = stat,
      p_value    = pval
    )
  })
}

#' Calculate partial correlations between toolmaking performance and brain measures
#'
#' Computes partial correlations controlling for a covariate (e.g., gross motor experience).
#' Uses residualization approach: residuals from regressing x and y on control are correlated.
#'
#' @param x Predictor variable (e.g., toolmaking_performance)
#' @param y Outcome variable(s) - can be vector or data frame of brain measures
#' @param control Control variable to partial out (e.g., gross_motor_experience)
#' @param method Correlation method (default: "pearson")
#' @param adjust_method Multiple comparison adjustment method (default: "BH" for Benjamini-Hochberg)
#' @param sig_q_value FDR threshold for significance (default: 0.25)
#' @return Data frame with columns: variable, r, n, p, p_adj, sig
calculate_toolmaking_brain_pcorr <- function(x, y, control, method = "pearson", adjust_method = "BH", sig_q_value = DEFAULT_FDR_THRESHOLD) {
  
  # Allow y to be vector or data frame (in case analyzing individual variable)
  if (is.data.frame(y)) {
    ys <- y
  } else {
    ys <- data.frame(outcome = y)
    names(ys) <- deparse(substitute(y))
  }
  
  outcome_names <- names(ys)
  
  one_outcome <- function(y_vec, outcome_name) {
    df <- data.frame(x = x, y = y_vec, control = control)
    df <- stats::na.omit(df)
    
    if (nrow(df) < 3 ||
        stats::sd(df$x) == 0 ||
        stats::sd(df$y) == 0 ||
        stats::sd(df$control) == 0) {
      return(tibble::tibble(
        variable = outcome_name,
        r = NA_real_,
        n = nrow(df),
        p = NA_real_
      ))
    }
    
    rx <- stats::resid(stats::lm(x ~ control, data = df))
    ry <- stats::resid(stats::lm(y ~ control, data = df))
    
    ct <- suppressWarnings(stats::cor.test(rx, ry, method = method))
    
    tibble::tibble(
      variable = outcome_name,
      r = unname(ct$estimate),
      n = length(rx),
      p = ct$p.value
    )
  }
  
  Map(one_outcome, ys, outcome_names) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      p_adj = stats::p.adjust(p, method = adjust_method),
      sig   = !is.na(p_adj) & p_adj < sig_q_value
    )
}

#' Calculate correlations between language measures and brain measures
#'
#' Computes Pearson correlations between language variables and brain variables.
#' Uses psych::corr.test for pairwise correlations with BH adjustment.
#'
#' @param df Data frame containing both language and brain variables
#' @param language_vars Character vector of language variable names
#' @param brain_vars Character vector of brain variable names
#' @return psych::corr.test object with correlation matrices, p-values, and adjusted p-values
calculate_language_brain_corr <- function(df, language_vars = NA, brain_vars = NA) {
  
  lang  <- df[, language_vars, drop = FALSE]
  brain <- df[, brain_vars, drop = FALSE] 
  
  psych::corr.test(lang, brain, use = "pairwise", method = "pearson", adjust = "BH")
}

#===============================================================================
# ANALYSIS
#===============================================================================

#' Run mixed-design ANOVA for a single measure
#'
#' Performs mixed-design ANOVA with training (within) and group (between) factors.
#' Returns F-statistics, p-values, and generalized eta-squared for main effects and interaction.
#'
#' @param data Data frame with subject, group, training, and measure columns
#' @param measure Name of the dependent variable to analyze
#' @return Data frame with ANOVA results: variable, F_group, p_group, etaG2_group, F_training, p_training, etaG2_training, F_interaction, n_interaction, p_interaction, etaG2_interaction
mixed_design_anova_function <- function(data, measure) {
  if (!requireNamespace("afex", quietly = TRUE)) stop("Please install afex.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr.")
  
  df <- data %>%
    dplyr::filter(!is.na(.data[[measure]])) %>%
    dplyr::mutate(
      subject  = factor(subject),
      training = factor(training, levels = c("pre", "post")),
      group    = factor(group,    levels = c("exp", "con"))
    )
  
  if (!all(c("pre","post") %in% unique(df$training)) ||
      !all(c("exp","con") %in% unique(df$group))) {
    return(data.frame(
      variable = measure,
      F_group=NA, p_group=NA, etaG2_group=NA,
      F_training=NA, p_training=NA, etaG2_training=NA,
      F_interaction=NA, n_interaction=NA,
      p_interaction=NA, etaG2_interaction=NA
    ))
  }
  
  complete_ids <- df %>%
    dplyr::group_by(group, subject) %>%
    dplyr::filter(dplyr::n_distinct(training) == 2) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(subject, group)
  
  n_interaction <- nrow(complete_ids)
  
  if (n_interaction < 4) {
    return(data.frame(
      variable = measure,
      F_group=NA, p_group=NA, etaG2_group=NA,
      F_training=NA, p_training=NA, etaG2_training=NA,
      F_interaction=NA, n_interaction=n_interaction,
      p_interaction=NA, etaG2_interaction=NA
    ))
  }
  
  df_model <- df %>% dplyr::semi_join(complete_ids, by = c("subject","group"))
  
  fit <- afex::aov_ez(
    id = "subject",
    dv = measure,
    data = df_model,
    within = "training",
    between = "group",
    type = 3
  )
  
  # Use raw numeric anova table
  tab <- as.data.frame(fit$anova_table)
  tab$term <- rownames(tab)
  
  get_val <- function(effect, col) {
    v <- tab[tab$term == effect, col]
    if (length(v) == 0) NA_real_ else as.numeric(v[1])
  }
  
  F_int <- get_val("group:training", "F")
  if (is.na(F_int)) F_int <- get_val("training:group", "F")
  
  p_int <- get_val("group:training", "Pr(>F)")
  if (is.na(p_int)) p_int <- get_val("training:group", "Pr(>F)")
  
  e_int <- get_val("group:training", "ges")
  if (is.na(e_int)) e_int <- get_val("training:group", "ges")
  
  data.frame(
    variable = measure,
    F_group = get_val("group", "F"),
    p_group = get_val("group", "Pr(>F)"),
    etaG2_group = get_val("group", "ges"),
    F_training = get_val("training", "F"),
    p_training = get_val("training", "Pr(>F)"),
    etaG2_training = get_val("training", "ges"),
    F_interaction = F_int,
    n_interaction = n_interaction,
    p_interaction = p_int,
    etaG2_interaction = e_int
  )
}



#' Run mixed-design ANOVA for multiple measures
#'
#' Wrapper function that runs mixed_design_anova_function for each measure and combines results.
#' Applies Benjamini-Hochberg FDR correction to interaction p-values.
#'
#' @param data Data frame with subject, group, training, and measure columns
#' @param measures Character vector of dependent variable names to analyze
#' @param fdr_level FDR threshold for significance (default: 0.25)
#' @param digits Number of decimal places for rounding (default: 3)
#' @return Data frame with ANOVA results for all measures, including adjusted p-values
run_mixed_design_anova <- function(data, measures, fdr_level = DEFAULT_FDR_THRESHOLD, digits = 3) {
  res <- purrr::map_dfr(measures, ~mixed_design_anova_function(data, .x))
  
  res <- res %>%
    dplyr::mutate(
      adj_p_value = p.adjust(p_interaction, method = "BH"),
      sig_25FDR   = !is.na(adj_p_value) & adj_p_value < fdr_level
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(dplyr::starts_with("F_"), dplyr::starts_with("p_"),
          dplyr::starts_with("eta2_"), adj_p_value),
        ~round(., digits)
      )
    )
  
  res
}

#===============================================================================
# ANALYSIS: delta (post - pre) follow-up for a mixed ANOVA interaction
#===============================================================================
anova_followup_delta <- function(data, variable,
                                 id_col = "subject",
                                 group_col = "group",
                                 time_col = "training",
                                 pre_level = "pre",
                                 post_level = "post") {
  
  df <- data %>%
    dplyr::select(
      subject = dplyr::all_of(id_col),
      group   = dplyr::all_of(group_col),
      time    = dplyr::all_of(time_col),
      value   = dplyr::all_of(variable)
    ) %>%
    dplyr::filter(!is.na(value), time %in% c(pre_level, post_level), group %in% c("exp", "con")) %>%
    dplyr::mutate(time = as.character(time))
  
  # Make wide (pre/post per subject)
  wide <- df %>%
    tidyr::pivot_wider(names_from = time, values_from = value) %>%
    dplyr::filter(!is.na(.data[[pre_level]]), !is.na(.data[[post_level]])) %>%
    dplyr::mutate(delta = .data[[post_level]] - .data[[pre_level]])
  
  # Main interaction follow-up: compare deltas between groups
  t_delta <- stats::t.test(delta ~ group, data = wide)
  
  # Helpful descriptive output
  summary_by_group <- wide %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_delta = mean(delta, na.rm = TRUE),
      sd_delta = stats::sd(delta, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    variable = variable,
    delta_df = wide,
    summary_by_group = summary_by_group,
    t_delta = t_delta
  )
}

#===============================================================================
# PRE-POST CORRELATION COMPARISON FUNCTIONS (WITH STEIGER'S Z)
#===============================================================================

.partial_residuals <- function(x, y, control) {
  df <- data.frame(x = x, y = y, control = control)
  df <- df[stats::complete.cases(df), , drop = FALSE]
  if (nrow(df) < 4L || stats::sd(df$x) == 0 || stats::sd(df$y) == 0) {
    return(list(x = numeric(0), y = numeric(0)))
  }
  list(
    x = stats::resid(stats::lm(x ~ control, data = df)),
    y = stats::resid(stats::lm(y ~ control, data = df))
  )
}

.get_z <- function(component) {
  if (is.null(component)) return(NA_real_)
  if (!is.null(component$statistic)) {
    zcand <- component$statistic
  } else if (!is.null(component$z.value)) {
    zcand <- component$z.value
  } else {
    zcand <- component$z
  }
  z <- suppressWarnings(as.numeric(zcand))
  if (length(z) == 1L && is.finite(z)) return(unname(z))
  NA_real_
}

partial_prepost_compare_one <- function(df, var, xvar = "toolmaking_performance",
                                        control = "gross_motor_experience",
                                        time = "training", pre_level = "pre", post_level = "post") {
  keep <- c("subject", time, xvar, control, var)
  df1 <- df[, keep, drop = FALSE]
  
  # Keep only pre/post rows and complete cases for required vars
  df1 <- df1 %>%
    dplyr::filter(.data[[time]] %in% c(pre_level, post_level)) %>%
    dplyr::filter(!is.na(.data[[xvar]]) & !is.na(.data[[var]]) & !is.na(.data[[control]]))
  
  # Subjects with BOTH pre and post
  have_both <- df1 %>%
    dplyr::group_by(subject) %>%
    dplyr::summarise(
      has_pre  = any(.data[[time]] == pre_level),
      has_post = any(.data[[time]] == post_level),
      .groups = "drop"
    ) %>%
    dplyr::filter(has_pre & has_post) %>%
    dplyr::pull(subject)
  
  n_use <- length(have_both)
  
  if (n_use < 4L) {
    return(tibble::tibble(
      variable = var,
      pre_r = NA_real_, pre_n = n_use, pre_p = NA_real_,
      post_r = NA_real_, post_n = n_use, post_p = NA_real_,
      z = NA_real_, p = NA_real_, n = n_use, test = NA_character_
    ))
  }
  
  # Split and align PRE/POST by subject
  dpre <- df1 %>%
    dplyr::filter(subject %in% have_both, .data[[time]] == pre_level) %>%
    dplyr::arrange(subject)
  
  dpost <- df1 %>%
    dplyr::filter(subject %in% have_both, .data[[time]] == post_level) %>%
    dplyr::arrange(subject)
  
  # Ensure exact same subject order in both
  common_subj <- intersect(dpre$subject, dpost$subject)
  dpre  <- dpre  %>% dplyr::filter(subject %in% common_subj) %>% dplyr::arrange(subject)
  dpost <- dpost %>% dplyr::filter(subject %in% common_subj) %>% dplyr::arrange(subject)
  
  n_use <- length(common_subj)
  if (n_use < 4L) {
    return(tibble::tibble(
      variable = var,
      pre_r = NA_real_, pre_n = n_use, pre_p = NA_real_,
      post_r = NA_real_, post_n = n_use, post_p = NA_real_,
      z = NA_real_, p = NA_real_, n = n_use, test = NA_character_
    ))
  }
  
  # Partial residuals (aligned subjects)
  pre_res  <- .partial_residuals(dpre[[xvar]],  dpre[[var]],  dpre[[control]])
  post_res <- .partial_residuals(dpost[[xvar]], dpost[[var]], dpost[[control]])
  
  # If residualization dropped rows (shouldn’t, but guard anyway)
  if (length(pre_res$x) < 4L || length(post_res$x) < 4L) {
    n_use2 <- min(length(pre_res$x), length(post_res$x), length(pre_res$y), length(post_res$y))
    return(tibble::tibble(
      variable = var,
      pre_r = NA_real_, pre_n = n_use2, pre_p = NA_real_,
      post_r = NA_real_, post_n = n_use2, post_p = NA_real_,
      z = NA_real_, p = NA_real_, n = n_use2, test = NA_character_
    ))
  }
  
  ax <- pre_res$x
  ay <- pre_res$y
  bx <- post_res$x
  by <- post_res$y
  
  # Correlation tests for pre and post (partial residual correlation)
  pre_ct  <- suppressWarnings(stats::cor.test(ax, ay))
  post_ct <- suppressWarnings(stats::cor.test(bx, by))
  pre_r  <- unname(pre_ct$estimate)
  pre_p  <- pre_ct$p.value
  post_r <- unname(post_ct$estimate)
  post_p <- post_ct$p.value
  
  # Six cross-correlations for Steiger non-overlap test
  r_ab <- suppressWarnings(stats::cor(ax, ay)) # pre target
  r_cd <- suppressWarnings(stats::cor(bx, by)) # post target
  r_ac <- suppressWarnings(stats::cor(ax, bx)) # Xpre-Xpost
  r_bd <- suppressWarnings(stats::cor(ay, by)) # Ypre-Ypost
  r_ad <- suppressWarnings(stats::cor(ax, by)) # Xpre-Ypost
  r_bc <- suppressWarnings(stats::cor(ay, bx)) # Ypre-Xpost
  
  z_val <- NA_real_
  p_val <- NA_real_
  test_used <- NA_character_
  
  if (requireNamespace("cocor", quietly = TRUE)) {
    cc <- try(
      cocor::cocor.dep.groups.nonoverlap(
        r.jk = r_ab, r.hm = r_cd,
        r.jh = r_ac, r.km = r_bd, r.jm = r_ad, r.kh = r_bc,
        n = n_use, alternative = "two.sided", test = "steiger1980"
      ),
      silent = TRUE
    )
    if (!inherits(cc, "try-error") && !is.null(cc@steiger1980$p.value) && !is.na(cc@steiger1980$p.value)) {
      z_val <- .get_z(cc@steiger1980)
      p_val <- cc@steiger1980$p.value
      test_used <- "Steiger z"
    } else {
      cc2 <- try(
        cocor::cocor.dep.groups.nonoverlap(
          r.jk = r_ab, r.hm = r_cd,
          r.jh = r_ac, r.km = r_bd, r.jm = r_ad, r.kh = r_bc,
          n = n_use, alternative = "two.sided", test = "fisher1925"
        ),
        silent = TRUE
      )
      if (!inherits(cc2, "try-error") && !is.null(cc2@fisher1925$p.value) && !is.na(cc2@fisher1925$p.value)) {
        z_val <- .get_z(cc2@fisher1925)
        p_val <- cc2@fisher1925$p.value
        test_used <- "Fisher z"
      }
    }
  }
  
  tibble::tibble(
    variable = var,
    pre_r = pre_r, pre_n = n_use, pre_p = pre_p,
    post_r = post_r, post_n = n_use, post_p = post_p,
    z = z_val, p = p_val, n = n_use, test = test_used
  )
}

partial_prepost_compare <- function(df, measures, xvar = "toolmaking_performance",
                                    control = "gross_motor_experience",
                                    time = "training", pre_level = "pre", post_level = "post") {
  measures <- as.character(unlist(measures))
  measures <- measures[measures %in% names(df)]
  
  if (length(measures) == 0) {
    return(tibble::tibble(
      variable = character(),
      pre_r = double(), pre_n = integer(), pre_p = double(),
      post_r = double(), post_n = integer(), post_p = double(),
      z = double(), p = double(), n = integer(), test = character()
    ))
  }
  
  purrr::map_dfr(measures, function(m) {
    out <- try(
      partial_prepost_compare_one(
        df = df, var = m, xvar = xvar, control = control,
        time = time, pre_level = pre_level, post_level = post_level
      ),
      silent = TRUE
    )
    if (inherits(out, "try-error")) {
      tibble::tibble(
        variable = m, pre_r = NA_real_, pre_n = NA_integer_, pre_p = NA_real_,
        post_r = NA_real_, post_n = NA_integer_, post_p = NA_real_,
        z = NA_real_, p = NA_real_, n = NA_integer_, test = NA_character_
      )
    } else {
      out
    }
  })
}


calculate_pca_brain_corr <- function(df, pca_vars = c("PC_thin","PC_long"), brain_vars = NA) {
  
  pca_vars   <- as.character(pca_vars)
  brain_vars <- as.character(brain_vars)
  
  pca_vars   <- pca_vars[pca_vars %in% names(df)]
  brain_vars <- brain_vars[brain_vars %in% names(df)]
  
  pca  <- df[, pca_vars,  drop = FALSE]
  brain <- df[, brain_vars, drop = FALSE]
  
  psych::corr.test(pca, brain, use = "pairwise", method = "pearson", adjust = "BH")
}


#===============================================================================
# analysis.R
# PCA (PC_thin / PC_long) × brain measures: PRE vs POST correlations + change test
# Change test: cocor.dep.groups.nonoverlap (Steiger 1980; Fisher fallback)
# Output: tidy long df with columns matching your plot/table pipeline
#===============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- helper: pull matched pre/post vectors for one predictor & one brain var ---
.match_prepost <- function(df, xvar, yvar,
                           id = "subject",
                           time = "training",
                           pre_level = "pre",
                           post_level = "post") {
  
  stopifnot(is.character(xvar), length(xvar) == 1,
            is.character(yvar), length(yvar) == 1,
            is.character(id),   length(id) == 1,
            is.character(time), length(time) == 1)
  
  need <- c(id, time, xvar, yvar)
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("Missing columns in df: ", paste(miss, collapse = ", "))
  }
  
  df1 <- df[, need, drop = FALSE]
  
  # subjects with BOTH pre & post and non-NA at BOTH timepoints for BOTH vars
  ok_ids <- df1 %>%
    dplyr::filter(.data[[time]] %in% c(pre_level, post_level)) %>%
    dplyr::group_by(.data[[id]]) %>%
    dplyr::filter(all(c(pre_level, post_level) %in% .data[[time]])) %>%
    dplyr::filter(!any(is.na(.data[[xvar]]) | is.na(.data[[yvar]]))) %>%
    dplyr::pull(.data[[id]]) %>%
    unique()
  
  if (!length(ok_ids)) {
    return(list(A = numeric(0), B = numeric(0), C = numeric(0), D = numeric(0)))
  }
  
  df2 <- df1 %>%
    dplyr::filter(.data[[id]] %in% ok_ids, .data[[time]] %in% c(pre_level, post_level)) %>%
    dplyr::arrange(.data[[id]], match(.data[[time]], c(pre_level, post_level)))
  
  # wide for clean pairing
  wide <- df2 %>%
    tidyr::pivot_wider(
      names_from  = dplyr::all_of(time),
      values_from = dplyr::all_of(c(xvar, yvar)),
      names_sep   = "_"
    )
  
  # A= x_pre, B= y_pre, C= x_post, D= y_post
  list(
    A = wide[[paste0(xvar, "_", pre_level)]],
    B = wide[[paste0(yvar, "_", pre_level)]],
    C = wide[[paste0(xvar, "_", post_level)]],
    D = wide[[paste0(yvar, "_", post_level)]]
  )
}


# --- helper: Pearson r and p (requires n >= 4 for stability) ------------------
.r_p <- function(u, v) {
  n <- min(length(u), length(v))
  if (!is.finite(n) || n < 4L) return(list(r = NA_real_, p = NA_real_))
  ct <- suppressWarnings(stats::cor.test(u, v, method = "pearson"))
  list(r = unname(ct$estimate), p = ct$p.value)
}

# --- core: Steiger dependent nonoverlap ---------------------------------------
# compares r(A,B) vs r(C,D) for same participants (pre vs post)
.steiger_dep_nonoverlap <- function(A, B, C, D) {
  n <- min(length(A), length(B), length(C), length(D))
  if (!is.finite(n) || n < 4L) return(list(z = NA_real_, p = NA_real_, test = "Steiger z", n = n))
  
  r_ab <- suppressWarnings(stats::cor(A, B, use = "pairwise.complete.obs"))
  r_cd <- suppressWarnings(stats::cor(C, D, use = "pairwise.complete.obs"))
  r_ac <- suppressWarnings(stats::cor(A, C, use = "pairwise.complete.obs"))
  r_bd <- suppressWarnings(stats::cor(B, D, use = "pairwise.complete.obs"))
  r_ad <- suppressWarnings(stats::cor(A, D, use = "pairwise.complete.obs"))
  r_bc <- suppressWarnings(stats::cor(B, C, use = "pairwise.complete.obs"))
  
  if (!requireNamespace("cocor", quietly = TRUE)) {
    return(list(z = NA_real_, p = NA_real_, test = "Steiger z", n = n))
  }
  
  cc <- try(
    cocor::cocor.dep.groups.nonoverlap(
      r.jk = r_ab, r.hm = r_cd,
      r.jh = r_ac, r.km = r_bd, r.jm = r_ad, r.kh = r_bc,
      n = n, alternative = "two.sided", test = "steiger1980"
    ),
    silent = TRUE
  )
  
  if (!inherits(cc, "try-error") &&
      !is.null(cc@steiger1980$p.value) &&
      is.finite(cc@steiger1980$p.value)) {
    return(list(z = .get_z(cc@steiger1980),
                p = cc@steiger1980$p.value,
                test = "Steiger z",
                n = n))
  }
  
  cc2 <- try(
    cocor::cocor.dep.groups.nonoverlap(
      r.jk = r_ab, r.hm = r_cd,
      r.jh = r_ac, r.km = r_bd, r.jm = r_ad, r.kh = r_bc,
      n = n, alternative = "two.sided", test = "fisher1925"
    ),
    silent = TRUE
  )
  
  if (!inherits(cc2, "try-error") &&
      !is.null(cc2@fisher1925$p.value) &&
      is.finite(cc2@fisher1925$p.value)) {
    return(list(z = .get_z(cc2@fisher1925),
                p = cc2@fisher1925$p.value,
                test = "Fisher z",
                n = n))
  }
  
  list(z = NA_real_, p = NA_real_, test = "Steiger z", n = n)
}

# --- MAIN: family-wise PC × brain change analysis (tidy long) -----------------
calculate_pca_brain_prepost_change <- function(df,
                                               pca_vars = c("PC_thin","PC_long"),
                                               brain_vars,
                                               id = "subject",
                                               time = "training",
                                               pre_level = "pre",
                                               post_level = "post",
                                               bh_thresh = DEFAULT_FDR_THRESHOLD) {
  
  pca_vars   <- as.character(pca_vars)
  pca_vars   <- pca_vars[pca_vars %in% names(df)]
  
  brain_vars <- as.character(brain_vars)
  brain_vars <- brain_vars[brain_vars %in% names(df)]
  
  if (!length(pca_vars) || !length(brain_vars)) return(tibble::tibble())
  
  res <- purrr::map_dfr(pca_vars, function(pc) {
    purrr::map_dfr(brain_vars, function(brain) {
      
      vecs <- .match_prepost(df, xvar = pc, yvar = brain,
                             id = id, time = time,
                             pre_level = pre_level, post_level = post_level)
      A <- vecs$A; B <- vecs$B; C <- vecs$C; D <- vecs$D
      
      n_subj <- min(length(A), length(B), length(C), length(D))
      
      pre  <- .r_p(A, B)
      post <- .r_p(C, D)
      cmp  <- .steiger_dep_nonoverlap(A, B, C, D)
      
      tibble::tibble(
        pc = pc,
        variable = brain,
        pre_r = pre$r,  pre_p = pre$p,  pre_n = n_subj,
        post_r = post$r, post_p = post$p, post_n = n_subj,
        delta_r = post$r - pre$r,
        n = n_subj,
        z = cmp$z,
        z_p = cmp$p,
        test = cmp$test
      )
    })
  })
  
  # BH within each PC (within the family call)
  res %>%
    dplyr::group_by(pc) %>%
    dplyr::mutate(
      z_p_BH = p.adjust(z_p, method = "BH"),
      sig = !is.na(z_p_BH) & z_p_BH < bh_thresh
    ) %>%
    dplyr::ungroup()
}

# convenience: run all named families, return one tidy long df
run_pca_families_prepost <- function(df, brain_families,
                                     pca_vars = c("PC_thin","PC_long"),
                                     id = "subject", time = "training",
                                     pre_level = "pre", post_level = "post",
                                     bh_thresh = DEFAULT_FDR_THRESHOLD) {
  stopifnot(is.list(brain_families), !is.null(names(brain_families)))
  
  purrr::imap_dfr(brain_families, function(vars, fam) {
    calculate_pca_brain_prepost_change(
      df = df,
      pca_vars = pca_vars,
      brain_vars = vars,
      id = id, time = time,
      pre_level = pre_level, post_level = post_level,
      bh_thresh = bh_thresh
    ) %>% dplyr::mutate(family = fam, .before = 1)
  })
}
