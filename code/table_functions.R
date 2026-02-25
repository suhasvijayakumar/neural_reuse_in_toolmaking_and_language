
# Constants
#===============================================================================
DEFAULT_FDR_THRESHOLD <- 0.25  # Benjamini-Hochberg FDR threshold

table_toolmaking_brain_pcorr <- function(x,
                                caption = NULL,
                                q_star = DEFAULT_FDR_THRESHOLD) {
  
  # If x is a list, combine it
  if (is.list(x) && !is.data.frame(x)) {
    df <- dplyr::bind_rows(x)
  } else {
    df <- x
  }
  
  # Format table
  out <- df %>%
    dplyr::mutate(
      `Brain Measure` = gsub("_", " ", variable),
      `Partial r` = round(r, 3),
      `p value` = ifelse(is.na(p), "NA", sprintf("%.3f", p)),
      `BH q value` = ifelse(
        is.na(p_adj),
        "NA",
        ifelse(sig,
               paste0(sprintf("%.3f", p_adj), " *"),
               sprintf("%.3f", p_adj))
      )
    )
  
  # Select columns
  out <- out[, c("Brain Measure", "Partial r", "n", "p value", "BH q value")]
  
  # Return tibble if no caption
  if (is.null(caption)) return(out)
  
  # Otherwise print kable
  knitr::kable(out, caption = caption, booktabs = TRUE, escape = FALSE) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("hover", "condensed", "responsive"),
      full_width = TRUE,
      position = "center"
    ) %>%
    kableExtra::footnote(
      general = paste0(
        "Asterisk (*) marks BH-corrected results (q < ",
        q_star,
        ")."
      )
    )
}



#===========================
corr_to_language_table_df <- function(corr_result, q_star = DEFAULT_FDR_THRESHOLD) {
  
  r   <- corr_result$r
  p   <- corr_result$p
  q   <- corr_result$p.adj
  n   <- corr_result$n
  
  # If n is a single number, replicate to a matrix like r
  if (length(n) == 1) {
    n <- matrix(n, nrow = nrow(r), ncol = ncol(r), dimnames = dimnames(r))
  }
  
  vars <- colnames(r)
  langs <- rownames(r)  # should be your 3 language measures
  
  # helper to format q with star
  fmt_q <- function(x) {
    ifelse(is.na(x), "NA",
           ifelse(x < q_star, paste0(sprintf("%.3f", x), " *"),
                  sprintf("%.3f", x)))
  }
  
  # start table
  out <- data.frame(Measure = vars, stringsAsFactors = FALSE)
  
  for (lang in langs) {
    out[[paste0(lang, "_r")]] <- round(r[lang, ], 3)
    out[[paste0(lang, "_n")]] <- as.integer(n[lang, ])
    out[[paste0(lang, "_p")]] <- ifelse(is.na(p[lang, ]), "NA", sprintf("%.3f", p[lang, ]))
    out[[paste0(lang, "_q")]] <- fmt_q(q[lang, ])
  }
  
  out$Measure <- gsub("_", " ", out$Measure)
  out
}


table_language_brain_corr <- function(corr_result,
                                      caption = "",
                                      q_star = DEFAULT_FDR_THRESHOLD,
                                      language_headers = NULL) {
  
  df <- corr_to_language_table_df(corr_result, q_star = q_star)
  
  # how many language measures are in corr_result (usually 3)
  langs <- rownames(corr_result$r)
  
  knitr::kable(
    df,
    format = "html",
    escape = FALSE,
    align = c("l","c","c","c","l","c","c","c","l","c","c","c","l"),
    col.names = c("Measure", rep(c("Pearson's r", "n", "p-value", "BH q-value"), length(langs))),
    caption = caption
  ) %>%
    kableExtra::add_header_above(
      c(" " = 1, stats::setNames(rep(4, length(langs)), language_headers)),
      bold = TRUE
    ) %>%
    kableExtra::kable_styling(bootstrap_options = c("hover", "responsive")) %>%
    kableExtra::kable_classic(full_width = FALSE) %>%
    kableExtra::footnote(general = paste0("Asterisk (*) marks BH-adjusted p < ", q_star, "."))
}


#===============================================================================
# TABLE
#===============================================================================

table_anova_results <- function(anova_df,
                                caption = NULL,
                                fdr_level = DEFAULT_FDR_THRESHOLD,
                                digits = 3) {
  
  df <- anova_df
  
  # Round numeric columns (including etaG2_* and p's)
  num_cols <- c("F_group", "p_group", "etaG2_group",
                "F_training", "p_training", "etaG2_training",
                "F_interaction", "p_interaction", "etaG2_interaction",
                "adj_p_value")
  
  df[num_cols] <- lapply(df[num_cols], function(x) round(x, digits))
  
  # Add star to adjusted p-values (after rounding)
  df$adj_p_value <- ifelse(
    is.na(df$adj_p_value),
    "NA",
    ifelse(df$sig_25FDR,
           paste0(sprintf(paste0("%.", digits, "f"), df$adj_p_value), " *"),
           sprintf(paste0("%.", digits, "f"), df$adj_p_value))
  )
  
  df$variable <- gsub("_", " ", df$variable)
  
  # Reorder columns to match desired layout
  df_out <- df[, c(
    "variable",
    "F_group", "p_group", "etaG2_group",
    "F_training", "p_training", "etaG2_training",
    "F_interaction", "n_interaction",
    "p_interaction", "etaG2_interaction",
    "adj_p_value"
  )]
  
  knitr::kable(
    df_out,
    format = "html",
    escape = FALSE,
    align = "r",
    col.names = c(
      "Variable",
      "F", "p-value", "$\\eta G^2$",
      "F", "p-value", "$\\eta G^2$",
      "F", "n", "p-value", "$\\eta G^2$", "adj p-value"
    ),
    caption = caption
  ) %>%
    kableExtra::add_header_above(
      c(" " = 1,
        "Effect: group" = 3,
        "Effect: training" = 3,
        "Effect: group × training" = 5),
      bold = TRUE
    ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("hover", "responsive", "condensed")
    ) %>%
    kableExtra::kable_classic(full_width = FALSE) %>%
    kableExtra::footnote(
      general = paste0("Asterisk (*) marks BH-adjusted p < ", fdr_level, ".")
    )
}


table_prepost_pcorr_changes <- function(x,
                                        caption = NULL,
                                        show_family = FALSE,
                                        q_star = DEFAULT_FDR_THRESHOLD) {
  
  # If x is a list, combine it
  if (is.list(x) && !is.data.frame(x)) {
    df <- dplyr::bind_rows(x)
  } else {
    df <- x
  }
  
  out <- df %>%
    dplyr::mutate(
      Family   = gsub("_", " ", family),
      Variable = gsub("_", " ", variable),
      
      `Partial r (pre)`  = round(pre_r, 3),
      `p-value (pre)`    = ifelse(is.na(pre_p), "NA", sprintf("%.3f", pre_p)),
      
      `Partial r (post)` = round(post_r, 3),
      `p-value (post)`   = ifelse(is.na(post_p), "NA", sprintf("%.3f", post_p)),
      
      `Δr (post−pre)`    = round(delta_r, 3),
      n                  = n,
      
      `Z (Steiger)`      = ifelse(is.na(z), "NA", sprintf("%.3f", z)),
      `p-value (change)` = ifelse(is.na(p), "NA", sprintf("%.3f", p)),
      
      `BH q-value` = ifelse(
        is.na(p_adj),
        "NA",
        ifelse(sig,
               paste0(sprintf("%.3f", p_adj), " *"),
               sprintf("%.3f", p_adj))
      )
    )
  
  # Select columns (preserves incoming order)
  if (show_family) {
    out <- out[, c("Family",
                   "Variable",
                   "Partial r (pre)", "p-value (pre)",
                   "Partial r (post)", "p-value (post)",
                   "Δr (post−pre)", "n", "Z (Steiger)", "p-value (change)",
                   "BH q-value")]
    align_like_this <- c("l","l","c","c","c","c","c","c","c","c","l")
    n_left <- 2
  } else {
    out <- out[, c("Variable",
                   "Partial r (pre)", "p-value (pre)",
                   "Partial r (post)", "p-value (post)",
                   "Δr (post−pre)", "n", "Z (Steiger)", "p-value (change)",
                   "BH q-value")]
    align_like_this <- c("l","c","c","c","c","c","c","c","c","l")
    n_left <- 1
  }
  
  # Return tibble if no caption
  if (is.null(caption)) return(out)
  
  # Alignment goes HERE (knitr::kable), not in kable_styling
  kb <- knitr::kable(out, caption = caption, booktabs = TRUE, escape = FALSE,
                     align = align_like_this)
  
  if (requireNamespace("kableExtra", quietly = TRUE)) {
    kb <- kb %>%
      kableExtra::kable_styling(
        bootstrap_options = c("hover", "condensed", "responsive"),
        full_width = TRUE,
        position = "center"
      ) %>%
      kableExtra::add_header_above(c(
        " " = n_left,
        "Pre-training" = 2,
        "Post-training" = 2,
        "Change in association" = 5
      )) %>%
      kableExtra::footnote(
        general = paste0("Asterisk (*) marks BH-corrected results (q < ", q_star, ").")
      )
    return(kb)
  }
  
  kb
}





# Helper
corr_test_to_long_pca <- function(corr,
                                  pca_order = c("PC_thin","PC_long"),
                                  brain_order = NULL) {
  
  rmat <- corr$r
  pmat <- corr$p
  qmat <- corr$p.adj
  nmat <- corr$n
  
  df_long <- as.data.frame(as.table(rmat), stringsAsFactors = FALSE)
  names(df_long) <- c("pc", "measure", "r")
  
  df_long$p     <- as.vector(pmat)
  df_long$p_adj <- as.vector(qmat)
  df_long$n     <- as.vector(nmat)
  
  df_long <- tibble::as_tibble(df_long) %>%
    dplyr::mutate(
      pc = as.character(pc),
      measure = as.character(measure),
      r = as.numeric(r),
      p = as.numeric(p),
      p_adj = as.numeric(p_adj),
      n = as.integer(n)
    )
  
  if (!is.null(pca_order)) {
    df_long <- df_long %>% dplyr::mutate(pc = factor(pc, levels = pca_order))
  }
  if (!is.null(brain_order)) {
    df_long <- df_long %>% dplyr::mutate(measure = factor(measure, levels = brain_order))
  }
  
  df_long
}

table_pca_brain_corr <- function(corr,
                                 caption = NULL,
                                 q_star = DEFAULT_FDR_THRESHOLD,
                                 brain_order = NULL,
                                 pca_order = c("PC_thin","PC_long")) {
  
  df_long <- corr_test_to_long_pca(corr, pca_order = pca_order, brain_order = brain_order) %>%
    dplyr::mutate(
      Measure = gsub("_"," ", as.character(measure)),
      r_str   = ifelse(is.na(r), "NA", sprintf("%.3f", as.numeric(r))),
      n_str   = as.character(n),
      p_str   = ifelse(is.na(p), "NA", sprintf("%.3f", as.numeric(p))),
      q_str   = ifelse(
        is.na(p_adj), "NA",
        ifelse(as.numeric(p_adj) < q_star,
               paste0(sprintf("%.3f", as.numeric(p_adj)), " *"),
               sprintf("%.3f", as.numeric(p_adj)))
      )
    )
  
  wide <- df_long %>%
    dplyr::select(Measure, pc, r_str, n_str, p_str, q_str) %>%
    tidyr::pivot_wider(
      names_from = pc,
      values_from = c(r_str, n_str, p_str, q_str),
      names_glue = "{pc}__{.value}"
    )
  
  # Rename columns for display (keep unique names; duplicates get trailing spaces)
  rename_map <- c(
    "PC_thin__r_str" = "Pearson's r", "PC_thin__n_str" = "n", "PC_thin__p_str" = "p-value", "PC_thin__q_str" = "BH q-value",
    "PC_long__r_str" = "Pearson's r ", "PC_long__n_str" = "n ", "PC_long__p_str" = "p-value ", "PC_long__q_str" = "BH q-value "
  )
  for (nm in names(rename_map)) if (nm %in% names(wide)) names(wide)[names(wide) == nm] <- rename_map[[nm]]
  
  out <- wide[, c("Measure",
                  "Pearson's r","n","p-value","BH q-value",
                  "Pearson's r ","n ","p-value ","BH q-value ")]
  
  if (is.null(caption)) return(out)
  
  kb <- knitr::kable(out, caption = caption, booktabs = TRUE, escape = FALSE,
                     align = c("l", rep("c", 8)))
  
  if (requireNamespace("kableExtra", quietly = TRUE)) {
    kb <- kb %>%
      kableExtra::kable_styling(
        bootstrap_options = c("hover","condensed","responsive"),
        full_width = TRUE,
        position = "center"
      ) %>%
      kableExtra::add_header_above(c(
        " " = 1,
        "PC thin (corss-sectional thinning)" = 4,
        "PC long (elongation and pointedness)" = 4
      )) %>%
      kableExtra::footnote(
        general = paste0("Asterisk (*) marks BH-corrected results (q < ", q_star, ").")
      )
    return(kb)
  }
  
  kb
}



#===============================================================================
# tables.R
# Wide PC table for PCA-prepost correlation change results
# Input: tidy long df from run_pca_families_prepost() or calculate_*()
#===============================================================================
# Build wide table: variable rows, PC blocks (pre/post/comparison)
pca_prepost_as_wide_by_pc <- function(tidy_res,
                                      pc_order = c("PC_thin","PC_long")) {
  
  pc_order <- pc_order[pc_order %in% unique(tidy_res$pc)]
  
  wide <- tidy_res %>%
    dplyr::select(family, variable, pc,
                  pre_r, pre_n, pre_p,
                  post_r, post_n, post_p,
                  z, z_p, z_p_BH) %>%
    tidyr::pivot_wider(
      names_from = pc,
      values_from = c(pre_r, pre_n, pre_p,
                      post_r, post_n, post_p,
                      z, z_p, z_p_BH),
      names_sep = "_"
    )
  
  # enforce column order: variable + blocks per PC
  block_cols <- unlist(lapply(pc_order, function(pc) {
    c(paste0("pre_r_", pc),  paste0("pre_n_", pc),  paste0("pre_p_", pc),
      paste0("post_r_", pc), paste0("post_n_", pc), paste0("post_p_", pc),
      paste0("z_", pc),      paste0("z_p_", pc),    paste0("z_p_BH_", pc))
  }))
  
  keep <- c("variable", block_cols)
  keep <- keep[keep %in% names(wide)]
  wide[, keep, drop = FALSE]
}

print_pca_prepost_wide_table <- function(wide_df,
                                         caption = "",
                                         pc_order = c("PC_thin","PC_long"),
                                         digits = 2,
                                         q_star = DEFAULT_FDR_THRESHOLD) {
  
  # pretty variable names
  wide_df$variable <- gsub("_", " ", wide_df$variable)
  
  # --- add * to BH q columns -------------------------------------------------
  bh_cols <- intersect(c("z_p_BH_PC_thin", "z_p_BH_PC_long"), names(wide_df))
  
  for (cc in bh_cols) {
    x_num <- suppressWarnings(as.numeric(wide_df[[cc]]))
    
    wide_df[[cc]] <- dplyr::if_else(
      is.na(x_num),
      "NA",
      dplyr::if_else(
        x_num < q_star,
        paste0(sprintf(paste0("%.", digits, "f"), x_num), " *"),
        sprintf(paste0("%.", digits, "f"), x_num)
      )
    )
  }
  
  # column names per block
  col_block <- rep(c("r","n","p",
                     "r","n","p",
                     "Z","p","BH q"),
                   length(pc_order))
  
  knitr::kable(
    wide_df,
    digits = digits,
    format = "html",
    escape = FALSE,
    booktabs = TRUE,
    align = c("l", rep("c", 8), "l", rep("c", 8), "l"),
    col.names = c("variable", col_block),
    caption = caption
  ) %>%
    kableExtra::add_header_above(
      bold = TRUE,
      c(" " = 1,
        "pre-training"  = 3, "post-training" = 3, "comparison" = 3,
        "pre-training"  = 3, "post-training" = 3, "comparison" = 3)
    ) %>%
    kableExtra::add_header_above(
      bold = TRUE,
      c(" " = 1, "PC thin (cross-sectional thinning)" = 9, "PC long (elongation and pointedness)" = 9)
    ) %>%
    kableExtra::kable_styling(bootstrap_options = c("hover","responsive")) %>%
    kableExtra::kable_classic(full_width = FALSE) %>%
    kableExtra::footnote(
      general = paste0("Asterisk (*) marks BH-corrected results (q < ", q_star, ").")
    )
}
