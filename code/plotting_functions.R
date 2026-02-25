#===============================================================================
# PLOTTING FUNCTIONS
#===============================================================================
# This script contains all visualization functions for the study
# Author: Suhas Vijayakumar

# Constants
#===============================================================================
DEFAULT_FDR_THRESHOLD <- 0.25  # Benjamini-Hochberg FDR threshold for significance
DEFAULT_BASE_SIZE <- 14        # Default base font size for plots

#===============================================================================
# ASYMMETRY QUOTIENT (AQ) VISUALIZATION FUNCTIONS
#===============================================================================

# Pivot wide AQ data to long and attach human-readable labels
pivot_AQ_long <- function(df_wide, vars, labels, out_name = "roi") {
  stopifnot(length(vars) == length(labels))
  out <- df_wide[, vars, drop = FALSE] |>
    pivot_longer(cols = everything(),
                 names_to = paste0("AQ_", out_name),
                 values_to = "AQ_values")
  out[[out_name]] <- factor(rep(labels, times = nrow(df_wide)), levels = labels)
  out
}


# Half-violin + half-point plot with mean±SEM, LR guides, consistent styling
plot_AQ <- function(dflong, x_var = "roi", y_var = "AQ_values",
                    title = "", fill_color = "gray50",
                    y_lim = c(-2.1, 2.1)) {
  
  n_groups <- nlevels(dflong[[x_var]])
  x_annot  <- 0.6
  
  ggplot(dflong, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#D3D4D8", linewidth = 1) +
    geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.9, shape = 4, size = 3.5, color = fill_color) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +
    stat_summary(fun = mean, geom = "point", colour = "black", shape = 10, size = 4, stroke = 1.2) +
    annotate("text", x = x_annot, y = y_lim[1] + 0.1, label = "L", fontface = "bold") +
    annotate("text", x = x_annot, y = y_lim[2] - 0.1, label = "R", fontface = "bold") +
    labs(title = title, x = "", y = "") +
    coord_flip() +
    ylim(y_lim) +
    theme_pubclean() +
    theme(legend.position = "none")
}

# Stack two aligned plots with proportional or fixed heights
stack_two <- function(p_top, p_bottom, top_n, bot_n, ratio = NULL) {
  aligned <- cowplot::align_plots(p_top, p_bottom, align = "v", axis = "l")
  rel_h   <- if (is.null(ratio)) c(top_n, bot_n) else c(1, ratio)
  cowplot::plot_grid(aligned[[1]], aligned[[2]], ncol = 1, rel_heights = rel_h)
}

#===============================================================================
# BRAIN MEASURE VISUALIZATION FUNCTIONS
#===============================================================================

# L/R rain + box + mean plot for white matter measures
plot_wm_measures <- function(df, vars, ylab, title, fill = "#D37676", ylim = NULL, size = NULL) {
  stopifnot(length(vars) == 2)
  
  long <- df %>%
    select(all_of(vars)) %>%
    pivot_longer(cols = all_of(vars), names_to = "colname", values_to = "value") %>%
    mutate(
      hemi = case_when(
        grepl("_L_", colname) | grepl("_L$", colname) ~ "L",
        grepl("_R_", colname) | grepl("_R$", colname) ~ "R",
        TRUE ~ colname
      ),
      hemi = factor(hemi, levels = c("L", "R"))
    )
  
  p <- ggplot(long, aes(hemi, value, fill = hemi)) +
    geom_rain(
      alpha = 0.9, rain.side = "f1x1", color = "#8D8680",
      boxplot.args = list(alpha = 1, color = "#1E1E1E", fill = "#FEFFFE", outlier.shape = NA),
      violin.args = list(color = NA)
    ) +
    geom_point(stat = "summary", fun = "mean", color = "#1E1E1E", shape = 10, size = 3, stroke = 1) +
    theme_pubclean(base_size = size) +
    scale_fill_manual(values = c(fill, fill)) +
    guides(fill = "none", color = "none") +
    labs(title = title, x = "hemisphere", y = ylab)
  
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# Box + mean plot for GM terminations L/R columns with ROI labels
plot_gm_terminations <- function(data_pre, left_cols, right_cols, roi_labels,
                                 title, fill_L, fill_R, base_size = 16) {
  stopifnot(length(left_cols) == length(right_cols),
            length(roi_labels) == length(left_cols))
  
  map <- tibble(
    col = c(left_cols, right_cols),
    hemi = c(rep("L", length(left_cols)), rep("R", length(right_cols))),
    roi  = rep(roi_labels, 2)
  )
  
  df_long <- data_pre %>%
    select(all_of(map$col)) %>%
    pivot_longer(cols = all_of(map$col), names_to = "col", values_to = "GM_terminations") %>%
    left_join(map, by = "col") %>%
    mutate(
      hemi = factor(hemi, levels = c("L","R")),
      roi  = factor(roi, levels = roi_labels)
    )
  
  ggplot(df_long, aes(x = roi, y = GM_terminations, fill = hemi)) +
    geom_boxplot(
      alpha = 0.85,
      position = position_dodge(width = 0.55),
      color = "#1E1E1E",
      width = 0.45,
      outlier.shape = NA
    ) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.55),
      alpha = 0.5, shape = 16, size = 2, stroke = 0.6, color = "#8D8680"
    ) +
    geom_point(
      stat = "summary", fun = "mean",
      colour = "#FEFFFE", shape = 10, size = 3, stroke = 1,
      position = position_dodge(width = 0.55)
    ) +
    labs(
      title = title,
      x = "regions of interest",
      y = expression(paste("volume (mm"^3,")"))
    ) +
    theme_pubclean(base_size = base_size) +
    scale_fill_manual(values = c("L" = fill_L, "R" = fill_R), name = "Hemisphere") +
    scale_color_manual(values = c("L" = fill_L, "R" = fill_R), guide = "none") + 
    guides(fill = guide_legend(override.aes = list(shape = NA))) +
    theme(legend.position = "top")
}

#=====================

plot_toolmaking_brain_heatmap <- function(df,
                               row_label = "toolmaking \nperformance",
                               title = "Partial correlations",
                               subtitle = NULL,
                               filename = NULL,
                               width = 8,
                               height = 2.5) {
  
  # Remove underscores and keep order
  df$variable <- gsub("_", " ", df$variable)
  df$variable <- factor(df$variable, levels = df$variable)
  
  # Label inside cell
  df$label <- paste0(
    sprintf("%.2f", df$r),
    "\n(n = ", df$n, ")"
  )
  
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = variable, y = row_label, fill = r)) +
    ggplot2::geom_tile(color = "grey90") +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      size = 4
    ) +
    ggplot2::geom_text(
      data = subset(df, sig),
      label = "*",
      nudge_y = 0.25,
      size = 6
    ) +
    ggplot2::scale_fill_gradient2(
      limits = c(-1, 1),
      low = color$correlation$low,
      mid = color$correlation$mid,
      high = color$correlation$high,
      midpoint = 0,
      name = "Partial r"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey30"),
      legend.title.position = "bottom"
    )
  
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = width, height = height, dpi = 300)
  }
  
  return(p)
}

plot_language_brain_heatmap <- function(corr,
                                        title = "Correlations",
                                        subtitle = NULL,
                                        filename = NULL,
                                        width = 10,
                                        height = 6,
                                        q_star = DEFAULT_FDR_THRESHOLD) {
  
  r <- corr$r
  p_adj <- corr$p.adj
  n <- corr$n
  
  # Make n a matrix if it is a single number
  if (length(n) == 1) {
    n <- matrix(n, nrow = nrow(r), ncol = ncol(r), dimnames = dimnames(r))
  }
  
  # Convert matrices -> long dataframe (keeps column order)
  df <- as.data.frame(as.table(r), stringsAsFactors = FALSE)
  names(df) <- c("language", "brain", "r")
  
  df$p_adj <- as.vector(p_adj)
  df$n     <- as.vector(n)
  df$sig   <- !is.na(df$p_adj) & df$p_adj < q_star
  
  # Remove underscores and keep brain order unchanged
  brain_levels <- colnames(r)
  df$brain <- factor(df$brain, levels = brain_levels)
  df$brain <- gsub("_", " ", as.character(df$brain))
  df$brain <- factor(df$brain, levels = gsub("_", " ", brain_levels))
  
  # Remove underscores from language labels too
  lang_levels <- rownames(r)
  df$language <- factor(df$language, levels = lang_levels)
  df$language <- gsub("_", " ", as.character(df$language))
  df$language <- factor(df$language, levels = gsub("_", " ", lang_levels))
  
  # Cell label (r + n)
  df$label <- paste0(sprintf("%.2f", df$r), "\n(n = ", df$n, ")")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = brain, y = language, fill = r)) +
    ggplot2::geom_tile(color = "grey90") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.5) +
    ggplot2::geom_text(
      data = subset(df, sig),
      label = "*",
      nudge_y = 0.25,
      size = 6
    ) +
    ggplot2::scale_fill_gradient2(
      limits = c(-1, 1),
      low = color$correlation$low,
      mid = color$correlation$mid,
      high = color$correlation$high,
      midpoint = 0,
      name = "Pearson's r"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey30"),
      legend.title.position = "bottom"
    )
  
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = width, height = height, dpi = 300)
  }
  
  return(p)
}


#===============================================================================
# PLOT
#===============================================================================

#===============================================================================
# FAMILY-WISE PLOTTING FUNCTIONS
#===============================================================================

# Create change correlation plots with arrows (dumbbell plots)
#===============================================================================
plot_prepost_pcorr_change_as_arrows <- function(df_family,
                                                title_prefix = "Pre-Post Correlation Changes:",
                                                subtitle = "Arrows show direction of change; * q<0.25 (BH)",
                                                lw_sig = 1.2,
                                                lw_nonsig = 0.7) {
  
  fam <- unique(df_family$family)
  if (length(fam) != 1) fam <- fam[1]
  
  ggplot2::ggplot(df_family, ggplot2::aes(y = variable_label)) +
    ggplot2::geom_point(ggplot2::aes(x = pre_r),
                        shape = 21, size = 3, fill = "gray70",
                        color = "black", stroke = 0.25, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(x = post_r),
                        shape = 16, size = 3,
                        color = "black", na.rm = TRUE) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = pre_r, xend = post_r, yend = variable_label,
        color = dir,
        linewidth = ifelse(sig, lw_sig, lw_nonsig)   # <-- weighted arrows
      ),
      arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed", ends = "last"),
      na.rm = TRUE
    ) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::geom_text(
      ggplot2::aes(x = (pre_r + post_r)/2, label = sprintf("Δr = %+.2f", delta_r)),
      size = 3.5, vjust = -0.6, show.legend = FALSE, na.rm = TRUE, color = "black"
    ) +
    ggplot2::geom_text(ggplot2::aes(x = -0.95, label = ifelse(sig, "*", "")),
                       size = 6, hjust = 0, color = "black", na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(x = 0.95, label = paste0("(n = ", n, ")")),
                       size = 3.5, hjust = 1, color = "gray35", na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::coord_cartesian(xlim = c(-1, 1)) +
    ggplot2::scale_color_manual(values = c(Increase = color$correlation$high,
                                           Decrease = color$correlation$low),
                                name = "Δr change") +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(
      title = paste(title_prefix, fam),
      subtitle = subtitle,
      x = "Partial r (–1 to 1)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "top",
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray50")
    )
}




#===============================================================================
# CORRELATION PLOTTING FUNCTIONS
#===============================================================================
# Heatmap with BRAIN MEASURES on the X-axis (left→right), predictors on Y
plot_corr_heatmap_faceted <- function(df,
                                      order_vec   = NULL,
                                      x_labels    = NULL,
                                      title       = "plot title",
                                      subtitle    = "plot subtitle",
                                      bh_thresh   = DEFAULT_FDR_THRESHOLD,
                                      show_vals   = TRUE,
                                      show_n      = TRUE,
                                      x_angle     = 60,
                                      palette     = c(low="#039EE3", mid="white", high="#FF0000")) {
  `%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
  
  .norm_corr_df <- function(df, bh_thresh = DEFAULT_FDR_THRESHOLD, predictor_fallback = "toolmaking") {
    df <- df %>%
      dplyr::mutate(
        r         = r %||% correlation %||% pre_r %||% post_r,
        p         = p %||% p_value %||% pre_p,
        p_adj     = p_adj %||% adj_p_value %||% p.adjust(p, method = "BH"),
        predictor = predictor %||% predictor_fallback,
        sig       = !is.na(p_adj) & p_adj < bh_thresh
      )
    df %>%
      dplyr::select(family, variable, predictor, r, n, p, p_adj, sig) %>%
      dplyr::distinct()
  }
  
  dd <- .norm_corr_df(df, bh_thresh = bh_thresh)
  
  if (!is.null(order_vec)) {
    dd <- dd %>% dplyr::mutate(variable = factor(variable, levels = order_vec))
  } else {
    # Only reorder if there are multiple variables in the family and no missing values
    dd <- dd %>%
      dplyr::group_by(family) %>%
      dplyr::mutate(
        variable = if (n() > 1 && all(!is.na(r))) {
          tryCatch({
            forcats::fct_reorder(variable, .x = abs(r), .fun = max)
          }, error = function(e) {
            variable
          })
        } else {
          variable
        }
      ) %>%
      dplyr::ungroup()
  }
  
  thm <- theme_minimal(base_size = 16)
  
  p <- ggplot2::ggplot(dd, ggplot2::aes(x = variable, y = predictor, fill = r)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.3) +
    { if (show_vals) ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 5) } +
    ggplot2::geom_text(data = dplyr::filter(dd, sig),
                       label = "*", nudge_x = 0, nudge_y = 0.25, size = 7, color = "black") +
    { if (show_n) ggplot2::geom_text(ggplot2::aes(label = sprintf("(n = %.0f)", n)), size = 4, nudge_y = -0.25) } +
    ggplot2::scale_fill_gradient2(limits = c(-1,1),
                                  low = palette["low"], mid = palette["mid"], high = palette["high"],
                                  name = "Partial r") +
    {
      if (!is.null(x_labels)) {
        ggplot2::scale_x_discrete(labels = x_labels)
      } else {
        ggplot2::scale_x_discrete(labels = ~ gsub("_"," ", .x))
      }
    } +
    ggplot2::scale_y_discrete(labels = ~ gsub("_","\n ", .x)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    thm +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = x_angle, vjust = 1, hjust = 1),
      legend.position = "right",
      legend.key.size = grid::unit(0.5, "cm"),
      legend.text     = element_text(size = 12),
      legend.title    = element_text(size = 10),
      plot.title = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle = element_text(size = rel(0.9), colour = "grey30", margin = margin(b = 6))
    )
  
  if ("family" %in% names(dd)) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(family), scales = "fixed")
  }
  
  p + ggplot2::coord_fixed()
}

#===============================================================================
# PRE-POST COMPARISON VISUALIZATION FUNCTIONS
#===============================================================================

# Plot pre-post change with paired t-tests + raincloud distribution
plot_prepost_change_raincloud <- function(data, measures, title = "Pre vs Post Training") {
  
  df_exp <- data %>% filter(group == "exp")
  
  df_long <- df_exp %>%
    dplyr::select(subject, training, all_of(measures)) %>%
    pivot_longer(
      cols = all_of(measures),
      names_to = "measure", values_to = "value"
    ) %>%
    dplyr::mutate(
      measure  = factor(measure, levels = measures, labels = gsub("_", " ", measures)),
      training = factor(training, levels = c("pre", "post"))
    )
  
  # Keep ONLY paired subjects
  paired_ids <- df_long %>%
    group_by(measure, subject) %>%
    summarise(
      has_pre  = any(training == "pre"  & !is.na(value)),
      has_post = any(training == "post" & !is.na(value)),
      .groups = "drop"
    ) %>%
    filter(has_pre & has_post) %>%
    dplyr::select(measure, subject)
  
  df_paired <- df_long %>%
    inner_join(paired_ids, by = c("measure", "subject"))
  
  # Paired t-tests (one per measure)
  ttests <- df_paired %>%
    dplyr::group_by(measure, subject) %>%
    dplyr::summarise(
      pre  = value[training == "pre"][1],
      post = value[training == "post"][1],
      .groups = "drop_last"
    ) %>%
    dplyr::summarise(
      t  = unname(t.test(post, pre, paired = TRUE)$statistic),
      df = unname(t.test(post, pre, paired = TRUE)$parameter),
      p  = t.test(post, pre, paired = TRUE)$p.value,
      mean_change = mean(post - pre, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_lab = formatC(p, format = "g", digits = 3))
  
  # Compute per-facet annotation y
  y_annot <- df_paired %>%
    group_by(measure) %>%
    summarise(y = max(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(y = y + 0.08 * ifelse(is.finite(y), abs(y), 1))
  
  ttests_annot <- ttests %>%
    left_join(y_annot, by = "measure") %>%
    mutate(
      label = paste0("paired t(", round(df, 0), ") = ", round(t, 2), ",\np = ", p_lab)
    )
  
  
  p <- ggplot(df_paired, aes(x = training, y = value, fill = training)) +
    ggdist::stat_halfeye(
      adjust = 0.7,
      width = 0.6,
      justification = -0.20,
      point_colour = NA,
      alpha = 0.5
    ) +
    
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
    geom_line(aes(group = subject), color = "grey50", alpha = 0.5) +
    geom_jitter(width = 0.10, alpha = 0.8, size = 2, shape = 25, fill = "gray30", color = "gray30") +
    stat_summary(fun = mean, geom = "point", size = 3, shape = 21,
                 fill = "white", color = "black") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, color = "black") +
    
    facet_wrap(~ measure, nrow = 1, scales = "free_y") +
    
    geom_text(
      data = ttests_annot,
      aes(x = 1.5, y = y, label = label),
      inherit.aes = FALSE,
      size = 4
    ) +
    
    labs(x = "", y = "Score", title = title) +
    scale_fill_manual(values = c("pre" = "#F7E6A6", "post" = "#E49806")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")
  
  return(p)
}

#===============================================================================
# DUMBBELL PLOTS FOR CORRELATION CHANGES
#===============================================================================

# Robust ordering helper
.vars_in_order <- function(tidy_one_family, prefer_pc = "PC_thin") {
  v <- tidy_one_family %>%
    dplyr::filter(pc == prefer_pc) %>%
    dplyr::pull(variable) %>% unique()
  if (length(v)) return(v)
  tidy_one_family %>% dplyr::pull(variable) %>% unique()
}

# Dumbbell plot for correlation changes
plot_family_dumbbell_dir_facet <- function(tidy, family_label,
                                           col_up = "#FF0000", col_down = "#039EE3") {
  dat_all <- tidy %>%
    dplyr::filter(family == family_label)
  
  var_lvls <- .vars_in_order(dat_all, prefer_pc = "PC_thin")
  if (!length(var_lvls)) var_lvls <- unique(dat_all$variable)
  
  dat <- dat_all %>%
    dplyr::mutate(
      variable = factor(variable, levels = rev(var_lvls)),
      delta_r  = post_r - pre_r,
      dir      = ifelse(delta_r >= 0, "Increase", "Decrease"),
      mid_r    = (pre_r + post_r) / 2
    ) 
  
  p <- ggplot(dat, aes(y = variable)) +
    geom_point(aes(x = pre_r), shape = 21, size = 4, stroke = 0.5,
               fill = "white", color = "black", na.rm = TRUE) +
    geom_point(aes(x = post_r), shape = 16, size = 4,
               color = "black", na.rm = TRUE) +
    geom_segment(aes(x = pre_r, xend = post_r, y = variable, yend = variable, color = dir),
                 linewidth = 1.2,
                 arrow = grid::arrow(length = unit(2.2, "mm"), type = "closed", ends = "last"),
                 na.rm = TRUE) +
    geom_text(aes(x = mid_r, label = sprintf("Δr = %+.2f", delta_r)),
              color = "black", size = 4.5, vjust = -0.6, show.legend = FALSE, na.rm = TRUE) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    coord_cartesian(xlim = c(-1, 1)) +
    facet_wrap(~ pc, nrow = 1) +
    scale_color_manual(values = c(Increase = col_up, Decrease = col_down), name = "Δr change") +
    scale_y_discrete(labels = ~ gsub("_", " ", .x)) +
    labs(title = paste0("Change in PCA-brain association | ",
                        gsub("_", " ", family_label)), 
         subtitle = "(\u25CB pre-training | \u25CF post-training | experimental group)",
         x = "Pearson r (–1 to 1)", y = NULL) +
    theme_test(base_size = 20) +
    theme(panel.grid.major.y = element_blank(),
          legend.position='top', 
          legend.justification='right',
          legend.direction='horizontal')
  
  p
    
  return(p)
  
}



#===============================================================================
# PLOT: delta (post - pre) by group
#===============================================================================
anova_followup_delta_plot <- function(res_delta,
                                      title = "Change score (post − pre)",
                                      subtitle = NULL,
                                      x_label = "",
                                      y_label = "Delta (post − pre)",
                                      filename = NULL,
                                      width = 6,
                                      height = 4) {
  
  df <- res_delta$delta_df
  if (!("group" %in% names(df)) || !("delta" %in% names(df))) {
    stop("res_delta must come from anova_followup_delta() and contain delta_df with columns: group, delta.")
  }
  
  # Pull p-value for subtitle (optional)
  if (is.null(subtitle) && !is.null(res_delta$t_delta)) {
    pval <- res_delta$t_delta$p.value
    subtitle <- paste0("t-test on delta by group: p = ", sprintf("%.3f", pval))
  }
  
  df$group <- factor(df$group)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = delta)) +
    
    # Mean bar
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      fun = mean,
      geom = "bar",
      width = 0.6,
      alpha = 0.8
    ) +
    
    # Mean ± SE
    ggplot2::stat_summary(
      fun.data = ggplot2::mean_se,
      geom = "errorbar",
      width = 0.15,
      linewidth = 0.6
    ) +
    
    # Individual data points
    ggplot2::geom_point(
      ggplot2::aes(fill = group, shape = group),
      size = 2.8,
      alpha = 0.9,
      position = ggplot2::position_jitter(width = 0.08, height = 0),
      color = "black",
      stroke = 0.4
    ) +
    
    # Zero reference line
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
    
    # Manual scales
    ggplot2::scale_fill_manual(values = c(
      "exp" = color$exp,
      "con" = color$con
    )) +
    ggplot2::scale_shape_manual(values = c(
      "exp" = 25,  # filled inverted triangle
      "con" = 21   # filled circle
    )) +
    
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = x_label, y = y_label) +
    
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(),
      legend.position = "none"
    )
  
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = width, height = height, dpi = 300)
  }
  
  p
}



anova_followup_detailed_plot <- function(data, variable,
                                         y_label = NULL,
                                         title = "Interaction",
                                         point_color = "black",
                                         point_shape = 16,
                                         ylim = NULL,
                                         filename = NULL,
                                         width = 12,
                                         height = 4) {
  
  if (is.null(y_label)) y_label <- gsub("_", " ", variable)
  
  df <- data %>%
    dplyr::filter(.data$training %in% c("pre", "post")) %>%
    dplyr::select(subject, group, training, value = dplyr::all_of(variable)) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::mutate(
      training = factor(training, levels = c("pre", "post")),
      group    = factor(group, levels = c("con", "exp"))
    )
  
  # Shared y-limits
  if (is.null(ylim)) ylim <- range(df$value, na.rm = TRUE)
  y_top <- ylim[2] - 0.05 * diff(ylim)   # text position
  y_bar <- ylim[2] - 0.10 * diff(ylim)   # bracket position
  
  exp_group <- dplyr::filter(df, group == "exp")
  con_group <- dplyr::filter(df, group == "con")
  
  # ---- paired t-tests (robust: match pre/post by subject) ----
  get_paired_p <- function(d) {
    wide <- d %>%
      dplyr::select(subject, training, value) %>%
      tidyr::pivot_wider(names_from = training, values_from = value) %>%
      dplyr::filter(!is.na(pre), !is.na(post))
    if (nrow(wide) < 2) return(NA_real_)
    stats::t.test(wide$pre, wide$post, paired = TRUE)$p.value
  }
  
  p_exp <- get_paired_p(exp_group)
  p_con <- get_paired_p(con_group)
  
  lab_exp <- if (is.na(p_exp)) "p = NA" else paste0("p = ", sprintf("%.3f", p_exp))
  lab_con <- if (is.na(p_con)) "p = NA" else paste0("p = ", sprintf("%.3f", p_con))
  
  # Interaction plot
  interaction_plot <- ggplot2::ggplot(df,
                                      ggplot2::aes(x = training, y = value, color = group, group = group)) +
    ggplot2::stat_summary(fun = mean, geom = "line",
                          ggplot2::aes(linetype = group), linewidth = 1) +
    ggplot2::stat_summary(fun = mean, geom = "point", size = 3) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se,
                          geom = "errorbar", width = 0.1) +
    ggplot2::scale_color_manual(values = c("con" = color$con, "exp" = color$exp)) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(title = title, x = "", y = y_label) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "left")
  
  # Helper to add bracket + p label (pre vs post)
  add_p_annot <- function(p, label_text) {
    p +
      ggplot2::geom_segment(x = 1, xend = 2, y = y_bar, yend = y_bar, linewidth = 0.5) +
      ggplot2::geom_segment(x = 1, xend = 1, y = y_bar, yend = y_bar - 0.02 * diff(ylim), linewidth = 0.5) +
      ggplot2::geom_segment(x = 2, xend = 2, y = y_bar, yend = y_bar - 0.02 * diff(ylim), linewidth = 0.5) +
      ggplot2::annotate("text", x = 1.5, y = y_top, label = label_text, size = 4)
  }
  
  # Experimental group plot
  exp_plot <- ggplot2::ggplot(exp_group, ggplot2::aes(x = training, y = value)) +
    ggplot2::geom_point(ggplot2::aes(group = subject),
                        alpha = 0.4,
                        color = color$exp,
                        shape = point_shape,
                        position = ggplot2::position_jitter(width = 0.08)) +
    ggplot2::geom_line(ggplot2::aes(group = subject),
                       alpha = 0.3,
                       color = point_color) +
    ggplot2::stat_summary(fun = mean, geom = "line",
                          color = color$pre$exp, linewidth = 1.2) +
    ggplot2::stat_summary(fun = mean, geom = "point",
                          size = 4, color = color$pre$exp) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(title = "Group: experimental", x = "", y = "") +
    ggplot2::theme_minimal(base_size = 14)
  
  exp_plot <- add_p_annot(exp_plot, lab_exp)
  
  # Control group plot
  con_plot <- ggplot2::ggplot(con_group, ggplot2::aes(x = training, y = value)) +
    ggplot2::geom_point(ggplot2::aes(group = subject),
                        alpha = 0.4,
                        color = color$con,
                        shape = point_shape,
                        position = ggplot2::position_jitter(width = 0.08)) +
    ggplot2::geom_line(ggplot2::aes(group = subject),
                       alpha = 0.3,
                       color = point_color) +
    ggplot2::stat_summary(fun = mean, geom = "line",
                          color = color$con, linewidth = 1.2) +
    ggplot2::stat_summary(fun = mean, geom = "point",
                          size = 4, color = color$con) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(title = "Group: control", x = "", y = "") +
    ggplot2::theme_minimal(base_size = 14)
  
  con_plot <- add_p_annot(con_plot, lab_con)
  
  combined_plot <- ggpubr::ggarrange(interaction_plot, exp_plot, con_plot, nrow = 1)
  
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, combined_plot, width = width, height = height, dpi = 300)
  }
  
  combined_plot
}



#===============================================================================
# PCA × Brain heatmap (PRE): styled to match plot_language_brain_heatmap()
# - Input: psych::corr.test object from calculate_pca_brain_corr()
# - x-axis: brain measures (keeps your brain_order if provided; else uses corr col order)
# - y-axis: PCA components (PC_thin, PC_long) in pca_order
# - Cell label: r + n
# - Star: BH q < q_star
#===============================================================================

plot_pca_brain_heatmap <- function(corr,
                                   title = "Correlations",
                                   subtitle = NULL,
                                   filename = NULL,
                                   q_star = DEFAULT_FDR_THRESHOLD,
                                   brain_order = NULL,
                                   pca_order = c("PC_thin", "PC_long"),
                                   digits_r = 2) {
  
  r     <- corr$r
  p_adj <- corr$p.adj
  n     <- corr$n
  
  # Make n a matrix if it is a single number
  if (length(n) == 1) {
    n <- matrix(n, nrow = nrow(r), ncol = ncol(r), dimnames = dimnames(r))
  }
  
  # Convert matrices -> long dataframe
  df <- as.data.frame(as.table(r), stringsAsFactors = FALSE)
  names(df) <- c("pca", "brain", "r")
  
  df$p_adj <- as.vector(p_adj)
  df$n     <- as.vector(n)
  df$sig   <- !is.na(df$p_adj) & df$p_adj < q_star
  
  # --- enforce PCA order (rows) ---
  pca_levels <- rownames(r)
  # prefer user-specified order if it matches what's present
  pca_levels_use <- pca_order[pca_order %in% pca_levels]
  if (length(pca_levels_use) == 0) pca_levels_use <- pca_levels
  df$pca <- factor(df$pca, levels = pca_levels_use)
  df$pca <- gsub("_", " ", as.character(df$pca))
  df$pca <- factor(df$pca, levels = gsub("_", " ", pca_levels_use))
  
  # --- enforce brain order (columns) ---
  brain_levels <- colnames(r)
  if (!is.null(brain_order)) {
    brain_levels_use <- as.character(brain_order)
    # keep only those actually present
    brain_levels_use <- brain_levels_use[brain_levels_use %in% brain_levels]
    if (length(brain_levels_use) == 0) brain_levels_use <- brain_levels
  } else {
    brain_levels_use <- brain_levels
  }
  
  df$brain <- factor(df$brain, levels = brain_levels_use)
  df$brain <- gsub("_", " ", as.character(df$brain))
  df$brain <- factor(df$brain, levels = gsub("_", " ", brain_levels_use))
  
  # Cell label (r + n)
  df$label <- ifelse(is.na(df$r),
                     paste0("NA\n(n = ", df$n, ")"),
                     paste0(sprintf(paste0("%.", digits_r, "f"), df$r), "\n(n = ", df$n, ")"))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = brain, y = pca, fill = r)) +
    ggplot2::geom_tile(color = "grey90") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.5) +
    ggplot2::geom_text(
      data = subset(df, sig),
      label = "*",
      nudge_y = 0.25,
      size = 6
    ) +
    ggplot2::scale_fill_gradient2(
      limits = c(-1, 1),
      low = color$correlation$low,
      mid = color$correlation$mid,
      high = color$correlation$high,
      midpoint = 0,
      name = "Pearson's r"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey30"),
      legend.title.position = "bottom"
    )
  
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 10, height = 6, dpi = 300)
  }
  
  p
}


# =============================
#===============================================================================
# plots.R
# Dumbbell/arrow plot for PCA-prepost correlation change results
# Expects tidy long df with: family, pc, variable, pre_r, post_r, z_p_BH
#===============================================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

plot_family_dumbbell_dir_facet <- function(
    tidy,
    family_label = NULL,
    order_vec    = NULL,
    title        = NULL,
    subtitle     = NULL,
    filename     = NULL,
    col_up       = color$correlation$high,
    col_down     = color$correlation$low,
    bh_thresh    = DEFAULT_FDR_THRESHOLD,
    facet_families = TRUE,
    star_x       = -0.98,
    star_size    = 8
) {
  
  dat_all <- tidy
  if (!is.null(family_label)) {
    dat_all <- dat_all %>% dplyr::filter(family %in% family_label)
  }
  
  dat_all <- dat_all %>%
    dplyr::mutate(
      delta_r = post_r - pre_r,
      dir     = ifelse(delta_r >= 0, "Increase", "Decrease"),
      mid_r   = (pre_r + post_r) / 2,
      sig     = !is.na(z_p_BH) & z_p_BH < bh_thresh
    )
  
  if (!is.null(order_vec)) {
    dat <- dat_all %>%
      dplyr::mutate(variable = factor(variable, levels = rev(order_vec))) %>%
      dplyr::arrange(family, desc(variable))
  } else {
    dat <- dat_all %>%
      dplyr::group_by(family) %>%
      dplyr::group_modify(function(.x, .k) {
        lv <- .vars_in_order(.x, prefer_pc = "PC_thin")
        if (!length(lv)) lv <- unique(.x$variable)
        .x %>% dplyr::mutate(variable = factor(variable, levels = rev(lv)))
      }) %>%
      dplyr::ungroup()
  }
  
  thm <- if (exists("theme_test")) theme_test(base_size = 16) else theme_minimal(base_size = 16)
  
  p <- ggplot(dat, aes(y = variable)) +
    geom_segment(aes(x = pre_r, xend = post_r, yend = variable, color = dir,
                     linewidth = ifelse(sig, 1.2, 0.7)),
                 arrow = grid::arrow(length = unit(2.2, "mm"), type = "closed", ends = "last"),
                 na.rm = TRUE) +
    geom_point(aes(x = pre_r), shape = 21, size = 3, stroke = 0.5,
               fill = "white", color = "black", na.rm = TRUE) +
    geom_point(aes(x = post_r, shape = dir, color = dir), size = 3, na.rm = TRUE) +
    scale_shape_manual(values = c(Increase = 62, Decrease = 60),
                       breaks = c("Increase","Decrease"), na.translate = FALSE, name = "Δr change") +
    scale_color_manual(values = c(Increase = col_up, Decrease = col_down),
                       breaks = c("Increase","Decrease"), na.translate = FALSE, name = "Δr change") +
    scale_linewidth_identity() +
    geom_text(aes(x = mid_r, label = sprintf("Δr %+.2f", delta_r)),
              color = "black", size = 3.5, vjust = -0.6, show.legend = FALSE, na.rm = TRUE) +
    geom_text(data = dplyr::filter(dat, sig),
              aes(x = star_x, label = "*"),
              inherit.aes = TRUE, color = "black",
              size = star_size, hjust = 1, vjust = 0.7, na.rm = TRUE) +
    geom_text(aes(x = 0.95, label = paste0("(n = ", n, ")")),
              inherit.aes = TRUE, color = "grey30",
              size = 3.5, hjust = 1, vjust = 0.7, show.legend = FALSE, na.rm = TRUE) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    coord_cartesian(xlim = c(-1, 1)) +
    { if (facet_families) facet_grid(rows = vars(family), cols = vars(pc), scales = "free_y")
      else facet_wrap(~ pc, nrow = 1) } +
    scale_y_discrete(labels = ~ gsub("_", " ", .x)) +
    labs(
      title    = title %||% "Change in PCA–brain association",
      subtitle = subtitle,
      x = "Pearson r (–1 to 1)", y = NULL
    ) +
    thm +
    theme(panel.grid.major.y = element_blank(),
          legend.position = "top",
          legend.justification = "right",
          legend.direction = "horizontal",
          plot.title = element_text(face = "bold", size = rel(1.1)),
          plot.subtitle = element_text(size = rel(0.9), colour = "grey30", margin = margin(b = 6)))
  
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 10, height = 6, dpi = 300)
  }
  
}

