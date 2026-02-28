library(ggplot2)
library(ggExtra)
library(gridExtra)
library(grid)
library(e1071)

################################################################################
#
#   Filename: fig4.R
#   Purpose: Visualize and compare factor scores between CFUSTFA and four
#            alternative models (FA, tFA, rSNFA, rSTFA), including marginal histograms
#            and annotated skewness values
#   Input data files: data/hcv.results_q4.Rdata
#   Output data files: results/fig4.eps
#   R Version: R-4.4.1
#   Required R packages: ggplot2, ggExtra, gridExtra, grid, e1071
#
################################################################################

load(paste(SPATH, "/data/hcv.results_q4.Rdata", sep=''))

cfustfa <- as.data.frame(scale(fit.ust$factor.score))
colnames(cfustfa) <- paste0("Factor", 1:4)
nfa  <- as.data.frame(t(fit.n$U))
tfa  <- as.data.frame(t(fit.t$U))
snfa <- as.data.frame(fit.sn$score)
stfa <- as.data.frame(fit.st$score)
colnames(nfa) <- colnames(tfa) <- colnames(snfa) <- colnames(stfa) <- paste0("Factor", 1:4)

model_list <- list(
  "FA"     = nfa,
  "tFA"    = tfa,
  "rSNFA"  = snfa,
  "rSTFA"  = stfa
)

factor_colors <- c(
  "Factor1" = "#e41a1c",
  "Factor2" = "#377eb8",
  "Factor3" = "tan4",
  "Factor4" = "#984ea3"
)


plot_compare <- function(cfust_df, model_df, factor_name, model_name, fill_x, model_skew, cfust_skew) {
  df <- data.frame(
    CFUSTFA = cfust_df[[factor_name]],
    Model = model_df[[factor_name]]
  )
  df <- na.omit(df)

  skew_text <- sprintf("skewness\n%s = %.2f\nCFUSTFA = %.2f", model_name, model_skew, cfust_skew)

  p <- ggplot(df, aes(x = CFUSTFA, y = Model)) +
    geom_point(color = "orange", alpha = 0.6, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    annotate("text", x = 4.8, y = -4.8
             , label = skew_text, size = 1, hjust = 1, vjust = 0, fontface = "italic") +
    coord_fixed(ratio = 1, xlim = c(-5, 5), ylim = c(-5, 5), clip = "on") +
    labs(x = "CFUSTFA", y = model_name) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      aspect.ratio = 1,
      plot.margin = unit(c(0.5, 0.2, 0.5, 0.2), "mm"),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 8)
    )

  ggMarginal(p,
             type = "histogram",
             margins = "both",
             size=5,
             xparams = list(fill = fill_x, alpha = 0.6, bins = 15, color = NA),
             yparams = list(fill = "#00BF7D", alpha = 0.5, bins = 15, color = NA))
}

plots_by_factor <- list()
for (i in 1:4) {
  factor_name <- paste0("Factor", i)
  fill_x <- factor_colors[[factor_name]]

  for (model_name_raw in names(model_list)) {
    model_name <- gsub("[^A-Za-z0-9]", "", model_name_raw)
    model_df <- model_list[[model_name_raw]]

    model_skew <- skewness(model_df[[factor_name]], na.rm = TRUE)
    cfust_skew <- skewness(cfustfa[[factor_name]], na.rm = TRUE)

    plots_by_factor[[paste(factor_name, model_name, sep = "_")]] <-
      arrangeGrob(
        plot_compare(cfustfa, model_df, factor_name, model_name, fill_x, model_skew, cfust_skew),
        ncol = 1
      )

  }
}

factor_titles <- lapply(1:4, function(i) {
  textGrob(paste("Factor", i), gp = gpar(fontsize = 12, fontface = "bold"))
})

plot_matrix <- matrix(plots_by_factor, nrow = 4, byrow = TRUE)

rows_with_titles <- lapply(1:4, function(i) {
  arrangeGrob(
    factor_titles[[i]],
    arrangeGrob(
      grobs = plot_matrix[i, ],
      ncol = 4,
      widths = unit(rep(1, 4), "null")
    ),
    ncol = 1,
    heights = c(0.8, 6)
  )
})

final_plot <- arrangeGrob(grobs = rows_with_titles, ncol = 1)

ggsave(paste(SPATH, "/results/fig4.eps", sep=''),
       plot = final_plot,
       device = cairo_ps,
       width = 7, height = 7, dpi = 300,
       fallback_resolution = 300)
