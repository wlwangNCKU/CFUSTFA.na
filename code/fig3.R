library(GGally)
library(VIM)
library(ggplot2)
library(Cairo)

################################################################################
#
#   Filename: fig3.R
#   Purpose: Visualize pairwise relationships among liver function variables
#            after KNN imputation, highlighting originally missing values
#   Input data files: data/source/hcv.csv
#   Output data files: results/fig3.eps
#   R Version: R-4.4.1
#   Required R packages: GGally, VIM, ggplot2, Cairo
#
################################################################################

df=read.csv(paste(SPATH, "/data/source/hcv.csv", sep=''))
selected_vars <- df[, c("ALB", "ALP", "ALT", "AST", "BIL", "CHE", "CHOL", "CREA", "GGT", "PROT")]
missing_idx <- is.na(selected_vars)
df_knn <- kNN(selected_vars, k = 5)
df_knn <- df_knn[, 1:ncol(selected_vars)]
missing_rows <- which(rowSums(missing_idx) > 0)
custom_diag <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, fill = "lightblue", color = "lightblue", alpha = 0.7) +
    geom_density(color = "steelblue", size = 0.4) +
    theme_classic() +
    theme(
      #panel.grid.major = element_line(color = "gray80"),
      #panel.grid.minor = element_line(color = "gray90"),
      #panel.background = element_rect(fill = "gray95", color = NA),  # **?霈滓??*
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}
custom_lower <- function(data, mapping, ...) {
  var_x <- rlang::as_label(mapping$x)
  var_y <- rlang::as_label(mapping$y)
  p <- ggplot(data, mapping) +
    geom_point(color = "lightblue", alpha = 0.8, size = 1.2, shape = 20) +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.4) +  # 靽格迤 size 霅血?
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  if (length(missing_rows) > 0 && var_x %in% colnames(missing_idx) && var_y %in% colnames(missing_idx)) {
    miss_idx_xy <- which(missing_idx[, var_x] | missing_idx[, var_y])
    if (length(miss_idx_xy) > 0) {
      knn_points_xy <- df_knn[miss_idx_xy, , drop = FALSE]
      p <- p + geom_point(data = knn_points_xy, mapping = mapping,
                          color = "darkorange", alpha = 1, shape = 4, size = 1.5)
    }
  }

  return(p)
}
custom_upper <- function(data, mapping, ...) {
  ggally_cor(data, mapping, size = 3, color = "steelblue") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}
p <- ggpairs(
  df_knn,
  lower = list(continuous = custom_lower),
  diag = list(continuous = custom_diag),    upper = list(continuous = custom_upper)
) +
  theme(
    text = element_text(size = 14),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(t = 10, r = 10, b = 80, l = 10)
  )

ggsave(paste(SPATH,'/results/fig3.eps', sep=''), plot = p,
       width = 7, height = 7, dpi = 300,
       device = cairo_ps, fallback_resolution = 300)
