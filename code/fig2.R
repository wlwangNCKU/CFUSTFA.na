library(dplyr)
library(stringr)
library(ggplot2)

################################################################################
#
#   Filename: fig2.R
#   Purpose: Compute and visualize RMSE of estimated parameters across different
#            sample sizes and missing rates in simulation experiments
#   Input data files: results/simu1/para/*.csv (parameter estimates)
#   Output data files: results/fig2.eps
#   R Version: R-4.4.1
#   Required R packages: dplyr, stringr, ggplot2
#
################################################################################

# === Step 1: Define true parameter values ===
true.para <- c(4, 2, 6, 8, 4, 0.13, 0.13, 0.97, 0.94, 0.54, 0.84, 0.84,
               0.27, 0.71, 0.77, 9, 10, 2, 3, 7, 7, 3, 1, 4, 4)

# Corresponding parameter names (used for labeling)
parameters <- c("mu1", "mu2","mu3","mu4","mu5",
                "B11", "B21", "B31", "B41", "B51",
                "B12", "B22", "B32", "B42", "B52",
                "d11", "d22", "d33", "d44", "d55",
                "la1", "la2", "la3", "la4", "nu")

# Define the folder containing estimated parameter CSV files

path <- paste(SPATH, "/results/simu1/para",sep="")
files <- list.files(path = path, pattern = "para", full.names = TRUE)

# === Step 2: Compute RMSE from simulation results ===
rmse_list <- lapply(files, function(file) {
  dat <- read.csv(file, header = FALSE)

  # Compute RMSE between estimates and true values
  rmse <- sqrt(colMeans((t(t(dat) - true.para))^2))

  # Extract sample size (N) and missing rate (dropout) from filename
  fname <- basename(file)
  fname_info <- str_match(fname, "\\((\\d+)_?(\\d*\\.?\\d*)\\)para")
  N <- as.numeric(fname_info[2])
  dropout <- ifelse(fname_info[3] == "", "0", fname_info[3])

  cat("read:", fname, "→ N =", N, "dropout =", dropout, "\n")

  # Return one data frame per file
  data.frame(
    parameter = parameters,
    RMSE = rmse,
    N = N,
    dropout = dropout
  )
})

# Combine all RMSE results into one data frame
rmse_df <- bind_rows(rmse_list)

# === Step 3: Format parameter names as math expressions for plotting ===
rmse_df$parameter <- factor(rmse_df$parameter,
                            levels = parameters,
                            labels = c(
                              expression(bold(mu)[1]), expression(bold(mu)[2]), expression(bold(mu)[3]),
                              expression(bold(mu)[4]), expression(bold(mu)[5]),

                              expression(b[11]), expression(b[21]), expression(b[31]),
                              expression(b[41]), expression(b[51]),

                              expression(b[12]), expression(b[22]), expression(b[32]),
                              expression(b[42]), expression(b[52]),

                              expression(d[11]), expression(d[22]), expression(d[33]),
                              expression(d[44]), expression(d[55]),

                              expression(Lambda[11]), expression(Lambda[21]),
                              expression(Lambda[12]), expression(Lambda[22]),

                              expression(nu)
                            ))

# === Step 4: Create RMSE plot ===
p <- ggplot(rmse_df, aes(x = N, y = RMSE, color = dropout, shape = dropout, linetype = dropout)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  facet_wrap(~ parameter, scales = "free", ncol = 5, labeller = label_parsed) +
  scale_x_continuous(breaks = c(150, 300, 600, 900, 1500)) +
  scale_color_manual(values = c("0" = "darkgrey", "0.1" = "#e41a1c", "0.2" = "#377eb8", "0.3" = "#4daf4a")) +
  scale_shape_manual(values = c("0" = 1, "0.1" = 16, "0.2" = 17, "0.3" = 15)) +
  scale_linetype_manual(values = c("0" = "solid", "0.1" = "solid", "0.2" = "dashed", "0.3" = "dotted")) +
  labs(
    x = "Sample Size (n)",
    y = "RMSE",
    color = "Missing Rate",
    shape = "Missing Rate",
    linetype = "Missing Rate"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal",
    aspect.ratio = 1
  )

# === Step 5: Save plot as EPS file ===
ggsave(paste(SPATH, "/results/fig2.eps",sep=""), p,
       width = 8, height = 8, device = cairo_ps)
