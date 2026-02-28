library(dplyr)

################################################################################
#
#   Filename: table1,S1-S3.R
#   Purpose: Summarize simulation results by merging parameter estimates and
#             standard errors. This script calculates the sample standard
#             deviation (STD) and the average estimated standard error (IMSE)
#             to evaluate estimator accuracy.
#   Input data files: results/simu1/para/*.csv (estimated parameters)
#                     results/simu1/sd/*.csv   (estimated standard errors)
#   Output data files: results/table1.csv,tableS1.csv,tableS2.csv,tableS3.csv
#   R Version: R-4.4.1
#   Required R packages: dplyr
#
################################################################################

para_path <- paste(SPATH, "/results/simu1/para", sep="")
sd_path   <- paste(SPATH, "/results/simu1/sd", sep="")

mu <- c(4, 2, 6, 8, 4)
B <- t(matrix(c(0.13, 0.13, 0.97, 0.94, 0.54,
                0.84, 0.84, 0.27, 0.71, 0.77), nrow = 2, byrow = TRUE))
D <- diag(c(9, 10, 2, 3, 7))
la <- c(7, 3, 1, 4)
nu <- 4

true.para <- c(mu, as.vector(B), diag(D), la, nu)
names(true.para) <- c(paste0("mu", 1:5),
                      paste0("b", 1:5, "1"),
                      paste0("b", 1:5, "2"),
                      paste0("d", 1:5, 1:5),
                      paste0("la", 1:2, "1"),
                      paste0("la", 1:2, "2"),
                      "nu")

para_files <- list.files(para_path, full.names = TRUE)
sd_files   <- list.files(sd_path,   full.names = TRUE)
para_files <- para_files[grepl("para", para_files)]
sd_files   <- sd_files[grepl("sd", sd_files)]

para_names <- trimws(gsub("para$", "", tools::file_path_sans_ext(basename(para_files))))
sd_names   <- trimws(gsub("sd$", "", tools::file_path_sans_ext(basename(sd_files))))
common_names <- intersect(para_names, sd_names)

result.list <- list()

for (name in common_names) {
  para_file <- para_files[para_names == name]
  sd_file   <- sd_files[sd_names == name]

  est <- read.csv(para_file, header = TRUE) %>% as.matrix()
  se  <- read.csv(sd_file, header = TRUE) %>% as.matrix()

  if (ncol(est) != length(true.para)) {
    warning(paste("Column mismatch, skipping:", name))
    next
  }

  std  <- apply(est, 2, sd, na.rm = TRUE)
  imse <- colMeans(se, na.rm = TRUE)

  result <- rbind(std  = round(std, 3),
                  imse = round(imse, 3))

  colnames(result) <- names(true.para)
  result.list[[name]] <- result
}

all_settings <- names(result.list)

n_all   <- as.numeric(gsub(".*\\(([^_]+)_.*", "\\1", all_settings))
r_all <- as.numeric(gsub(".*_([^\\)]+)\\).*", "\\1", all_settings))

unique_rs <- sort(unique(r_all))

table_mapping <- c(
  "0"   = "tableS1",
  "0.1" = "tableS2",
  "0.2" = "tableS3",
  "0.3" = "table1"
)

for (current_r in unique_rs) {

  idx <- which(r_all == current_r)
  sub_settings <- all_settings[idx]
  sub_n_values <- n_all[idx]

  order_idx <- order(sub_n_values)
  sorted_sub_names <- sub_settings[order_idx]

  r_blocks <- lapply(sorted_sub_names, function(name) {
    n_val <- gsub(".*\\(([^_]+)_.*", "\\1", name)
    temp_df <- as.data.frame(result.list[[name]])
    formatted_df <- data.frame(
      n = c(n_val, ""),
      Measure = c("STD", "IMSE"),
      temp_df,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    return(formatted_df)
  })

  final_r_table <- do.call(rbind, r_blocks)
  target_name <- table_mapping[as.character(current_r)]

  if (is.na(target_name)) {
    target_name <- paste0("table1_r_", current_r)
  }

  file_output <- paste0(SPATH, "/results/", target_name, ".csv")
  write.csv(final_r_table, file = file_output, row.names = FALSE)
}
