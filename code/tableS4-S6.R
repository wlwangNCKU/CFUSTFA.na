library(dplyr)

################################################################################
#
#   Filename: tableS4-6.R
#   Purpose: Compare performance of 5 factor analysis models across different
#             sample sizes (n = 300, 600, 900) based on 100 Monte Carlo trials.
#             Performance is evaluated using AIC, BIC, MSPE, and CT.
#   Input data files: results/simu2/*.csv
#   Output data files:
#             - results/tableS4.csv: Results for n = 300
#             - results/tableS5.csv: Results for n = 600
#             - results/tableS6.csv: Results for n = 900
#   R Version: R-4.4.1
#   Required R packages: dplyr
#
################################################################################

RESULT_PATH <- paste0(SPATH, '/results/simu2/')

get_freq <- function(data_matrix) {
  best_model_idx <- apply(data_matrix, 1, which.min)
  freq_counts <- table(factor(best_model_idx, levels = 1:5))
  return(as.vector(freq_counts))
}

process_simulation_block <- function(n, mr, path) {
  models <- c("n", "t", "rsn", "rst", "cfust")
  data_list <- list()

  for (m in models) {
    file_name <- paste0(path, "(", n, "_", mr, "_", m, ").csv")
    if (file.exists(file_name)) {
      data_list[[m]] <- read.csv(file_name, header = TRUE)
    } else {
      warning(paste("file not found:", file_name))
      return(NULL)
    }
  }

  x1 <- data_list[["n"]]; x2 <- data_list[["t"]]; x3 <- data_list[["rsn"]]
  x4 <- data_list[["rst"]]; x5 <- data_list[["cfust"]]

  metrics <- list(
    AIC  = cbind(x1[,3], x2[,3], x3[,3], x4[,3], x5[,3]),
    BIC  = cbind(x1[,4], x2[,4], x3[,4], x4[,4], x5[,4]),
    MSPE = cbind(x1[,5], x2[,5], x3[,5], x4[,5], x5[,5]),
    CT   = cbind(x1[,6], x2[,6], x3[,6], x4[,6], x5[,6])
  )

  out_matrix <- matrix(NA, 12, 5)
  rownames(out_matrix) <- c('AIC.mean','AIC.sd','AIC.freq',
                            'BIC.mean','BIC.sd','BIC.freq',
                            'MSPE.mean','MSPE.sd','MSPE.freq',
                            'CT.mean','CT.sd','CT.freq')
  colnames(out_matrix) <- c("FA", "tFA", "rSNFA", "rSTFA", "CFUSTFA")

  #AIC
  out_matrix[1,] <- round(colMeans(metrics$AIC), 2)
  out_matrix[2,] <- round(apply(metrics$AIC, 2, sd), 2)
  out_matrix[3,] <- get_freq(metrics$AIC)

  #BIC
  out_matrix[4,] <- round(colMeans(metrics$BIC), 2)
  out_matrix[5,] <- round(apply(metrics$BIC, 2, sd), 2)
  out_matrix[6,] <- get_freq(metrics$BIC)

  #MSPE
  out_matrix[7,] <- round(colMeans(metrics$MSPE), 3)
  out_matrix[8,] <- round(apply(metrics$MSPE, 2, sd), 3)
  out_matrix[9,] <- get_freq(metrics$MSPE)

  #CT
  out_matrix[10,] <- round(colMeans(metrics$CT), 3)
  out_matrix[11,] <- round(apply(metrics$CT, 2, sd), 3)
  out_matrix[12,] <- get_freq(metrics$CT)

  res_df <- data.frame(
    n = c(as.character(n), rep("", 11)),
    Missing_rate = c(as.character(mr), rep("", 11)),
    Criterion = rownames(out_matrix),
    out_matrix,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  return(res_df)
}

n_list <- c(300, 600, 900)
mr_list <- c(0.1, 0.2, 0.3)
table_names <- c("300" = "tableS4", "600" = "tableS5", "900" = "tableS6")

for (n in n_list) {
  current_n_tables <- list()
  for (mr in mr_list) {
    block <- process_simulation_block(n, mr, RESULT_PATH)
    if (!is.null(block)) {
      current_n_tables[[as.character(mr)]] <- block
    }
  }

  if (length(current_n_tables) > 0) {
    final_table_n <- do.call(rbind, current_n_tables)
    target_name   <- table_names[as.character(n)]
    output_filename <- paste0(SPATH, "/results/", target_name, ".csv")
    write.csv(final_table_n, file = output_filename, row.names = FALSE)
  }
}
