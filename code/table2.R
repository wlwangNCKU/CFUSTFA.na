library(readr)
library(dplyr)
library(moments)
library(xtable)

################################################################################
#
#   Filename: table2.R
#   Purpose: Generate summary statistics (range, mean, SD, skewness, kurtosis,
#            missing rate) for each variable in the dataset and export as CSV
#   Input data files: data/source/hcv.csv
#   Output data files: results/table2.csv
#   R Version: R-4.4.1
#   Required R packages: readr, dplyr, moments, xtable
#
################################################################################

data <- read_csv(paste(SPATH, "/data/source/hcv.csv",sep=""))
data <- data[,-1]

summary_table <- data.frame(
  Variable = names(data),
  Range = sapply(data, function(x) if (is.numeric(x)) {
    paste0(round(min(x, na.rm = TRUE), 2), " - ", round(max(x, na.rm = TRUE), 2))
  } else { NA }),
  Mean = sapply(data, function(x) if (is.numeric(x)) round(mean(x, na.rm = TRUE), 2) else NA),
  SD = sapply(data, function(x) if (is.numeric(x)) round(sd(x, na.rm = TRUE), 3) else NA),
  Skewness = sapply(data, function(x) if (is.numeric(x)) round(skewness(x, na.rm = TRUE), 4) else NA),
  Kurtosis = sapply(data, function(x) if (is.numeric(x)) round(kurtosis(x, na.rm = TRUE), 4) else NA),
  Missing = sapply(data, function(x) round(mean(is.na(x)) * 100, 2))
)

summary_table <- summary_table %>%
  filter(!is.na(Mean))

write.csv(summary_table, file = paste(SPATH, "/results/table2.csv",sep=""), row.names = FALSE)
