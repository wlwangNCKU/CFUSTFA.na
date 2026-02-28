library(cubature)
library(MomTrunc)
library(mvtnorm)

################################################################################
#
#   Filename: table4.R
#   Purpose: Extract parameter estimates and standard errors from fitted
#            CFUSTFA model (q = 4), format results for reporting, and export
#            as a structured table (including μ, B, d, λ, ν)
#   Input data files: data/hcv.results_q4.Rdata;
#                     data/source/hcv.csv;
#                     function/cfustfa.na.se.hcv.R
#   Output data files: results/table4.csv
#   R Version: R-4.4.1
#   Required R packages: cubature, MomTrunc, mvtnorm
#
################################################################################

Data = as.matrix(read.csv(paste(SPATH, "/data/source/hcv.csv", sep = ""), header=T))
Y.na = scale( Data[,-1])
load(paste(SPATH, "/data/hcv.results_q4.Rdata", sep = ""))
source(paste(SPATH, "/function/cfustfa.na.se.hcv.R", sep = ""))
#sderr=cfustna.se(Y.na, q=4, s=4, para.est=fit.ust$para)
#round(sderr$SE, 3)

# Estimate CFUSTFA model parameters (with standard errors)
result <- cfustna.se(Y.na, q=4, s=4, para.est=fit.ust$para)
res_df <- as.data.frame(result)
res_df$param <- rownames(res_df)
rownames(res_df) <- NULL
colnames(res_df) <- c("est", "se", "param")

# Extract parameters: μ, B, d, λ, ν
# μ
mu_df <- res_df[grep("^mu", res_df$param), ]
mu_df$var <- c("ALB", "ALP", "ALT", "AST", "BIL", "CHE", "CHOL", "CREA", "GGT", "PROT")

# B
B_df <- res_df[grep("^B", res_df$param), ]
B_df$idx <- as.numeric(sub("B (\\d+),(\\d+)", "\\1", B_df$param))
B_df$col <- as.numeric(sub("B (\\d+),(\\d+)", "\\2", B_df$param))
B_df$col <- as.numeric(gsub(".*,", "", gsub("\\)", "", B_df$param)))
B_mat <- reshape(B_df[, c("idx", "col", "est", "se")],
                 timevar = "col", idvar = "idx", direction = "wide")
colnames(B_mat)[1] <- "var_idx"
B_mat$var <- mu_df$var

# d
d_df <- res_df[grep("^dd", res_df$param), ]
d_df$var <- mu_df$var

# λ and ν
lambda <- res_df[grep("^la", res_df$param), ]
nu <- res_df[grep("^nu", res_df$param), ]
lambda_row1 <- sprintf("%.3f", lambda$est)
lambda_row2 <- sprintf("(%.3f)", lambda$se)
nu_row1 <- sprintf("%.3f", nu$est)
nu_row2 <- sprintf("(%.3f)", nu$se)

# Format results for reporting
final_df <- data.frame(
  Variable = mu_df$var,
  mu = sprintf("%.3f (%.3f)", mu_df$est, mu_df$se),
  col1_B = sprintf("%.3f (%.3f)", B_mat$est.1, B_mat$se.1),
  col2_B = sprintf("%.3f (%.3f)", B_mat$est.2, B_mat$se.2),
  col3_B = sprintf("%.3f (%.3f)", B_mat$est.3, B_mat$se.3),
  col4_B = sprintf("%.3f (%.3f)", B_mat$est.4, B_mat$se.4),
  d = sprintf("%.3f (%.3f)", d_df$est, d_df$se),
  stringsAsFactors = FALSE
)

# Add lambda (la)
lambda_row <- data.frame(
  Variable = "la",
  mu = lambda_row1[1],
  col1_B = lambda_row1[2],
  col2_B = lambda_row1[3],
  col3_B = lambda_row1[4],
  col4_B = lambda_row1[5],
  d = "",
  stringsAsFactors = FALSE
)

lambda_se_row <- data.frame(
  Variable = "",
  mu = lambda_row2[1],
  col1_B = lambda_row2[2],
  col2_B = lambda_row2[3],
  col3_B = lambda_row2[4],
  col4_B = lambda_row2[5],
  d = "",
  stringsAsFactors = FALSE
)

# Add nu
nu_row <- data.frame(
  Variable = "nu",
  mu = "", col1_B = "", col2_B = "", col3_B = "", col4_B = "",
  d = nu_row1,
  stringsAsFactors = FALSE
)

nu_se_row <- data.frame(
  Variable = "",
  mu = "", col1_B = "", col2_B = "", col3_B = "", col4_B = "",
  d = nu_row2,
  stringsAsFactors = FALSE
)

# Combine all into final table
final_table <- rbind(final_df, lambda_row, lambda_se_row, nu_row, nu_se_row)

#  Export cleaned parameter table to CSV
write.csv(final_table, file = paste(SPATH, "/results/table4.csv", sep = ""), row.names = FALSE, quote = FALSE)
