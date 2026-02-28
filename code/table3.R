library(dplyr)

################################################################################
#
#   Filename: table3.R
#   Purpose: Extract and summarize model selection criteria (log-likelihood,
#            number of parameters, AIC, BIC) for five models across q = 1 to 5
#   Input data files: data/results_q1.Rdata
#                     data/results_q2.Rdata
#                     data/results_q3.Rdata
#                     data/results_q4.Rdata
#                     data/results_q5.Rdata
#   Output data files: results/table3.csv
#   R Version: R-4.4.1
#   Required R packages: dplyr
#
################################################################################

models <- list(
  FA = "fit.n",
  tFA = "fit.t",
  rSNFA = "fit.sn",
  rSTFA = "fit.st",
  CFUSTFA = "fit.ust"
)

results <- data.frame()
for (q in 1:5) {
  load(file.path(SPATH, paste0("data/hcv.results_q", q, ".Rdata")))
  for (model_name in names(models)) {
    fit <- get(models[[model_name]])
    mds <- fit$mds

    results <- rbind(results, data.frame(
      Model = model_name,
      q = q,
      logLik = round(mds["logli"], 3),
      p = mds["no.para"],
      AIC = round(mds["AIC"], 3),
      BIC = round(mds["BIC"], 3)
    ))
  }
}

results$Model <- factor(results$Model, levels = c("FA", "tFA", "rSNFA", "rSTFA", "CFUSTFA"))
results <- results %>% arrange(Model, q)
rownames(results) <- NULL

write.csv(results, file = paste(SPATH, "/results/table3.csv",sep=""), row.names = FALSE)
