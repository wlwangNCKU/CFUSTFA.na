rm(list = ls())
SPATH = paste(getwd(),"/CFUSTFA.na",sep="")

# Re-produce Figure 1
source(paste(SPATH, '/code/fig1.R', sep=''))

# Run Simulation 1
source(paste(SPATH, '/code/simulation1.R', sep=''))

# Run Simulation 2
source(paste(SPATH, '/code/simulation2.R', sep=''))

# fit hcv.RData
source(paste(SPATH, '/code/fit_hcvdata.R', sep=''))

# Re-produce Figure 2
source(paste(SPATH, '/code/fig2.R', sep=''))

# Re-produce Figure 3
source(paste(SPATH, '/code/fig3.R', sep=''))

# Re-produce Figure 4
source(paste(SPATH, '/code/fig4.R', sep=''))

# Re-produce Table 1 and Tables S1-S3
source(paste(SPATH, '/code/table1,S1-S3.R', sep=''))

# Re-produce Table 2
source(paste(SPATH, '/code/table2.R', sep=''))

# Re-produce Table 3
source(paste(SPATH, '/code/table3.R', sep=''))

# Re-produce Table 4
source(paste(SPATH, '/code/table4.R', sep=''))

# Re-produce Table S4-S6
source(paste(SPATH, '/code/tableS4-S6.R', sep=''))
