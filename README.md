# CFUSTFA.na
Supplement: Data and Code for "Handling missing data, skewness, and outliers in medical research: a robust factor analysis approach using the canonical fundamental skew-t distribution" by Wan-Lun Wang, Luis M. Castro and Tsung-I Lin*

# Author responsible for the code #
For questions, comments or remarks about the code please contact responsible author, Wan-Lun Wang (wangwl@gs.ncku.edu.tw), Luis M. Castro (lmcastro@uc.cl) and Dr. Tsung-I Lin (tilin@nchu.edu.tw).

# Configurations #
The code was written/evaluated in R with the following software versions:
R version 4.2.1 (2022-06-23 ucrt) 
Platform: x86_64-w64-mingw32/x64 (64-bit)  
Running under: Windows 10 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.utf8 LC_CTYPE=Chinese (Traditional)_Taiwan.utf8
[3] LC_MONETARY=Chinese (Traditional)_Taiwan.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Traditional)_Taiwan.utf8   

attached base packages:
[1] stats4 stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] mvtnorm_1.1-3     mnormt_2.1.1      relliptical_1.0   rgl_0.110.2       MASS_7.3-58.1  
[6] ks_1.13.5         misc3d_0.9-1      GGally_2.1.2      VIM_6.2.2         ggplot2_3.4.0  
[11] Cairo_1.6-0       ggExtra_0.10.1    gridExtra_2.3     grid_4.2.1        e1071_1.7-13  
[16] moments_0.14.1    MomTrunc_6.0      tmvtnorm_1.5      cubature_2.0.4.5  dplyr_1.0.10  
[21] tidyr_1.2.1       stringr_1.5.0     writexl_1.4.0     readxl_1.4.1      readr_2.1.3  
[26] xtable_1.8-4

loaded via a namespace (and not attached):
[1] tools_4.2.1      utf8_1.2.2       R6_2.5.1         compiler_4.2.1   cli_3.4.1  
[6] lifecycle_1.0.3  magrittr_2.0.3   rlang_1.0.6      purrr_0.3.5      vctrs_0.5.0  
[11] generics_0.1.3   tibble_3.1.8     pillar_1.8.1     glue_1.6.2       fansi_1.0.3  
[16] xfun_0.35        cachem_1.0.6     jsonlite_1.8.3   sessioninfo_1.2.2 knitr_1.41

# Descriptions of the codes # 
Please extract the file "Data_and_Code_CFUSTFA.na-main.zip" to the "current working directory" of the R package.
The getwd() function shall determine an absolute pathname of the "current working directory".

Before running the codes **fig1.R**, **fig2.R**, **fig3.R**, **fig4.R**, **fit_hcvdata.R**, **simulation1.R**, **simulation2.R**, **table1,S1-S3.R**, **table2.R**, **table3.R**, **table4.R** and **tableS4-S6.R**, one needs to install the following R packages:
  
install.packages("mvtnorm") Version 1.1-3  
install.packages("tmvtnorm") Version 1.5  
install.packages("cubature") Version 2.0.4.5  
install.packages("MomTrunc") Version 6.0  
install.packages("moments") Version 0.14.1  
install.packages("gclus") Version 1.3.2  
install.packages("mclust") Version 6.0.0  
install.packages("rgl") Version 0.110.2  
install.packages("misc3d") Version 0.9-1  
install.packages("plot3D") Version 1.4  
install.packages("matrixcalc") Version 1.0-6  
install.packages("devtools") Version 2.4.5  
install.packages("ggplot2") Version 3.4.0

R codes for the implementation of the proposed methodology are provided.

## Subfolder: ./function ##
`./function`
contains the program (function) of 
- (1) **cfustfa.na.R**: main script for calculating the parameter estimates via the AECM algorithm under the CFUSTFA model with missing information;
- (2) **cfustfa.na.se.hcv.R**: main script for computing the standard errors of the parameter estimates under the CFUSTFA model via the observed information matrix, specifically for the Hepatitis C Virus (HCV) dataset;
- (3) **cfustfa.na.se.simu.R** main script for computing the standard errors of the parameter estimates under the CFUSTFA model via the observed information matrix for the simulation dataset;
- (4) **cfustfn.R** main script for generating  data from the CFUSTFA model (rcfustfa), simulating samples from the restricted skew-t distribution (rumst), and computing the CFUST density (dCFUST) to support likelihood-based calculations and simulation studies;
- (5) **FA.na.R**: for calculating the parameter estimates via the EM algorithm under the FA model with missing information; 
- (6) **rSNFA.na.R**: main script for calculating the parameter estimates via the EM algorithm under the rSNFA model with missing information;
- (7) **rSTFA.na.R**: main script for calculating the parameter estimates via the EM algorithm under the rSTFA model with missing information; 
- (8) **tFA.R**: main script for calculating the parameter estimates via the EM algorithm under the tFA model with missing information; and 
- (9) *gener_na.R*: main script for producting missing values.

## Subfolder: ./code ##
`./code`
contains
- (1) **fig1.R**: main script for drawing the 3-dimensional density contours of the rCFUST distribution (source **cfustfn.R**, and then run **fig1.R**);
- (2) **simulation1.R**: main script for conducting Experiment 1 to re-generate simulated datasets and estimating model parameters under different sample sizes and missing rates; 
- (3) **simulation2.R**: main script for conducting Experiment 2 to examine the performance of the CFUSTFA model relative to the FA, tFA, rSNFA and rSTFA approaches across various sample sizes and missing rates;
- (4) **table1,S1-S3.R**: main script for summarizing estimation accuracy using 2 criteria — STD and IMSE — based on the simulation results (load files in the `./simu1/para/*.csv` and `./simu1/sd/*.csv` subsubfolders, and then run **table1,S1-S3.R**);
- (5) **tableS4-S6.R**: main script for performance comparison of 5 models (n=300, 600, 900) based on 100 Monte Carlo trials using AIC, BIC, MSPE, and CT (load files in the `./simu2/*.csv' subsubfolder, and then run **tableS4-S6.R**);
- (6) **fig2.R**: main script for creating RMSE comparison plots across different parameters, sample sizes, and missing data rates (load files in the `./simu1/para/*.csv' subsubfolder first, and then run **fig2.R**);
- (7) **table2.R**: main script for summarizing the range, mean, SD, skewness, kurtosis, and missingness rate of biochemical markers (read **hcv.csv** in the './data/source/' subsubfoloder first, and then run **table2.R**);
- (8) **table3.R**: main script for computing and comparing AIC/BIC values across 5 models and 5 latent dimensions (read **hcv.results_q1.Rdata** to **hcv.results_q5.Rdata** in the `./data/` subfolder first, and then run **table3.R**);
- (9) **table4.R**: main script for reporting ML estimates and standard errors of CFUSTFA model parameters for q = 4 (read **hcv.csv** in the './data/source/' subsubfoloder first, then load **hcv.results_q4.Rdata** in the `./data/` subfoloder, and finally run **cfustfa.na.se.hcv.R**);
- (10) **fig3.R**: main script for creating pairwise scatter plots with histograms and correlation coefficients using data (read **hcv.csv** in the './data/source/' subsubfoloder first, and then run **fig3.R**);
- (11) **fig4.R**: main script for plotting model-wise factor score comparisons between CFUSTFA and 4 benchmark models (load **hcv.results_q4.Rdata** in the `./data/` subfolder, and then run **fig4.R**);
- (12) **fit_hcvdata.R**: main script for model fitting (FA, tFA, rSNFA, rSTFA, CFUSTFA) to the hepatitis C dataset across q = 1–5; saved results are stored in `./data/`.

## Subfolder: ./data ##
`./data`
contains
- (1) **hcv.results_q1.Rdata** to **hcv.results_q5.Rdata**: store the fitted results of the five competing models (FA, tFA, rSNFA, rSTFA, CFUSTFA) under different factor numbers q=1,2,3,4,5. They are used to directly generate the outputs (tables and figures) without re-running the full estimation process.

`./data/source`
subsubfolder contains
- (2) **hcv.csv** contains the hepatitis C dataset used for real-data analysis in Section 6.

## Subfolder: ./results ##
`./results`
contains
- (1) **fig1.png**: 3D density contour and scatter plot illustration for the CFUST distribution;
- (2) **fig2.eps**: (Figure 2) RMSE line plots across different parameters, sample sizes, and missing data rates based on simulation study;
- (3) **fig3.eps**: (Figure 3) pairwise scatter plots, marginal distributions, and Pearson correlations of the biochemical markers from the Hepatitis C dataset;
- (4) **fig4.eps**: (Figure 4) comparison plots of estimated factor scores between CFUSTFA and alternative models across four latent factors;
- (5) **table1.csv**: (Table 1) simulation results for assessing the precision of parameter estimates with r = 30% missing values;
- (6) **table2.csv**: (Table 2) summary statistics of 10 biochemical variables in the Hepatitis C dataset, including range, mean, skewness, kurtosis, and missing percentage;
- (7) **table3.csv**: (Table 3) model comparison results across 5 factor numbers for FA, tFA, rSNFA, rSTFA, and CFUSTFA, evaluated using log-likelihood, AIC, and BIC;
- (8) **table4.csv**: (Table 4) maximum likelihood estimates and associated standard errors of parameters under the best-fitted CFUSTFA model (q = 4);
- (9) **tableS1.csv**: (Table S.1) simulation summary statistics (STD, IMSE) of parameter estimates under different sample sizes with r = 0% missing values;
- (10) **tableS2.csv**: (Table S.2) simulation summary statistics (STD, IMSE) of parameter estimates under different sample sizes with r = 10% missing values;
- (11) **tableS3.csv**: (Table S.3) simulation summary statistics (STD, IMSE) of parameter estimates under different sample sizes with r = 20% missing values;
- (12) **tableS4.csv**: (Table S.4) performance comparison of 5 factor analysis models with n = 300 based on 100 Monte Carlo trials;
- (13) **tableS5.csv**: (Table S.5) performance comparison of 5 factor analysis models with n = 600 based on 100 Monte Carlo trials;
- (14) **tableS6.csv**: (Table S.6) performance comparison of 5 factor analysis models with n = 900 based on 100 Monte Carlo trials;
- (15) `./simu1/` subsubfolder: storing intermediate simulation results, including estimated parameters (`./results/simu1/para/*.csv`) and standard errors (`./results/simu1/sd/*.csv`) for various settings of sample sizes and missing rates in Experiment 1;
- (16) `./simu2/` subsubfolder: storing intermediate simulation outputs, specifically various model selection criteria and performance metrics such as log-likelihood, AIC, BIC, MSPE, and CT for various settings of sample sizes and missing rates in Experiment 2.

## Additional Remark ##
- Note (1): One can directly run each "source(.)" described in **master.R** file in the seperate R session to obtain the results.
- Note (2): Because the estimation procedures in **fit_hcvdata.R** involve fitting five different models across multiple factor dimensions, the computations are time-consuming. Therefore, we have saved the fitted results in **hcv.results_q1.Rdata** to **hcv.results_q5.Rdata** under the `./data/` subfolder.
- Note (3): Since **simulation1.R** takes a long time to run, to reproduce Tables 1 and S1–S3 and Figure 2 in Section 5,  we record these intermediately numerical results in `./results/simu1/para/` and `./results/simu1/sd/` subsubfolders so that one can use the R codes *table1,S1-S3.R** script in subfolder `./code/` to obtain the final results based on the pre-saved .csv files.
- Note (4): Since **simulation2.R** takes a long time to run, to reproduce Tables S4-S6 in the Supplementary Material, we record these intermediately numerical outputs in the `./results/simu2/` subsubfolder so that one can use the R codes **tableS4-S6.R** to obtain the final results based on the pre-saved .csv files.
