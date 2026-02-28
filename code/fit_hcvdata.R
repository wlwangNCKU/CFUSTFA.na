library(mvtnorm)
library(moments)
library(MomTrunc)
library(VIM)

################################################################################
#
#   Filename: fit_hcvdata.R
#   Purpose: Fit multiple factor analysis models (FA, tFA, rSNFA, rSTFA, CFUSTFA)
#            with missing data handling using KNN imputation and AECM/EM algorithms;
#            Save results for different numbers of latent factors (q = 1 to 5)
#   Input data files: data/source/hcv.csv;
#                     function/function/gener_na.r;
#                     function/rSTFA.na.r;
#                     function/cfustfa.na.r;
#                     function/FA.na.r;
#                     function/rSNFA.na.r;
#                     function/tFA.na.r
#   Output data files: data/hcv.results_q1.Rdata
#                      data/hcv.results_q2.Rdata
#                      data/hcv.results_q3.Rdata
#                      data/hcv.results_q4.Rdata
#                      data/hcv.results_q5.Rdata
#   R Version: R-4.4.1
#   Required R packages: mvtnorm, moments, MomTrunc, VIM
#
################################################################################

source(paste(SPATH, '/function/gener_na.r', sep=''))
source(paste(SPATH, '/function/rSTFA.na.r', sep=''))
source(paste(SPATH, '/function/cfustfa.na.r', sep=''))
source(paste(SPATH, '/function/FA.na.r', sep=''))
source(paste(SPATH, '/function/rSNFA.na.r', sep=''))
source(paste(SPATH, '/function/tFA.na.r', sep=''))

Data = as.matrix(read.csv(paste(SPATH, "/data/source/hcv.csv", sep = ""), header=T))
Y.raw = Data[,-1]
p = ncol(Y.raw)
Y.knn = kNN(Y.raw, k = 5)[,1:p]
Y.std = scale(Y.knn)
mu.int = colMeans(Y.std)
V = cov(Y.std)
eig.V = eigen(V)
Y.na = scale(Y.raw)

# fit the model
test=function(q)
{
  B.int = eig.V$ve[,1:q] %*% diag(sqrt(eig.V$va[1:q]), q)
  D.int = diag(diag(V - B.int %*% t(B.int)))
  FF = factanal(Y.knn, q, scores='regression')$scores
  la.int =skewness(FF)
  La.int=diag(la.int,q)
  init.para = list(mu = mu.int, B = B.int, D = D.int, la = la.int, La = La.int, nu = 5)
  fit.st=rSTFA.na.EM(Y.na, q=q, init.para, zero.mu = F, zero.la = F, tol = 1e-6, max.iter=2000, per=100)
  fit.ust=CFUSTFANA.AECM(Y.na, q=q, s=q, init.para, zero.mu = F, zero.la = F, uST=T, tol = 0.0001, max.iter = 10, per = 10)
  fit.n=FA.na.EM(Y.na, q=q, init.para, eta=0.005, tol = 1e-6, max.iter=2000, per=100)
  fit.sn=rSNFA.na.EM(Y.na, q=q, init.para, zero.mu = F, zero.la = F, tol = 1e-6, max.iter=2000, per=100)
  fit.t=tFA.na.EM(Y.na, q=q, init.para, eta=0.005, tol = 1e-6, max.iter=2000, per=100)
  save(fit.st, fit.ust, fit.n, fit.sn, fit.t, file = paste0(SPATH,"/data/hcv.results_q", q, ".Rdata"))
}
test(q=1)
test(q=2)
test(q=3)
test(q=4)
test(q=5)
