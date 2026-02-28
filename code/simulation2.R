library(mvtnorm)
library(MomTrunc)
library(tmvtnorm)
library(moments)
library(relliptical)
library(VIM)

################################################################################
#
#   Filename: simulation2.R
#   Purpose: Conduct simulation study for  examining the performance of the CFUSTFA model
#   relative to the FA, tFA, rSNFA and rSTFA approaches across various sample sizes and missing rates
#   Input data files: function/cfustfn.R;
#                     function/cfustfa.na.R;  function/FA.na.R; function/tFA.na.R;
#                     function/rSNFA.na.R; function/rSTFA.na.R; function/cfustfa.na.R
#   Output data files: results/simu2/<M>/(n_dropout_model).csv
#   R Version: R-4.4.1
#   Required R packages: moments, mvtnorm, MomTrunc, relliptical, tmvtnorm, VIM
#
################################################################################

source(paste0(SPATH, "/function/gener_na.R"))
source(paste0(SPATH, "/function/cfustfn.R"))
source(paste0(SPATH, '/function/FA.na.R'))
source(paste0(SPATH, '/function/tFA.na.R'))
source(paste0(SPATH, '/function/rSNFA.na.R'))
source(paste0(SPATH, '/function/rSTFA.na.R'))
source(paste0(SPATH, "/function/cfustfa.na.R"))

simu2=function(n)
{
  p = 5; q = 2
  mu=c(1.47,-0.04,0.96,1.01,-0.14)
  B = matrix(c(0.05,0.87,0.12,0.69,0.01,0.46,0.07,0.02,0.59,0.66),p,q)
  nu=4
  dd = runif(p)
  U1=rchisq(n, df=1)
  U2=rchisq(n, df=1)
  U=cbind(U1, U2)
  E=t(mu+t(rmvt(n,  sigma = diag(dd), df=nu)))
  Y=t(B%*%t(U))+E
  FF = factanal(Y, q, scores='regression')$scores
  la =skewness(FF)
  list(data=Y, mu=mu, B=B, D=diag(dd), la=la, nu=nu)
}

MC=function(n, q=2, na.rate, M)
{
  fit.n=fit.t=fit.rsn=fit.rst=fit.ust=matrix(NA, nrow=M, ncol=6)
  for(j in 1:M) {
  cat('**************************************', '\n')
  cat('Replication=', j, '\n')
  cat('**************************************', '\n')
  simu=simu2(n=n)
  Y=simu$data
  Y.na=gener.na(Y, na.rate)
  na.posi = is.na(Y.na)
  na.ind = which(as.vector(t(na.posi))==T)
  y.true = as.vector(t(Y))[na.ind]
  p = ncol(Y.na)
  Y.knn = kNN(Y.na, k = 5)[,1:p]
  mu.int = simu$mu; B.int = simu$B; D.int = simu$D; la.int = simu$la; La.int = diag(simu$la, q); nu.int = simu$nu
  init.para = list(mu = mu.int, B = B.int, D = D.int, la = la.int, La = La.int, nu = nu.int)
## Fit FA
   fit.n=FA.na.EM(Y.na, q, init.para, eta=0.005, tol = 1e-6, max.iter=200, per=100)
   mspe.n=mean((y.true-fit.n$Y.mis)^2)
   write.table(t(c(fit.n$mds, mspe.n, fit.n$cpu.time)), file=paste0(SPATH, '/results/simu2/','(',n,'_',na.rate,'_n).csv'), sep=",",append=T, row.names=F, col.names=F, na = "NA")
## Fit tFA
   fit.t=tFA.na.EM(Y.na, q, init.para, eta=0.005, tol = 1e-6, max.iter=200, per=100)
   mspe.t=mean((y.true-fit.t$Y.mis)^2)
   write.table(t(c(fit.t$mds, mspe.t, fit.t$cpu.time)), file=paste0(SPATH, '/results/simu2/','(',n,'_',na.rate,'_t).csv'), sep=",",append=T, row.names=F, col.names=F, na = "NA")
## Fit rSN
   fit.rsn=rSNFA.na.EM(Y.na, q, init.para, zero.mu = F, zero.la = F, tol = 1e-6, max.iter=200, per=100)
   mspe.rsn=mean((y.true-fit.rsn$Y.mis)^2)
   write.table(t(c(fit.rsn$mds, mspe.rsn, fit.rsn$cpu.time)), file=paste0(SPATH, '/results/simu2/','(',n,'_',na.rate,'_rsn).csv'), sep=",",append=T, row.names=F, col.names=F, na = "NA")
## Fit rST
  fit.rst=rSTFA.na.EM(Y.na, q, init.para, zero.mu = F, zero.la = F, tol = 1e-6, max.iter=200, per=100)
  mspe.rst=mean((y.true-fit.rst$Y.mis)^2)
  write.table(t(c(fit.rst$mds, mspe.rst, fit.rst$cpu.time)), file=paste0(SPATH, '/results/simu2/','(',n,'_',na.rate,'_rst).csv'), sep=",",append=T, row.names=F, col.names=F, na = "NA")
## Fit CFUST
  fit.cfust=CFUSTFANA.AECM(Y.na, q, s=q, init.para, zero.mu = F, zero.la = F, uST=T, tol = 1e-6, max.iter = 200, per = 10)
  mspe.cfust=mean((y.true-fit.cfust$Y.mis)^2)
  write.table(t(c(fit.cfust$mds, mspe.cfust, fit.cfust$cpu.time)), file=paste0(SPATH, '/results/simu2/','(',n,'_',na.rate,'_cfust).csv'), sep=",",append=T, row.names=F, col.names=F, na = "NA")
 }
}
M=3
MC(n=300, q=2, na.rate=0.1, M)
MC(n=300, q=2, na.rate=0.2, M)
MC(n=300, q=2, na.rate=0.3, M)
MC(n=600, q=2, na.rate=0.1, M)
MC(n=600, q=2, na.rate=0.2, M)
MC(n=600, q=2, na.rate=0.3, M)
MC(n=900, q=2, na.rate=0.1, M)
MC(n=900, q=2, na.rate=0.2, M)
MC(n=900, q=2, na.rate=0.3, M)
