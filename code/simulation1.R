library(cubature)
library(MomTrunc)
library(mvtnorm)
library(relliptical)
library(tmvtnorm)

################################################################################
#
#   Filename: simulation1.R
#   Purpose: Conduct simulation study for the CFUSTFA model with missing data,
#            estimate parameters using AECM algorithm, and store estimates and
#            standard errors across various sample sizes and missing rates
#   Input data files: function/cfustfn.R;
#                     function/cfustfa.na.R;
#                     function/cfustfa.na.se.simu.R
#   Output data files: results/simu1/para/<M>/(n_dropout)para.csv (parameter estimates)
#                      results/simu1/sd/<M>/(n_dropout)sd.csv   (standard errors)
#   R Version: R-4.4.1
#   Required R packages: cubature, relliptical, mvtnorm, MomTrunc, tmvtnorm
#
################################################################################

source(paste(SPATH, "/function/gener_na.R",sep=""))
source(paste(SPATH, "/function/cfustfn.R",sep=""))
source(paste(SPATH, "/function/cfustfa.na.R",sep=""))
source(paste(SPATH, "/function/cfustfa.na.se.simu.R",sep=""))

##------------------------------------------------------------------------
simu1 = function(M, n, na.rate)
{
  ##--------------------- pamareter values ----------------------------------
  p=5
  mu = c(4, 2, 6, 8, 4)
  B = t(matrix(c(0.13, 0.13, 0.97, 0.94, 0.54,0.84, 0.84, 0.27, 0.71, 0.77), nrow = 2, ncol = 5, byrow = TRUE))
  D = diag(c(9, 10, 2, 3, 7))
  La= matrix(c(7,3,1,4),2,2)
  nu = 4
  #true.para = c(mu=mu, B=B, d=diag(D), la=la, nu=nu)

  ##--------------------- generate data ----------------------------------
  for(i in 1:M){
    cat("Replication=", i, "\n")
    repeat{
      gen = rcfustfa(n,mu,B,D,La,nu)
      Y = gen$data
      Y.na = gener.na(Y, na.rate)
      init.para = list(mu=mu, B=B, D=D, La=La+matrix(runif(4, -1, 1),2,2), nu=nu)
      fit1 = try(CFUSTFANA.AECM(Y.na, q=2, s=2, init.para, zero.mu = F, zero.la = F, uST=F, tol = 0.001, max.iter = 150, per = 10),silent = TRUE)
      if(class(fit1)!="try-error") break
    }
    se = cfustna.se.new(Y.na, q=2, s=2, para.est=fit1$para)
    est.para = c(se[,1])
    est.se = c(se[,2])
    write.table(t(est.para),paste0(SPATH,"/results/simu1/para/","(",n,"_",na.rate,").csv"),sep=",",append=T,row.names=F,col.names=F, na = "NA")
    write.table(t(est.se),paste0(SPATH,"/results/simu1/sd/","(",n,"_",na.rate,").csv"),sep=",",append=T,row.names=F,col.names=F, na = "NA")
  }
}
simu1(100,150, na.rate=0)
simu1(100,150, na.rate=0.1)
simu1(100,150, na.rate=0.2)
simu1(100,150, na.rate=0.3)
simu1(100,300, na.rate=0)
simu1(100,300, na.rate=0.1)
simu1(100,300, na.rate=0.2)
simu1(100,300, na.rate=0.3)
simu1(100,600, na.rate=0)
simu1(100,600, na.rate=0.1)
simu1(100,600, na.rate=0.2)
simu1(100,600, na.rate=0.3)
simu1(100,900, na.rate=0)
simu1(100,900, na.rate=0.1)
simu1(100,900, na.rate=0.2)
simu1(100,900, na.rate=0.3)
simu1(100,1500, na.rate=0)
simu1(100,1500, na.rate=0.1)
simu1(100,1500, na.rate=0.2)
simu1(100,1500, na.rate=0.3)
