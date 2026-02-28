#library(mvtnorm)
#library(mnormt)
library(relliptical)
#library(rgl)
#library(MASS)
#library(ks)
#library(misc3d) #countour3D

sqrt.mt = function(S)
{
  p = ncol(S)
  if(p == 1) S.sqrt = as.matrix(sqrt(S))
  else{
    eig.S = eigen(S)
    S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve) }
}

rumst=function(n, xi, S, Lambda, nu) {
  p=length(xi)
  s = ncol(Lambda) 
  tau=rgamma(n,nu/2 , nu/2) 
  zeta = t(rtelliptical(n ,mu=rep(0,s), Sigma=diag(s), lower=rep(0,s), upper=rep(Inf,s), dist="Normal"))
  U= t(rmvnorm(n, rep(0,p), S)) 
  W= Lambda %*% zeta + U 
  y=t(xi + t(t(W)/sqrt(tau))) 
  return(y)
}

rcfustfa=function(n,mu,B,D,La,nu)
{ 
  p = nrow(B)
  q = ncol(B)
  #La = diag(la,q,q)
  a_nu = sqrt(nu/pi)*exp(lgamma((nu-1)/2)-lgamma(nu/2))
  alpha = B %*% La
  xi = c(mu - a_nu*alpha %*% rep(1,q))
  BB = B %*% t(B)
  S = BB + D
  data = rumst(n=n, xi=xi, S=S, Lambda=alpha, nu=nu)
  return(list(data=data,xi=xi,BB=BB,S=S,alpha=alpha,a_nu=a_nu))
}


dCFUST = function(Y, mu.ast, Sigma, alpha, nu) {
  p = ncol(Y)
  s = ncol(alpha)
  Omega = Sigma + alpha %*% t(alpha)
  Omega.inv = solve(Omega)
  y.cent = t(Y) - mu.ast
  Delta_Sigma = diag(1, s) - t(alpha) %*% Omega.inv %*% alpha
  h = t(alpha) %*% Omega.inv %*% y.cent
  M = colSums((Omega.inv %*% y.cent) * y.cent)
  A = t(h) * sqrt((nu + p) / (nu + M))
  NU = round(nu, 0)
  cdf = apply(t(A), 2, pmvt, lower = rep(-Inf, s), delta = rep(0, s),
              sigma = Delta_Sigma, df = NU + p)
  density_1 = 2^s * dmvt(Y, delta = mu.ast, sigma = Omega, df = nu, log = FALSE) * cdf
  return(density_1)
}



