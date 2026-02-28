#rSTFA.na
sqrt.mt = function(S)
{
  p = ncol(S)
  if(p == 1) S.sqrt = as.matrix(sqrt(S))
  else{
    eig.S = eigen(S)
    S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)}
}


drMST = function(Y, mu, S, la, nu, log = F)
{
  require(mvtnorm)
  if(is.vector(Y) == T) Y = as.matrix(Y)
  if(is.matrix(S) == F) S = as.matrix(S)
  p = length(la)
  cent = t(Y) - mu
  om = S + la %*% t(la)
  eigom=eigen(om)
  omh.inv = eigom$ve %*% diag(1/sqrt(eigom$va), p) %*% t(eigom$ve)
  e = omh.inv %*% cent
  Aj = colSums(e^2)
  tmp = c(omh.inv %*% la)
  xi = colSums(tmp * e)
  si2 = 1 - sum(tmp^2)
  r0 = (nu + p) / (nu + Aj)
  A = xi / sqrt(si2)
  c0 = A * sqrt(r0) ; Tc0 = pt(c0, df = nu + p)
  logtp = lgamma((nu+p) / 2) - lgamma(nu / 2) - log(det(om)) / 2 - p / 2 * log(pi*nu) - (nu + p) / 2 * log(1 + Aj / nu)
  if(log == F){den =2*exp(logtp)*Tc0}
  else {den = log(2) + logtp + log(Tc0)}
  return(den)
}

#rSTFA.na-EM
rSTFA.na.EM =function(Y.na, q, init.para, zero.mu = F, zero.la=F, tol, max.iter, per)
{
  begin = proc.time()[1]
  n=nrow(Y.na)
  p=ncol(Y.na)
  la = init.para$la
  mu = init.para$mu
  B = init.para$B
  D = init.para$D
  nu=init.para$nu
  na.posi = is.na(Y.na)                          
  na.class = colSums(t(na.posi) * 2 ^ (0:(p-1)))     
  uni.na.class = unique(na.class)                    
  num.na.class = length(uni.na.class)                
  Ip = diag(p)
  O.list = ind.list = as.list(num.na.class)                                                                      
  for(i in 1:num.na.class)
  {
    ind.list[[i]]= which(na.class == uni.na.class[i])
    O.list[[i]] = matrix(diag(p)[!na.posi[ind.list[[i]][1],],], ncol = p)
  }
  Y=Y.na
  Y[is.na(Y.na)] = 99999
  if (zero.mu)   mu = rep(0, p)
  if (zero.la)   la = rep(0, q)
  anu = sqrt(nu / pi) * exp(lgamma((nu-1) / 2) - lgamma(nu/2))
  Delta = diag(1, q) + (1 - (nu - 2) / nu * anu^2) * la %*% t(la)
  Dh = sqrt.mt(Delta)
  eigDelta=eigen(Delta)
  Dh.inv = eigDelta$ve %*% diag(1/sqrt(eigDelta$va), q) %*% t(eigDelta$ve)
  B=B%*%Dh.inv
  S = B %*% t(B) + D
  d=c(B%*%la) 
  mu.ast=mu - anu * d
  logsum=c()
  for(i in 1:num.na.class)
  {   
    O = O.list[[i]]
    ind=ind.list[[i]]  
    Y.pat = matrix(Y[ind,] , ncol = p)
    Y.o=Y.pat%*%t(O)
    OSO = O %*% S %*% t(O)       #Sigma.oo#
    mu.ast.o=as.vector(O%*%mu.ast) #mu.o-cd.o#
    d.o=O%*%d  #do#
    den = drMST(Y.o, mu.ast.o, OSO, d.o, nu, log = F) 
    logsum[i]=sum(log(den))
  }
  init.LL = sum(logsum)
  iter = 1
  cat(paste("rSTFA.na.EM is running...", sep = ""), "\n")
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  cat("iter =", iter-1, "\t init.LL =", init.LL, "\n")
  
  ell.fn = function(Y, mu.ast, S, d, nu)
  {
    logsum=c()
    for(i in 1:num.na.class)
    {   
      O = O.list[[i]]
      ind=ind.list[[i]]  
      Y.pat = matrix(Y[ind,] , ncol = p)
      Y.o=Y.pat%*%t(O)
      OSO = O %*% S %*% t(O)       #Sigma.oo#
      mu.ast.o=as.vector(O%*%mu.ast) #mu.o-cd.o#
      d.o=O%*%d  #do#
      den = drMST(Y.o, mu.ast.o, OSO, d.o, nu, log = F) 
      logsum[i]=sum(log(den))
    }
    -sum(logsum)
  }
  
  LL =  lk = init.LL
  epsilon= Inf
  
  while((iter <max.iter) && (epsilon > tol))
  {
    Y.hat=Reta=b.hat =DPS =matrix(NA, p, n)
    hu=u1j=matrix(NA,q,n)
    tau=ga=g1j=g2j=matrix(NA,1,n)
    Phi.temp.Sum=matrix(0,p,p)
    RPsi.sum=matrix(0,p,q)
    hu2.sum=matrix(0,q,q)
    zeta1.sum=rep(0,q)
    for(i in 1:num.na.class)
    {  
      ind = ind.list[[i]]
      O = O.list[[i]]
      po=nrow(O)
      om=S+d%*%t(d)
      OMO= O %*% om %*% t(O)
      ODO= O %*% D %*% t(O)
      Y.pat = matrix(Y[ind,] , ncol = p)
      Y.o=Y.pat%*%t(O)
      mu.ast.o=as.vector(O%*%mu.ast)
      d.o=as.vector(O%*%d)  
      mu.o=as.vector(O%*%mu)
      
      OMoo=t(O) %*% solve(OMO) %*% O  
      Coo=t(O) %*% solve(ODO) %*% O
      DCoo=D%*%Coo
      IBoo=(diag(p)-DCoo)%*%B 
      IDoo=(diag(p)-DCoo)%*%D
      cent.o =t(Y.o) - mu.ast.o
      cent= Y.pat - mu.ast  
      
      eigom=eigen(OMO)
      om.sq.inv = eigom$ve %*% diag(1/sqrt(eigom$va), po) %*% t(eigom$ve)
      e = om.sq.inv %*% cent.o      
      Aj = colSums(e^2)
      tmp = c(om.sq.inv %*% d.o)  
      xi= colSums(tmp * e) 
      si2 = 1 - sum(tmp^2)
      A = xi / sqrt(si2)    
      rn2 = (nu + po - 2) / (nu + Aj); cn2 = A * sqrt(rn2)
      r0 = (nu + po) / (nu + Aj); c0 = A * sqrt(r0)
      r2 = (nu + po + 2) / (nu + Aj); c2 = A * sqrt(r2)
      Tc0 = pt(c0, df = nu + po); Tc2 = pt(c2, df = nu + po + 2)
      tc0 = dt(c0, df = nu + po); tcn2 = dt(cn2, df = nu + po - 2)
      tau[,ind] = r0 * Tc2 / Tc0
      g1j[,ind] = xi * tau[,ind] + sqrt(si2) * sqrt(r0) * tc0 / Tc0
      g2j[,ind] = si2 + xi * g1j[,ind]
      ga[,ind] = xi + sqrt(si2) * sqrt(rn2) * tcn2 / Tc0
      Woo = solve(diag(1, q) + t(B) %*% Coo %*% B) 
      v.o = t(B) %*% Coo %*% (t(Y.pat)-mu)
      hu[,ind] = Woo %*% (v.o + (la) %*% t(ga[,ind] - anu))
      u1j[,ind]=Woo%*%(t(t(v.o)*tau[,ind])+(la)%*%t(g1j[,ind]-anu*tau[,ind]))
      u2j = Woo%*%(t(t(v.o)*g1j[,ind])+(la)%*%t(g2j[,ind]-anu*g1j[,ind]))
      if(q==1)
        zeta1 = rowSums(u2j) - anu * rowSums(as.matrix(t(u1j[,ind])))
      else
        zeta1 = rowSums(u2j) - anu * rowSums(as.matrix(u1j[,ind]))
      zeta1.sum=zeta1.sum+zeta1
      hu2=(zeta1%*%t(la)+u1j[,ind]%*%t(v.o)+diag(length(ind),q))%*%Woo
      hu2.sum=hu2.sum+hu2
      RPsi.temp=IBoo%*%hu2    
      RPsi.sum=RPsi.sum+RPsi.temp
      Reta[,ind]= DCoo%*%B%*%(u1j[,ind])
      b.hat[,ind] = mu + DCoo %*% (t(Y.pat)- mu)
      DPS[,ind]=IBoo%*%u1j[,ind]
      ze.hat = mu + B %*% hu[, ind]
      Y.hat[,ind] = ze.hat + DCoo %*% (t(Y.pat)- ze.hat)
      Phi.temp =  length(ind) * IDoo + IBoo%*%hu2%*%t(IBoo)
      Phi.temp.Sum = Phi.temp.Sum + Phi.temp  
    }  
    zeta2 = sum(g2j) - 2 * anu * sum(g1j) + anu^2 * sum(tau)
    mu= rowSums(t(c(tau)*t(b.hat))- Reta)/sum(tau)
    B=((b.hat - mu)%*%t(u1j)+ RPsi.sum)%*% solve(hu2.sum)
    A.b = b.hat - mu
    ttmp=Phi.temp.Sum- RPsi.sum%*%t(B) -B%*% t(RPsi.sum)+B%*% hu2.sum%*%t(B)
    Tmp=ttmp+A.b%*%t(DPS) -A.b%*%t(B%*%u1j)+DPS%*%t(A.b)-B%*%u1j%*%t(A.b)+(A.b)%*%(c(tau)*t(A.b))
    D = diag(diag(Tmp)/n)               
    if (!zero.la) {la=zeta1.sum/zeta2}
    S = B %*% t(B) + D
    d=c(B%*%la) 
    mu.ast=mu-anu*d
    nu = optim(par = nu, fn = ell.fn, method = "L-BFGS-B", lower = 1e-3, upper = Inf, Y = Y, mu = mu.ast, S = S, d = d)$par
    anu = sqrt(nu / pi) * exp(lgamma((nu-1) / 2) - lgamma(nu/2))
    for(i in 1:num.na.class)
    {   
      O = O.list[[i]]
      OSO = O %*% S %*% t(O)
      ind = ind.list[[i]]
      Y.pat = matrix(Y[ind,] , ncol = p)
      Y.o=Y.pat%*%t(O)
      mu.ast.o=as.vector(O%*%mu.ast)
      d.o=O%*%d     
      den = drMST(Y.o, mu.ast.o, OSO, d.o, nu, log = F) 
      logsum[i]=sum(log(den))
    }
    newLL =sum(logsum)
    lk = c(lk, newLL)
    if (iter < 2) epsilon = abs(LL-newLL)/abs(newLL)
    else 
    {
      tmp = (newLL - LL)/(LL-lk[iter-1])
      tmp2 =  LL + (newLL-LL)/(1-tmp)
      epsilon = abs(tmp2-newLL)
    }
    diff = newLL - LL
    if(diff < 0)  cat('iter=', iter, "logL is decreasing!\n")
    LL= newLL
    if (iter%%per == 0) 
      cat("iter =", iter, "\t newLL =", newLL, "\t diff =", diff,  "\t Aikten's diff =", epsilon, "\tla=", la, "\tnu=", nu, "\n")
    iter = iter + 1
  }
  cat(paste(rep("=", 60), sep = "", collapse = ""), "\n")
  cat("iter =", iter, "\t logli =", newLL, "\t diff =", diff, "\t Aikten's diff =", epsilon, "\tla=", la, "\tnu=", nu, "\n")
  
  Delta = diag(1, q) + (1 - (nu - 2) / nu * anu^2) * la %*% t(la)
  Dh=sqrt.mt(Delta)
  eigDelta=eigen(Delta)
  Dh.inv = eigDelta$ve %*% diag(1/sqrt(eigDelta$va), q) %*% t(eigDelta$ve)
  B = B %*% Dh
  hu = Dh.inv %*% hu
  no.para = 2 * p + sum(la != 0) + p * q - q * (q - 1)/2 + 1 
  AIC = -2 * newLL + no.para * 2
  BIC = -2 * newLL + no.para * log(n)    
  mds=c(no.para=no.para, logli =newLL, AIC = AIC, BIC = BIC)
  na.ind = which(as.vector(t(na.posi))==T)
  Y.mis = as.vector(Y.hat)[na.ind]
  para = list(mu = mu,  B = B,  dd = diag(D), la = la, nu=nu, S = S)
  end = proc.time()[1]
  cat("rSTFA.na.ECM takes", end - begin, "seconds\n\n")
  list(mds=mds, iter=iter, para=para, score=t(hu), Yhat=t(Y.hat), Y.mis = Y.mis, cpu.time=end - begin)
}
