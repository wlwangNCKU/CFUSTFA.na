#rSNFA.na
sqrt.mt = function(S)
{
  p = ncol(S)
  if(p == 1) S.sqrt = as.matrix(sqrt(S))
  else{
    eig.S = eigen(S)
    S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve)}
}

drMSN = function(Y, mu, S, la, log = F)
{
  require(mvtnorm)
  if(is.vector(Y) == T) Y = as.matrix(Y)
  if(is.matrix(S) == F) S = as.matrix(S)
  p = ncol(Y)
  cent = t(Y) - mu
  om=S+la%*%t(la)
  xi=as.vector(t(la)%*%solve(om)%*%cent)
  si2=as.vector(1-t(la)%*%solve(om)%*%la)
  Phi = pnorm(xi/sqrt(si2))
  phi.p = dmvnorm(t(cent), mean = rep(0,p), sigma = om, log=FALSE)
  if(log == F) {den = 2 * phi.p * Phi}
  else {den = log(2) + log(phi.p) + log(Phi)}
  return(den)
}

#rSNFA.na-EM
rSNFA.na.EM = function(Y.na, q, init.para, zero.mu = F, zero.la = F, tol, max.iter, per)
{
  begin = proc.time()[1]
  n=nrow(Y.na)
  p=ncol(Y.na)
  la = init.para$la
  mu = init.para$mu
  B = init.para$B
  D = init.para$D
  # O matrices
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
  cc=sqrt(2/pi)
  Dh.inv=solve(sqrt.mt(diag(1,q)+(1-cc^2)*la%*%t(la)))
  B=B%*%Dh.inv
  S = B %*% t(B) + D
  d=c(B%*%la)
  mu.ast=mu-cc*d
  om=S+d%*%t(d)
  logsum=c()
  for(i in 1:num.na.class)
  {
    O = O.list[[i]]
    ind=ind.list[[i]]
    Y.pat = matrix(Y[ind,] , ncol = p)
    Y.o=Y.pat%*%t(O)
    OSO = O %*% S %*% t(O)       
    mu.ast.o=as.vector(O%*%mu.ast) 
    d.o=O%*%d  
    den = drMSN(Y.o, mu.ast.o, OSO, d.o, log = F)
    logsum[i]=sum(log(den))
  }
  logli.old=sum(logsum)
  cat(paste("rSNFA.na.EM is running...", sep = ""), "\n")
  iter = 0
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  cat("iter = ", iter, ",\t logli = ", logli.old, sep = "", "\n")
  repeat{
    iter=iter+1
    Y.hat = matrix(NA, p, n)
    hu=matrix(NA,q,n)
    hg1=hg2=matrix(NA,1,n)
    Phi.temp.Sum=matrix(0,p,p)
    Reta.sum=RPsi.sum=matrix(0,p,q)
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
      tmp = c(om.sq.inv %*% d.o)   
      a= colSums(tmp * e)
      ka2 = 1 - sum(tmp^2)
      ka = sqrt(ka2)
      A = a/ka
      hg1[,ind] = a + ka * dnorm(A)/pnorm(A)
      hg2[,ind]= ka2 + a * hg1[,ind]
      Woo = solve(diag(1,q) + t(B) %*% Coo %*% B)
      v.o = t(B) %*% Coo %*% (t(Y.pat)-mu)
      hu[,ind]=  Woo %*% (v.o + (la %*% t(hg1[,ind] - cc)))    
      xi = Woo %*%( t(t(v.o)*hg1[,ind]) + la %*% t(hg2[,ind] - cc * hg1[,ind])) 
      if(q==1)
        zeta1 = rowSums(xi) - cc * rowSums(as.matrix(t(hu[,ind])))
      else
        zeta1 = rowSums(xi) - cc * rowSums(as.matrix(hu[,ind]))
      zeta1.sum=zeta1.sum+zeta1
      hu2=(diag(length(ind),q)+hu[,ind]%*%t(v.o)+zeta1%*%t(la))%*%Woo
      hu2.sum=hu2.sum+hu2
      RPsi.temp=IBoo%*%hu2    
      RPsi.sum=RPsi.sum+RPsi.temp
      if(q==1)
        Reta.temp=IBoo%*%t(hu[,ind])%*%(hu[,ind])
      else
        Reta.temp=IBoo%*%hu[,ind]%*%t(hu[,ind])
      Reta.sum=Reta.sum+Reta.temp
      ze.hat = mu + B %*% hu[, ind]
      Y.hat[,ind] = ze.hat + DCoo %*% (t(Y.pat)- ze.hat)
      if(q==1)
        Lacov.sum =(hu2 - t(hu[, ind]) %*% (hu[, ind]))
      else
        Lacov.sum =(hu2 - hu[, ind] %*% t(hu[, ind]))
      Phi.temp =  length(ind) * IDoo + (IBoo-B)%*%Lacov.sum %*%t(IBoo-B)
      Phi.temp.Sum = Phi.temp.Sum + Phi.temp
    }
    zeta2 = sum(hg2) - 2 * cc * sum(hg1) + n * cc^2
    mu=rowMeans((Y.hat-(B%*%hu)))
    
    B=((Y.hat - mu)%*%t(hu)+ RPsi.sum -Reta.sum)%*% solve(hu2.sum)
    Y.mu = Y.hat - mu
    Y.temp = Y.mu - B %*% hu
    ttmp=Phi.temp.Sum +Y.temp%*%t(Y.temp)
    D = diag(diag(ttmp)/n)
    if (!zero.la) {
      la=zeta1.sum/zeta2
    }
    S = B %*% t(B) + D
    d=c(B%*%la)
    mu.ast=mu-cc*d
    for(i in 1:num.na.class)
    {
      O = O.list[[i]]
      OSO = O %*% S %*% t(O)
      ind = ind.list[[i]]
      Y.pat = matrix(Y[ind,] , ncol = p)
      Y.o=Y.pat%*%t(O)
      mu.ast.o=as.vector(O%*%mu.ast)
      d.o=O%*%d
      den = drMSN(Y.o, mu.ast.o, OSO, d.o, log = F)
      logsum[i]=sum(log(den))
    }
    logli.new=sum(logsum)
    diff = logli.new - logli.old
    if (iter%%per == 0) cat("iter =", iter, "\tlogli =", logli.new, "\tdiff =", diff, "\n")
    if (diff < tol | iter == max.iter) break
    logli.old = logli.new
  }
  Delta = diag(1,q) + (1-cc^2) * la %*% t(la)
  Dh=sqrt.mt(Delta)
  B = B %*% Dh
  Dh.inv=solve(Dh)
  hu=t( Dh.inv %*% hu)
  m = 2 * p + sum(la != 0) + p * q - q * (q - 1)/2
  AIC = -2 * logli.new + m * 2
  BIC = -2 * logli.new + m * log(n) 
  mds = c(no.para = m, logli = logli.new, AIC = AIC, BIC = BIC)
  na.ind = which(as.vector(t(na.posi))==T)
  Y.mis = as.vector(Y.hat)[na.ind]
  para = list(mu = mu,  B = B,  dd = diag(D), la = la, S = S)
  end = proc.time()[1]
  cat(paste(rep("=",60),sep="",collapse=""),"\n")
  cat("iter =", iter, "\tlogli =", logli.new, "\tdiff =", diff, "\n")
  cat('rSNFA.na.EM takes', end - begin, 'seconds.\n')
  list(mds=mds, iter=iter, para=para, Y.hat=t(Y.hat), score=hu, Y.mis = Y.mis, cpu.time=end - begin)
}
