sqrt.mt = function(S)
{
  p = ncol(S)
  if(p == 1) S.sqrt = as.matrix(sqrt(S))
  else{
    eig.S = eigen(S)
    S.sqrt = eig.S$ve %*% diag(sqrt(eig.S$va)) %*% t(eig.S$ve) }
}

dCFUST = function(Y, mu.ast, Sigma, alpha, nu)
{
  p = ncol(Y)
  s = ncol(alpha) 
  Omega = Sigma + alpha %*% t(alpha) 
  Omega.inv = solve(Omega)  
  y.cent= t(Y) - mu.ast  
  Delta_Sigma = diag(1,s) -  t(alpha) %*% Omega.inv %*% alpha 
  h = t(alpha) %*% Omega.inv %*% y.cent  
  M = colSums((Omega.inv %*% y.cent)*y.cent) 
  A = t(h) * sqrt((nu+p)/(nu+M))
  NU=round(nu,0) 
  cdf = apply(t(A), 2, pmvt, lower = rep(-Inf, s), delta = rep(0, s), sigma = Delta_Sigma, df = NU+p) 
  denCFUST = 2^s * dmvt(Y, delta = mu.ast, sigma = Omega, df = nu, log = FALSE) * cdf
  return(denCFUST)
}

CFUSTFANA.AECM = function (Y.na, q, s, init.para, zero.mu = F, zero.la = F, uST=T, tol = 0.001, max.iter = 200, per = 5)
{
  method='AECM'
  begin = proc.time()[1]
  n = nrow(Y.na) 
  p = ncol(Y.na) 
  mu = init.para$mu  
  B  = init.para$B 
  D  = init.para$D   
  La =  init.para$La  
  nu = init.para$nu 
  iter = 0
  if (zero.mu)  mu = rep(0, p)
  if (zero.la)  la = rep(0, q)
  
  # O matrices
  na.posi = is.na(Y.na) 
  na.class = colSums(t(na.posi) * 2 ^ (0:(p-1))) 
  uni.na.class = unique(na.class) 
  num.na.class = length(uni.na.class)
  O.list = ind.list = as.list(num.na.class)
  for(k in 1:num.na.class)
  {
    ind.list[[k]]= which(na.class == uni.na.class[k])  
    O.list[[k]] = matrix(diag(p)[!na.posi[ind.list[[k]][1],],], ncol = p)  
  }
  Y.na[na.posi] = 99999
  
  #Observed likelihood function
  den = rep(NA,n)
  a_nu = sqrt(nu/pi)*exp(lgamma((nu-1)/2)-lgamma(nu/2)) 
  Js = matrix(1,s,s)
  vec1 = matrix(rep(1,s),nrow=s) 
  S = B %*% t(B) + D  
  S.inv = solve(S)
  alpha = B %*% La    
  mu.ast = mu - a_nu * alpha %*% vec1 
  
  for(i in 1:num.na.class)
  {
    O = O.list[[i]]
    ind = ind.list[[i]] 
    Y.pat = matrix(Y.na[ind,],ncol=p) 
    Y.o = Y.pat %*% t(O)  
    OSO = O %*% S %*% t(O) 
    mu.ast.o = c(O %*% mu.ast)  
    alpha.o = O %*% alpha  
    den[ind] = dCFUST (Y=Y.o, mu.ast.o, OSO,alpha.o, nu=nu)
  }  
  logli.old = sum(log(den)) 
  cat(paste("CFUSTFA.na.EM is running...", sep = ""), "\n")
  Qnu = function(Y.na, mu.ast, S, alpha, nu)
  {
    n = nrow(Y.na)
    den = rep(NA, n)
    for(i in 1:num.na.class)
    {
      O = O.list[[i]]
      ind = ind.list[[i]]
      Y.pat = matrix(Y.na[ind,] , ncol = p) 
      Y.o = Y.pat %*% t(O)  
      OSO = O %*% S %*% t(O)  
      mu.ast.o = c(O %*% mu.ast)  
      alpha.o = O %*% alpha   
      den[ind] = dCFUST(Y=Y.o, mu.ast.o, OSO, alpha.o, nu)
    } 
    indv.den = sum(log(den))
    return(-indv.den) 
  }
  if (zero.la == F) {
    cat( paste(rep('=', 10),sep='',collapse=''), method, ': CFUSTFAna: (p=', p, '; q=', q, '; s=', s, ')', '\t La_i~=0 for all i', paste(rep('=',10),sep='',collapse=''), '\n')
  } else {
    cat( paste(rep('=', 10),sep='',collapse=''), method, ': TFAna: (p=', p, '; q=', q, '; s=', s, ')', '\t La_i=0 for all i', paste(rep('=',10),sep='',collapse=''), '\n')
  }
  cat('iter = ', iter, '\t init.log.like = ', logli.old, sep = ' ', '\n')
  #----------------------------------------------------------------------------------------- 
  repeat {
    iter = iter + 1 
    tau = rep(NA,n)
    gamma1 = matrix(NA,n,q)         
    gamma2 = array(NA,dim=c(q,q,n)) 
    phi = matrix(NA, n, p)            
    zeta = array(NA, dim=c(p,q,n))   
    eta = matrix(NA,n,q)              
    H = Psi = array(NA,dim=c(q,q,n))  
    psi = array(NA, dim=c(p,q,n))     
    Phi = array(NA, dim=c(p,p,n))    
    om = S + alpha %*% t(alpha) 
    om.inv = solve(om)

    for(k in 1:num.na.class){   
      O = O.list[[k]]
      ind = ind.list[[k]]  
      po = nrow(O)  
      num = length(ind)
      
      mu.ast.o = O %*% mu.ast 
      alpha.o = O %*% alpha 
      OSO = O %*% S %*% t(O) 
      Oom = O %*% om %*% t(O)
      Oom.inv = solve(Oom)	
      Y.pat = matrix(Y.na[ind,], ncol = p) 
      Y.o = Y.pat %*% t(O) 
      cent.o = t(Y.o) - c(mu.ast.o)  
      Soo = t(O) %*% solve(OSO) %*% O 
      ISoo = diag(1,p) - S %*% Soo

      Delta_Sig = diag(1,s) - t(alpha.o) %*% Oom.inv %*% alpha.o
      Delta_Sigma = ( Delta_Sig + t(Delta_Sig) )/2 
      h = t(alpha.o) %*% Oom.inv %*% cent.o  
      M = colSums((Oom.inv %*% cent.o)*cent.o)   
      A1 = sqrt((nu+po+2)/(nu+M)) * t(h)  
      A2 = sqrt((nu+po)/(nu+M)) * t(h)   
      NU = round(nu,0) 
      cdf1 = apply(t(A1), 2, pmvt, lower = rep(-Inf, s), delta = rep(0, s), sigma = Delta_Sigma, df = NU+po+2)  
      cdf2 = apply(t(A2), 2, pmvt, lower = rep(-Inf, s), delta = rep(0, s), sigma = Delta_Sigma, df = NU+po) 
      tau[ind] = ((po+nu)/(M+nu))*(cdf1/cdf2) 
      
      for (j in 1:num){
        inx = ind[j]
        ka2 = (nu+M[j])/(nu+po+2) * Delta_Sigma 
        tmt.mom = meanvarTMD(lower = rep(0 ,s), upper = rep(Inf, s), mu=h[,j], Sigma=ka2, nu=nu+po+2, dist="t")
        hg1 = tmt.mom$mean
        hg2 = tmt.mom$EYY
        gamma2[,,inx] = ( tau[inx]*hg2 + t(tau[inx]*hg2) )/2
        if ( q > 1 ){
          if ( !isSymmetric(gamma2[,,inx]) ) {cat('num=',inx,'\t gamma2=', '\n');print(gamma2[,,inx])}
        }
        gamma1[inx,] = tau[inx] * hg1
        phi[inx,] = t( tau[inx]*S%*%Soo%*%Y.pat[j,] +  ISoo%*%(mu.ast*tau[inx]+ alpha%*%gamma1[inx,]))
        zeta[,,inx] = S%*%Soo%*%Y.pat[j,]%*%t(gamma1[inx,]) + ISoo%*%( mu.ast%*%t(gamma1[inx,]) + alpha%*%as.matrix(gamma2[,,inx]))
      }
    }
    total_tau = sum(tau)
    total_gammal = colSums(gamma1)
    total_gamma2 = rowSums(gamma2, dims=2)
    total_phi = colSums(phi)
    total_zeta = rowSums(zeta, dims=2)
    
    
    if (!zero.mu) mu = ( total_phi - c(alpha %*% total_gammal) + c(a_nu*total_tau*t(alpha %*% vec1)) ) / total_tau
    if (q==1) alpha_a = a_nu*total_tau*mu%*%t(vec1) - a_nu*total_phi%*%t(vec1) - mu*total_gammal + total_zeta
    if (q!=1) alpha_a = a_nu*total_tau*mu%*%t(vec1) - a_nu*total_phi%*%t(vec1) - mu%*%t(total_gammal) + total_zeta
    alpha_b = -a_nu^2*total_tau*Js  +  a_nu*total_gammal%*%t(vec1) + a_nu*vec1%*%t(total_gammal) - total_gamma2
    alpha = alpha_a %*% solve(-alpha_b)
    
    mu.ast = mu - a_nu * alpha %*% vec1 
  if(uST==T)   
  {
   if (!zero.la && q==1)  La = solve(t(B)%*%B) %*% t(B) %*% alpha
   else if (!zero.la && q>1)   La = diag( diag( solve(t(B)%*%B) %*% t(B) %*% alpha ), q)
   }
  if(uST==F) La = solve(t(B)%*%B) %*% t(B) %*% alpha
  
    #---------------------------------------------
    for(k in 1:num.na.class){   
      O = O.list[[k]]
      ind = ind.list[[k]]  
      po = nrow(O)  
      num = length(ind)
      
      OSO = O %*% S %*% t(O)
      ODO = O %*% D %*% t(O)
      Y.pat = matrix(Y.na[ind,], ncol = p) 
      cent.mu = t(Y.pat) - mu
      Soo = t(O) %*% solve(OSO) %*% O 
      Coo = t(O) %*% solve(ODO) %*% O
      ISoo = diag(1,p) - S %*% Soo
      ICoo = diag(1,p) - D %*% Coo
      
      W = solve( diag(1,q) + t(B) %*% Coo %*% B)  
      V = t(B) %*% Coo %*% cent.mu  
      for (j in 1:num){
        inx = ind[j]
        b = V[,j] - a_nu * La %*% vec1
        eta[inx,] = W %*% ( b * tau[inx] + La %*% gamma1[inx,] )
        H[,,inx] = W %*% ( b %*% t(gamma1[inx,]) + La %*% gamma2[,,inx] )
        Psi[,,inx] = ( diag(1,q) + eta[inx,] %*% t(b) + H[,,inx] %*% t(La) ) %*% W
        Phi[,,inx] = phi[inx,]%*%t(Y.pat[j,])%*%Soo%*%S + ( phi[inx,]%*%t(mu.ast) + zeta[,,inx]%*%t(alpha) ) %*% t(ISoo) + ISoo %*% S
        psi[,,inx] = ICoo %*% ( mu%*%t(eta[inx,]) + B%*%Psi[,,inx] ) + D %*% Coo %*% Y.pat[j,] %*% t(eta[inx,])
      }
    }      
    total_eta = colSums(eta)
    total_Psi = ( rowSums(Psi, dims=2) + t(rowSums(Psi, dims=2)) )/2
    total_Phi = ( rowSums(Phi, dims=2) + t(rowSums(Phi, dims=2)) )/2
    total_psi = rowSums(psi, dims=2)
    
    
    B = (total_psi - mu %*% t(total_eta) ) %*% solve(total_Psi)
    
    
    D.part_phi = total_Phi - total_phi%*%t(mu) - mu%*%t(total_phi) + total_tau*mu%*%t(mu) 
    D.part_psi = (mu%*%t(total_eta)-total_psi)%*%t(B) + B%*%t(mu%*%t(total_eta)-total_psi) + B%*%total_Psi%*%t(B)
    D =diag( diag(D.part_phi + D.part_psi) / n)
    S = B %*% t(B) + D
    alpha = B %*% La
    mu.ast = mu - a_nu * alpha %*% vec1 
    
    
    nu = optim(par=nu, fn=Qnu, method="L-BFGS-B", lower=1e-3, upper=Inf, Y.na=Y.na, mu.ast=mu.ast, S=S, alpha=alpha)$par   
    
    a_nu = sqrt(nu/pi)*exp(lgamma((nu-1)/2)-lgamma(nu/2))
    for(i in 1:num.na.class)
    {   
      O = O.list[[i]]
      ind = ind.list[[i]]
      Y.pat = matrix(Y.na[ind,] , ncol = p) 
      Y.o = Y.pat %*% t(O)  
      OSO = O %*% S %*% t(O)  
      mu.ast.o = c(O %*% mu.ast) 
      alpha.o = O %*% alpha 
      den[ind] = dCFUST (Y=Y.o, mu=mu.ast.o,Sigma= OSO, alpha=alpha.o, nu)
    } 
    
    logli.new = sum(log(den)) 
    
    diff = logli.new - logli.old
    if (iter%%per == 0) cat("iter =", iter, "\tlogli =", logli.new, "\tdiff =", diff, "\tLa=", La, "\tnu=", nu, "\n")
    if (diff < tol | iter == max.iter) break
    logli.old = logli.new
  }
  if(uST==T) la=diag(La)
  else la=c(La)
  para = list(mu = mu,  B = B,  dd = diag(D), la = la, nu = nu)
  cat(paste(rep("=", 60), sep = "", collapse = ""), "\n")
  cat("iter =", iter, "\tlogli =", logli.new, "\tdiff =", diff, "\t la=", la, "\t nu=", nu, "\n")

  U = matrix(NA, n, q)
  Y.est = matrix(NA, n, p)
  Y.hat = matrix(NA, n, p)
  
  om = S + alpha %*% t(alpha) 
  om.inv = solve(om)
  
  for(k in 1:num.na.class){   
    O = O.list[[k]]
    ind = ind.list[[k]] 
    po = nrow(O)   
    num = length(ind)
    
    mu.ast.o = O %*% mu.ast 
    alpha.o = O %*% alpha
    OSO = O %*% S %*% t(O)
    ODO = O %*% D %*% t(O)
    Oom = O %*% om %*% t(O)
    Oom.inv = solve(Oom)	
    Y.pat = matrix(Y.na[ind,], ncol = p) 
    Y.o = Y.pat %*% t(O) 
    
    cent.o = t(Y.o) - c(mu.ast.o) 
    cent.mu = t(Y.pat) - mu
    
    Soo = t(O) %*% solve(OSO) %*% O 
    Coo = t(O) %*% solve(ODO) %*% O
    ISoo = diag(1,p) - S %*% Soo
    ICoo = diag(1,p) - D %*% Coo
    
    W = solve( diag(1,q) + t(B) %*% Coo %*% B)  
    V = t(B) %*% Coo %*% cent.mu  
    Delta_Sig = diag(1,q) - t(alpha.o) %*% Oom.inv %*% alpha.o
    Delta_Sigma = ( Delta_Sig + t(Delta_Sig) )/2  
    h = t(alpha.o) %*% Oom.inv %*% cent.o  
    M = colSums((Oom.inv %*% cent.o)*cent.o)  
    A1 = sqrt((nu+po+2)/(nu+M)) * t(h)  
    A2 = sqrt((nu+po)/(nu+M)) * t(h)  
    NU = round(nu,0) 
    cdf1 = apply(t(A1), 2, pmvt, lower = rep(-Inf, q), delta = rep(0, q), sigma = Delta_Sigma, df = NU+po+2)  
    cdf2 = apply(t(A2), 2, pmvt, lower = rep(-Inf, q), delta = rep(0, q), sigma = Delta_Sigma, df = NU+po)     
    tau[ind] = ((po+nu)/(M+nu))*(cdf1/cdf2)  
  
    for (j in 1:num){
      inx = ind[j]
      gy.ka2 = (nu+M[j])/(nu+po+2) * Delta_Sigma
      gy.tmt.mom = meanvarTMD(lower = rep(0 ,q), upper = rep(Inf, q), mu=h[,j], Sigma=gy.ka2, nu=nu+po, dist="t")
      gy.hg1 = gy.tmt.mom$mean
      U[inx,] = t(W %*% ( V[,j] + La %*% gy.hg1 - a_nu*La%*%vec1 ) )
      Y.hat[inx,] = mu + B %*% U[inx,] + D %*% Coo %*% Y.pat[j,] - D %*% Coo %*% mu - D %*% Coo %*% B %*% U[inx,]
      Y.est[inx,] = t( mu + B%*% U[inx,] ) 
    }
  }    
  U = U 
  na.ind = which(as.vector(t(na.posi))==T)
  Y.mis = as.vector(t(Y.hat))[na.ind]

  nona.ind = which(as.vector(t(na.posi))==F)
  Y.obs = as.vector(t(Y.na))[nona.ind]
  y.est.obs = as.vector(t(Y.est))[nona.ind]  
  mse = mean((Y.obs - y.est.obs)^2, na.rm=T)
  ##=======
  #m = p + ( p*q-q*(q-1)/2 ) + p + sum(la != 0) + 1 #num( mu+B+D+La+nu )
  m = p + ( p*q-q*(q-1)/2 ) + p + ifelse(uST, q, q*s) + 1 #num( mu+B+D+La+nu )
  AIC = -2 * logli.new + m * 2
  BIC = -2 * logli.new + m * log(n) 
  #mds = c(no.para = m, logli = logli.new, AIC = AIC, BIC = BIC,mse = mse)
  mds = c(no.para = m, logli = logli.new, AIC = AIC, BIC = BIC)
  ##=====
  end = proc.time()[1]
  cpu.time=end - begin
  cat("CFUSTFA.na.AECM takes", cpu.time, "seconds\n\n")
  return(list(method = method, iter = iter, q = q, s=s, mds = mds, para = para, mse = mse, factor.score = U, y.est = Y.est, Y.mis = Y.mis, cpu.time=cpu.time))
}
