cfustna.se = function(Y.na, q, s, para.est)
{
  n = nrow(Y.na)
  p = ncol(Y.na)
  mu=para.est$mu
  B=para.est$B
  dd=para.est$dd
  D=diag(dd)
  la=para.est$la
  La=diag(la, q)
  nu=para.est$nu
  na.posi = is.na(Y.na) #遺失值位置T/F
  po=p-rowSums(na.posi)#每列觀察到的變數個數
  na.class = colSums(t(na.posi) * 2 ^ (0:(p-1))) #每個觀測值遺失pattern
  uni.na.class = unique(na.class) #遺失pattern
  num.na.class = length(uni.na.class) #遺失pattern個數
  O.list = ind.list = as.list(num.na.class)
  for(k in 1:num.na.class)
  {
    ind.list[[k]]= which(na.class == uni.na.class[k])  
    O.list[[k]] = matrix(diag(p)[!na.posi[ind.list[[k]][1],],], ncol = p) 
  }
  #Y=Y.na
  Y.na[na.posi] = 99999
  
  #---------------------------------------------------------------------------
  a_nu = sqrt(nu/pi)*exp(lgamma((nu-1)/2)-lgamma(nu/2)) 
  Js = matrix(1,s,s)
  vec1 = matrix(rep(1,s),nrow=s) 
  S = B %*% t(B) + D  #pxp
  S.inv = solve(S)    #pxp
  alpha = B %*% La    #pxs
  mu.ast = mu - a_nu * alpha %*% vec1 #mu-a_nu*alpha*1s
  tau = rep(NA,n)
  gamma1 = matrix(NA,n,q)           #E(tau*ga|yjo)=E(tau|yjo)E(ga.tilde|yjo) 
  gamma2 = array(NA,dim=c(q,q,n))   #E(tau*ga*ga'|yjo)=E(tau|yjo)E(ga.tilde*ga.tilde'|yjo) #qxqxn
  phi = matrix(NA, n, p)            #E(tau Y|yjo)
  zeta = array(NA, dim=c(p,q,n))    #E(tau Y gamma'|yjo) 
  eta = matrix(NA,n,q)              #E(tau U|yjo)
  H = Psi = array(NA,dim=c(q,q,n))  #E(tau U gamma'|yjo),E(tau U U'|yjo)  dim=s*s*n
  psi = array(NA, dim=c(p,q,n))     #E(tau Y U'|yjo)  dim=s*s*n
  Phi = array(NA, dim=c(p,p,n))      #E(tau Y Y'|yjo)
  log.tau = rep(NA,n)
  om = S + alpha %*% t(alpha) #Y~cfustp(mu-a_nu*alpha*1q=mu.ast, Sigma=S, Lambda=alpha, nu=nu) #pxp
  om.inv = solve(om)
  
  for(k in 1:num.na.class){   
    O = O.list[[k]]
    ind = ind.list[[k]]  
    po = nrow(O)  #p
    num = length(ind)
    
    mu.ast.o = O %*% mu.ast #dim=pjo*1
    alpha.o = O %*% alpha #dim=pjo*q
    OSO = O %*% S %*% t(O) #pxp
    Oom = O %*% om %*% t(O)
    Oom.inv = solve(Oom)	
    Y.pat = matrix(Y.na[ind,], ncol = p) 
    Y.o = Y.pat %*% t(O) #dim=n*pjo
    
    cent.o = t(Y.o) - c(mu.ast.o)  # yj-(mu-a_nu*alpha*1q) #dim=pjo*n
    Soo = t(O) %*% solve(OSO) %*% O 
    ISoo = diag(1,p) - S %*% Soo
    
    #tau|yjo~gamma()
    Delta_Sig = diag(1,s) - t(alpha.o) %*% Oom.inv %*% alpha.o
    Delta_Sigma = ( Delta_Sig + t(Delta_Sig) )/2  #Gammai
    h = t(alpha.o) %*% Oom.inv %*% cent.o  #q*n
    M = colSums((Oom.inv %*% cent.o)*cent.o)  #n 
    A1 = sqrt((nu+po+2)/(nu+M)) * t(h)  #n*s  #sqrt((nu+pjo+2)/(nu+M)),vector:500
    A2 = sqrt((nu+po)/(nu+M)) * t(h)   #n*s
    NU = round(nu,0) #pmvt_df --(n,0>c<F
    cdf1 = apply(t(A1), 2, pmvt, lower = rep(-Inf, s), delta = rep(0, s), sigma = Delta_Sigma, df = NU+po+2)  #$@$l n
    cdf2 = apply(t(A2), 2, pmvt, lower = rep(-Inf, s), delta = rep(0, s), sigma = Delta_Sigma, df = NU+po)  #$@%@ n    
    tau[ind] = ((po+nu)/(M+nu))*(cdf1/cdf2)  #n*g
    
    f=NULL
    #ga.tilde|yjo~Tts(mu=hj,cov=(nu+Mj/nu+pjo+2)*Delta_Sigma,nu+pjo+2)
    for (j in 1:num){  
      #E(log_tau)
      inx = ind[j]
      con = digamma( (po + nu + s)/2 ) - digamma((po + nu)/2) - s/(nu + M[j]) #g_nu() function
      f[j] = adaptIntegrate( function(x) dmvt( x, delta = rep(0,s), sigma = (nu+M[j])/(nu+po)*Delta_Sigma, df = nu+p, log = FALSE)* (con-log(1+t(x)%*%solve(Delta_Sigma)%*%x/(M[j]+nu))+(nu+po+s)*t(x)%*%solve(Delta_Sigma)%*%x/((M[j]+nu+t(x)%*%solve(Delta_Sigma)%*%x)*(M[j]+nu)) ), rep(-100000,q), c(h[,j]), maxEval=2000000,absError=10e-5,tol=1e-5 )$integral
      # as.vector(unlist(f))
      log.tau[inx] = tau[inx] - log((nu + M[j]) / 2) - ((nu+po)/(nu+M[j])) + digamma((nu + po) / 2) + 1 / cdf2[j] * f[j]
       
      ka2 = (nu+M[j])/(nu+po+2) * Delta_Sigma  #sxs
      ###NU=round(nu,0) pmvt_df --(n,0>c<F ###meanvarTMD_Student's t-distribution_nu --(n,0%?9j<F
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
    
    #U|yjo,gamma,tau~Nq()
    W = solve( diag(1,q) + t(B) %*% Coo %*% B)  #q*q
    V = t(B) %*% Coo %*% cent.mu   #dim=s*n
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
  
  ################################################################################
  da_nu = 1/2*sqrt(1/(pi*nu))*exp(lgamma((nu-1)/2)-lgamma(nu/2))+1/2*sqrt(nu/pi)*exp(lgamma((nu-1)/2)-lgamma(nu/2))*(digamma((nu-1)/2)-digamma(nu/2))
  
  vech = outer(1:p,1:q) >= 0
  
  smu = sb = sdr = sla = snu = NULL
  for(j in 1:n) 
  {
    tmp_m = phi[j,] - mu * tau[j] - B %*% eta[j,] 
    mj = solve(D) %*% tmp_m   
    smu = rbind(smu,t(mj)) 
    
    Bj= solve(D) %*% (psi[,,j] - mu %*% t(eta[j,])- B%*% Psi[,,j])
    sb = rbind(sb,Bj[vech]) #p*q 
    
    #dim(tmp_d)=p*p
    tmp_d = Phi[,,j] - phi[j,] %*% t(mu) - mu %*% t(phi[j,]) - psi[,,j] %*% t(B) - B %*% t(psi[,,j]) + mu %*% t(eta[j,]) %*% t(B) + B%*% eta[j,] %*% t(mu) + tau[j] * mu %*% t(mu) + B %*% Psi[,,j] %*% t(B)
    Dj = -1/2 * (solve(D) - solve(D) %*% tmp_d %*% solve(D) )
    sdr = rbind(sdr,diag(Dj)) #n*p
    
    Laj = -a_nu * eta[j,] %*% t(rep(1,s)) + H[,,j] - a_nu^2 * tau[j] * La %*% Js + a_nu * La %*% rep(1,s) %*% t(gamma1[j,]) + a_nu * La %*% gamma1[j,] %*% t(rep(1,s)) - La %*% gamma2[,,j]
    sla = rbind(sla, diag(Laj)) #n*p
    
    nu1 = log.tau[j] -tau[j] + log(nu/2) + 1 - digamma(nu/2)
    #nu2 = t(eta[,j]) %*% la - t(gamma1[,j]) %*% La %*% la + t(la) %*% eta[,j] - t(la) %*% La %*% gamma1[,j] + 2*a_nu*tau[j] %*% t(la) %*% la
    
    nuu2 = (eta[j] - La %*% gamma1[j,]) %*% t(rep(1,s)) %*% t(La)+ La %*% rep(1,s) %*% t(eta[j] - La %*% gamma1[j,]) + 2*a_nu *tau[j] * La  %*% rep(1,s) %*% t(rep(1,s)) %*% t(La)
    nu2 = sum(diag(nuu2))
    
    nuj = 1/2 * nu1 - 1/2 * da_nu * nu2
    snu = rbind(snu,nuj)
  }
   #SEmu = sqrt(diag(solve(t(smu)%*%smu)))
   #SEb  = sqrt(diag(solve(t(sb)%*%sb)))
   SEla  = sqrt(diag(solve(t(sla)%*%sla)))
   ss = cbind(smu, sb, sdr, snu)  
   SE.tmp = sqrt(diag(solve(t(ss)%*%ss))) 
  
  para = c(la, mu, B[vech], diag(D), nu) 
  para.name = c(paste("la ",1:q,sep=""))  
  para.name = c(para.name,paste("mu ",1:p,sep=""))
  para.name = c(para.name,paste("B ",outer(1:p,1:q,paste,sep=",")[vech],sep=""))
  para.name = c(para.name,paste("dd ",1:p,sep=""))
  para.name = c(para.name,paste("nu ",1,sep=""))
  #SE = cbind(est=para, se1=c(SEmu, SEb, SE.tmp))
  SE = cbind(est=para, se1=c(SEla, SE.tmp))
  rownames(SE)=para.name
  list(SE=SE)
}
