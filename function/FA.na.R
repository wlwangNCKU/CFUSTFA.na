#FA.na

#FA.na-EM
FA.na.EM = function(Y.na, q, init.para, eta=0.005, tol, max.iter, per)
{
  begin = proc.time()[1]
  n = nrow(Y.na)
  p = ncol(Y.na)
  iter = 0
  na.posi = is.na(Y.na)
  mu = init.para$mu  
  B  = init.para$B  
  D  = init.para$D  
  nu=init.para$nu
  s2 = diag(D)
  na.class = colSums(t(na.posi) * 2 ^ (0:(p-1)))
  uni.na.class = unique(na.class)
  num.na.class = length(uni.na.class)
  Ip = diag(p)
  O.list = ind.list = as.list(num.na.class)
  tmp.na.class = rep(NA, n)
  for(i in 1:num.na.class)
  {
    ind = which(na.class == uni.na.class[i])
    ind.list[[i]] = ind
    tmp.na.class[ind] = i
    O.list[[i]] = matrix(Ip[!na.posi[ind[1],],], ncol = p)
  }
  na.class = tmp.na.class
  size.na.class = as.vector(table(na.class))
  det.o = Delta.o = rep(NA, n)
  uni.det.o = rep(NA, num.na.class)
  uni.Soo = array(NA, dim = c(num.na.class, p, p))
  Y.h = t(Y.na)
  Y.h[t(na.posi)] = 99999
  poj = rowSums(!na.posi)
  po = sum(poj)
  S = B %*% t(B) + diag(s2)
  y.cent = Y.h - mu
  for(i in 1:num.na.class)
  {
    O = O.list[[i]]
    OSO = O %*% S %*% t(O)
    uni.det.o[i] = det(OSO)
    uni.Soo[i,,] = t(O) %*% solve(OSO) %*% O
    ind = ind.list[[i]]
    ind.e = matrix(y.cent[,ind], nrow = p)
    Delta.o[ind] = colSums(ind.e * (uni.Soo[i,,] %*% ind.e))
  }
  det.o = uni.det.o[na.class]
  logli.old = - (po*log(2*pi) + sum(log(det.o)) + sum(Delta.o)) / 2
  cat(paste("FA.na.EM is running...", sep = ""), "\n")
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  cat("iter = ", iter, ",\t logli = ", logli.old, sep = "", "\n")
  Soo.ycent = matrix(NA, p, n)
  repeat{
    iter = iter + 1
    for(i in 1:num.na.class)
    {
      ind = ind.list[[i]]
      ind.e = matrix(y.cent[,ind], nrow = p)
      Soo.ycent[,ind] = uni.Soo[i,,] %*% ind.e
    }
    Y.h = mu + S %*% Soo.ycent
    sum.Soo = matrix(colSums(size.na.class *t(apply(uni.Soo, 1, as.vector))), p)
    S.inv = solve(S)
    ga = S.inv %*% B
    om = diag(q) - t(ga) %*% B
    u.j=t(ga)%*%(Y.h - mu)
    sum.om=n * S - S %*% sum.Soo %*% S
    #M-step
    mu=rowMeans((Y.h-(B%*%u.j )))
    B=( (Y.h - mu)%*%t(u.j)+ sum.om%*%ga)%*% solve( u.j%*%t(u.j)+ t(ga)%*%sum.om%*%ga+  (n*om) )
    s2=diag(  (  ((Y.h - mu)-(B%*%u.j))%*%t((Y.h - mu)-(B%*%u.j))  )+ sum.om -  (sum.om%*%ga%*%t(B))-(B%*%t(ga)%*%sum.om)+ (B%*%t(ga)%*%sum.om%*%ga%*%t(B)) +  ((B%*%om%*%t(B))*n)     )/n
    s2[s2<eta] = eta
    S = B %*% t(B) + diag(s2)
    y.cent = Y.h - mu
    for(i in 1:num.na.class)
    {
      O = O.list[[i]]
      OSO = O %*% S %*% t(O)
      uni.det.o[i] = det(OSO)
      uni.Soo[i,,] = t(O) %*% solve(OSO) %*% O
      ind = ind.list[[i]]
      ind.e = matrix(y.cent[,ind], nrow = p)
      Delta.o[ind] = colSums(ind.e * (uni.Soo[i,,] %*% ind.e))
    }
    det.o = uni.det.o[na.class]
    logli.new = - (po*log(2*pi) + sum(log(det.o)) + sum(Delta.o)) / 2
    diff = logli.new - logli.old
    if (iter%%per == 0) cat("iter =", iter, "\tlogli =", logli.new, "\tdiff =", diff, "\n")
    if (diff < tol | iter == max.iter) break
    logli.old = logli.new
  }
  m =  p + ( p*q-q*(q-1)/2 ) + p
  AIC = -2 * logli.new + m * 2
  BIC = -2 * logli.new + m * log(n) 
  mds = c(no.para = m, logli = logli.new, AIC = AIC, BIC = BIC)
  na.ind = which(as.vector(t(na.posi))==T)
  Y.mis = as.vector(Y.h)[na.ind]
  para = list(mu = mu,  B = B,  s2 = diag(s2), nu = nu, S = S)
  end = proc.time()[1]
  cat(paste(rep("=",60),sep="",collapse=""),"\n")
  cat("iter =", iter, "\tlogli =", logli.new, "\tdiff =", diff, "\n")
  cat('FA.na.EM takes', end - begin, 'seconds.\n')
  list(mds=mds, iter=iter, para=para, Y.hat=t(Y.h), U=u.j, Y.mis=Y.mis, cpu.time=end - begin)
}

