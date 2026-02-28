#generate missing values
gener.na = function(X, na.rate)
{
  n = nrow(X)
  p = ncol(X)
  num.na = floor(n*p*na.rate)
  keep.part = sample(p, n, replace = T) + p * 0:(n-1)
  na.posi = tabulate(sample((1:(n*p))[-keep.part], num.na), nbins = n*p)
  na.posi = matrix(na.posi, ncol = p, byrow = T)
  X[na.posi == 1] = NA
  return(X)
}