logfactorial <- Vectorize(function(n){
  if(n<25) return(log(factorial(n)))
  return(1/2*log(2*pi) + (n+1/2)*log(n) - n)
})

H.ghype <- function(ensemble = NULL, xi = NULL, omega = NULL, directed = NULL, selfloops = NULL, m = NULL){
  # entropy of the multinomial approximation

  try(if(is.null(ensemble) & (is.null(xi) | is.null(omega) | is.null(directed) | is.null(selfloops)) )
    stop('specify ensemble'))

  if(!is.null(ensemble)){
    if(is.null(xi))    xi <- ensemble$xi
    if(is.null(omega))    omega <- ensemble$omega
    directed <- ensemble$directed
    selfloops <- ensemble$selfloops
  }
  ix <- mat2vec.ix(xi, directed, selfloops)
  xi <- xi[ix]
  omega <- omega[ix]

  pp <- sum(xi*omega)
  ps <- xi*omega/pp
  nnzeros <- ps!=0
  ps <- ps[nnzeros]

  pstar <- max(ps)
  k <- sum(ix)
  if(is.null(m))  m <- ensemble$m
  x <- seq(10, 140, by=10)
  xmax <- x[which(log10(k) + log10(choose(m, x)) + log10(x) + log10(log(x)) + x*log10(pstar) < -2)[1]]
  if (is.na(xmax))
      xmax <- 140

  x <- t(matrix(2:xmax, xmax-1, sum(nnzeros)))

  Hval <- - logfactorial(m) - m*sum(ps*log(ps)) + sum( choose(m, x) * ps^x * (1-ps)^(m-x) * matrix(logfactorial(x), nrow=sum(nnzeros)) )
  return('H' = Hval)
}

maxEntropy <- function(ensemble=NULL, m=NULL, N=NULL, directed=NULL, selfloops=NULL){
  if(is.null(ensemble)){
    mat <- matrix(0, N, N)
  } else{
    mat <- ensemble$xi
    m <- ensemble$m
    N <- nrow(mat)
    directed <- ensemble$directed
    selfloops <- ensemble$selfloops
  }

  ix <- mat2vec.ix(mat, directed, selfloops)

  k <- sum(ix)

  ps <- 1/k

  x <- 2:20
  sums <-  k * choose(m,x) * ps^x * (1-ps)^(m-x) * logfactorial(x)
  while(sums[length(sums)]>1e-2){
    x <- (max(x)+1):(max(x)+10)
    sums <- c(sums, k*choose(m,x)*ps^x * (1-ps)^(m-x) * logfactorial(x))
  }

  Hval <- - logfactorial(m) - m*log(ps) + sum(sums)
  return('H' = Hval)
}

entropyRatio <- function(model){
  return(H.ghype(model)/maxEntropy(model))
}

# H.ghype(model0)
# H.ghype(ensemble = model0, xi=matrix(1,10,10))

# adj <- vec2mat(rmultinom(1,180,rep(1/499500,499500)), F, F, 1000)
# adj <- adj + t(adj)
# randmodel <- scm(adj, F, F)
#
# maxEntropy(randmodel)
# H.ghype(randmodel)
#
# entropyRatio(randmodel)

