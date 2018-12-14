logfactorial <- Vectorize(function(n){
  if(n<170) return(log(factorial(n)))
  return(1/2*log(2*pi) + (n+1/2)*log(n) - n)
})


logapproxchoose <- Vectorize(function(n, k) {
    if(n <= 6000 && k <= 140)  # chosen by hand
        return(log(choose(n, k)))

    # Stirling approximation of binomial coefficient
    return( (0.5 + n)*log(n) - (0.5 + k)*log(k) - (0.5 + n - k)*log(n - k) - 0.5 * log(2 * pi) )
}, vectorize.args = "k")


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

  # pstar <- max(ps)
  # k <- sum(ix)
  if(is.null(m))  m <- ensemble$m

  xmax <- m %/% 2
  x <- t(matrix(2:xmax, xmax-1, sum(nnzeros)))
  sums <-  sum( exp(matrix(logapproxchoose(m, x), sum(nnzeros)) + x*log(ps) + log(1-ps)*(m-x) + log(matrix(logfactorial(x), nrow=sum(nnzeros)))) )
  # while(sums[length(sums)]>1e-200){
  #     x <- t(matrix((max(x)+1):(max(x)+1), 1, sum(nnzeros)))
  #     sums <- c(sums, sum( exp(matrix(logapproxchoose(m, x), sum(nnzeros)) + x*log(ps) + log(1-ps)*(m-x) + log(matrix(logfactorial(x), nrow=sum(nnzeros)))) ))
  # }

  Hval <- - logfactorial(m) - m*sum(ps*log(ps)) + sum(sums)
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

  x <- 2:(m %/% 2)
  sums <-  k * exp(logapproxchoose(m,x) + x * log(ps) + (m-x) * log(1-ps) + log(logfactorial(x)))
  # while(sums[length(sums)]>1e-2){
  #   x <- (max(x)+1):(max(x)+10)
  #   sums <- c(sums, k*approxchoose(m,x)*ps^x * (1-ps)^(m-x) * logfactorial(x))
  # }

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

