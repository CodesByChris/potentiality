# Potentiality
#
# Authors:
#     Christian Zingg and Giona Casiraghi
#
# Source:
#     Function computing Potentiality as in the Potentiality-Project
#     (https://github.com/sg-dev/entropy-social-organisation). However, unlike
#     there the current implementation uses functions from the
#     entropy-adaptation project (https://github.com/sg-dev/entropy-adaptation),
#     which compute the multinomial entropy approximation directly in R.

library(magrittr)
library(dplyr)
library(tibble)
library(ghypernet)



################################################################################
### Potentiality
################################################################################

.logfactorial <- Vectorize(function(n){
    if (n == 0)
        return(1)
    return(sum(sapply(1:n, log)))
})


.logapproxchoose <- Vectorize(function(n, k) {
    stopifnot(0 <= k && k <= n)
    if (k == n || k == 0)
        return(0)
    return(.logfactorial(n) - .logfactorial(k) - .logfactorial(n-k))
}, vectorize.args = "k")


.H.ghype <- function(ensemble = NULL, xi = NULL, omega = NULL, directed = NULL, selfloops = NULL, m = NULL){
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

    x <- t(matrix(2:m, m-1, sum(nnzeros)))
    sums <-  sum( exp(matrix(.logapproxchoose(m, x), sum(nnzeros)) + x*log(ps) + log(1-ps)*(m-x) + log(matrix(.logfactorial(x), nrow=sum(nnzeros)))) )
    # while(sums[length(sums)]>1e-200){
    #     x <- t(matrix((max(x)+1):(max(x)+1), 1, sum(nnzeros)))
    #     sums <- c(sums, sum( exp(matrix(.logapproxchoose(m, x), sum(nnzeros)) + x*log(ps) + log(1-ps)*(m-x) + log(matrix(.logfactorial(x), nrow=sum(nnzeros)))) ))
    # }

    Hval <- - .logfactorial(m) - m*sum(ps*log(ps)) + sum(sums)
    return('H' = Hval)
}


.maxEntropy <- function(ensemble=NULL, m=NULL, N=NULL, directed=NULL, selfloops=NULL){
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
    sums <-  k * exp(.logapproxchoose(m,x) + x * log(ps) + (m-x) * log(1-ps) + log(.logfactorial(x)))
    # while(sums[length(sums)]>1e-2){
    #   x <- (max(x)+1):(max(x)+10)
    #   sums <- c(sums, k*approxchoose(m,x)*ps^x * (1-ps)^(m-x) * .logfactorial(x))
    # }

    Hval <- - .logfactorial(m) - m*log(ps) + sum(sums)
    return('H' = Hval)
}


.entropyRatio <- function(model){
    return(.H.ghype(model)/.maxEntropy(model))
}


Potentiality <- function(network,
                         directed = is_directed(network),
                         has_selfloops = any(which_loop(network))) {
    # Computes the potentiality according to a MLE fit.
    #
    # The MLE fit has the property that network is preserved as the expected
    # network over the ensemble.
    #
    # For the special cases where network has no nodes or no edges a
    # potentiality of 0 is computed, because there is always only this *one*
    # network which fulfils these constraints.
    #
    # Args:
    #     network: igraph graph of which to compute the potentiality
    #     directed: (optional) Whether network is directed. If omitted, this
    #         is detected from base_network.
    #     has_selfloops: (optional) Whether base_network is allowed to have
    #         self-loops. If omitted, this is detected from base_network.
    #
    # Returns:
    #     The computed potentiality.

    # Handle empty networks
    if (vcount(network) == 0 || ecount(network) == 0)
        return(0)

    # Compute params
    if(has_selfloops)
        stop("ERROR: Currently potentiality of selfloop-network not implemented.")

    # Compute Potentiality
    ghype(network, directed = directed, selfloops = has_selfloops) %>%
        .entropyRatio() %>%
        return()
}
