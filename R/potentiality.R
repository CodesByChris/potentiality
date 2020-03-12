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

.logfactorialtable <- Vectorize(function(n){
    if (n == 0)
        return(0)  # log(0!) := log(1) = 0
    return(cumsum(sapply(1:n, log)))
})


.logapproxchoose <- Vectorize(function(n, k, logfactorial_table) {
    stopifnot(0 <= k && k <= n)
    if (k == n || k == 0)
        return(0)
    return(logfactorial_table[n] - logfactorial_table[k] - logfactorial_table[n-k])
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
    logfactorial_table <- .logfactorialtable(m)
    sums <-  sum( exp(matrix(.logapproxchoose(m, x, logfactorial_table), sum(nnzeros)) + x*log(ps) + log(1-ps)*(m-x) + log(matrix(logfactorial_table[x], nrow=sum(nnzeros)))) )

    Hval <- - logfactorial_table[m] - m*sum(ps*log(ps)) + sum(sums)
    return('H' = Hval)
}


.maxEntropy <- function(ensemble=NULL, m=NULL, N=NULL, directed=NULL, selfloops=NULL){
    if (is.null(ensemble)) {
        mat <- matrix(0, N, N)
    } else {
        mat <- ensemble$xi
        m <- ensemble$m
        N <- nrow(mat)
        directed <- ensemble$directed
        selfloops <- ensemble$selfloops
    }

    ix <- mat2vec.ix(mat, directed, selfloops)

    k <- sum(ix)

    ps <- 1/k

    x <- 2:m
    logfactorial_table <- .logfactorialtable(m)
    sums <-  k * exp(.logapproxchoose(m, x, logfactorial_table) + x * log(ps) + (m-x) * log(1-ps) + log(logfactorial_table[x]))

    Hval <- - logfactorial_table[m] - m*log(ps) + sum(sums)
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
