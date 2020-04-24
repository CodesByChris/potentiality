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

    # Get xi and omega
    if (!is.null(ensemble)){
        if (is.null(xi))
            xi <- ensemble$xi
        if (is.null(omega))
            omega <- ensemble$omega
        directed <- ensemble$directed
        selfloops <- ensemble$selfloops
    }
    ix <- mat2vec.ix(xi, directed, selfloops)
    xi <- xi[ix]
    omega <- omega[ix]

    # Get m
    if (is.null(m))
        m <- ensemble$m

    # Compute the ps (according to Eq.(8) in DOI:10.3390/e21090901)
    ps <- xi*omega/sum(xi*omega)
    ps <- ps[ps != 0]
    # Note: For empirical systems with many disconnected components in the
    #     interaction network this filtering tremendously speeds up the
    #     computation, given that most entries of ps are zero due to the
    #     disconnected components.


    # Special Case: 1 Probability
    if (length(ps) == 1) {
        if (ps != 1)
            stop("Single probability != 1 encountered: ", ps)
        return(0)
    }

    # Compute last term (using a trick from scipy with a nested binomial PMF, see https://github.com/scipy/scipy/blob/v1.4.1/scipy/stats/_multivariate.py#L3158)
    logfactorial_table <- .logfactorialtable(m)
    sapply(2:m,
           function(x, ps, logfactorial_table) {
               logfactorial_table[x] * sum(
                   Vectorize(dbinom, "prob")(x = x, size = m, prob = ps)
               )
           },
           ps = ps, logfactorial_table = logfactorial_table) %>%
        sum() ->
        last_term

    # Compute the entropy (Eq.(7) in DOI:10.3390/e21090901)
    Hval <- -logfactorial_table[m] - m*sum(ps*log(ps)) + last_term
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

    # Compute the ps (according to Eq.(9), i.e. all node-pairs equally likely)
    k <- sum(ix)
    ps <- 1/k

    # Compute the entropy (using a numpy-trick as in .H.ghype above)
    x <- 2:m
    logfactorial_table <- .logfactorialtable(m)
    last_term <- k * sum(logfactorial_table[x] * Vectorize(dbinom, "x")(x = x, size = m, prob = ps))
    Hval <- -logfactorial_table[m] - m*log(ps) + last_term
    return('H' = Hval)
}


.entropyRatio <- function(model){
    observed_H <- .H.ghype(model)
    if (observed_H == 0)
        return(0)
    return(observed_H / .maxEntropy(model))
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
