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
    sapply(2:m,
           function(x, ps) {
               lfactorial(x) * sum(
                   Vectorize(dbinom, "prob")(x = x, size = m, prob = ps)
               )
           },
           ps = ps) %>%
        sum() ->
        last_term

    # Compute the entropy (Eq.(7) in DOI:10.3390/e21090901)
    return(-lfactorial(m) - m*sum(ps*log(ps)) + last_term)
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
    last_term <- k * sum(lfactorial(x) * Vectorize(dbinom, "x")(x = x, size = m, prob = ps))
    return(-lfactorial(m) - m*log(ps) + last_term)
}


.entropyRatio <- function(model){
    observed_H <- .H.ghype(model)
    if (observed_H == 0)
        return(0)
    return(observed_H / .maxEntropy(model))
}


Potentiality <- function(network,
                         directed = is_directed(network),
                         has_selfloops = any(which_loop(network)),
                         full_model = TRUE) {
    # Computes the potentiality according to a MLE fit.
    #
    # The MLE fit has the property that network is preserved as the expected
    # network over the ensemble.
    #
    # For the special cases where network has no nodes or no edges a
    # potentiality of 0 is computed, because there is always only this *one*
    # network which fulfils these constraints.
    #
    # Args: network: igraph graph of which to compute the potentiality directed:
    # (optional) Whether network is directed. If omitted, this is detected from
    # base_network. has_selfloops: (optional) Whether base_network is allowed to
    # have self-loops. If omitted, this is detected from base_network.
    # full_model: (optional) whether to use full ghype or bccm. If omitted it
    # defaults to TRUE
    #
    # Returns: The computed potentiality.

    # Handle empty networks
    if (vcount(network) == 0 || ecount(network) == 0)
        return(0)

    # Compute params
    if(has_selfloops)
        stop("ERROR: Currently potentiality of selfloop-network not implemented.")

    # Compute Potentiality
    ## if full_model -> use full ghype propensities
    if(isTRUE(full_model)){
        suppressWarnings(
          ghype(network, directed = directed, selfloops = has_selfloops, unbiased = FALSE) %>%
            .entropyRatio() -> pot_val
        )
    }

    ## if !full_model -> use bccm with inferred blocks
    if(isFALSE(full_model)){
        net <- igraph::graph_from_adjacency_matrix(get.adjacency(network), weighted = TRUE)
        labs <- igraph::membership(igraph::cluster_fast_greedy(graph = igraph::as.undirected(net), modularity = FALSE))
        suppressWarnings(
          bccm(igraph::get.adjacency(network, sparse = F), labels = labs, directed = directed, selfloops = has_selfloops, ignore_pvals = TRUE) %>%
            .entropyRatio() -> pot_val
        )
    }
  return(pot_val)
}
