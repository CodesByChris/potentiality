# Potentiality implementation.


#' Compute vector v, where v_i := sum_{p in probs} dbinom(xs_i, size, p)
.dbinom_sum_over_probs <- function(xs, size, probs) {

    # Get unique probabilities and their frequencies
    probs_counts <- plyr::count(tibble::tibble(.probs = probs), vars = ".probs")

    # Get vector of binom sums over p_ij (entries correspond to xs)
    binom_sums_per_x <- numeric(length = length(xs))
    for (i in seq_along(probs_counts$.probs)) {
        freq <- probs_counts$freq[[i]]
        prob <- probs_counts$.probs[[i]]
        binoms <- stats::dbinom(x = xs, size = size, prob = prob)
        binom_sums_per_x <- binom_sums_per_x + freq * binoms
    }
    return(binom_sums_per_x)
}


#' Entropy of multinomial approximation.
.entropy_ghype <- function(ensemble = NULL, xi = NULL, omega = NULL,
                           directed = NULL, selfloops = NULL, m = NULL) {

    # Validate args
    if (is.null(ensemble) && (is.null(xi) || is.null(omega) ||
                              is.null(directed) || is.null(selfloops)))
        stop("Either `ensemble` or separate params must be given.")

    if (!is.null(ensemble)) {
        xi <- ensemble$xi
        omega <- ensemble$omega
        directed <- ensemble$directed
        selfloops <- ensemble$selfloops
    }

    if (is.null(m))
        m <- ensemble$m

    # Get xi and omega
    ix <- ghypernet::mat2vec.ix(xi, directed, selfloops)
    xi <- xi[ix]
    omega <- omega[ix]

    # Compute the ps (according to Eq.(8) in DOI:10.3390/e21090901)
    ps <- xi * omega / sum(xi * omega)
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

    # Compute last term (using a trick from scipy with a nested binomial PMF,
    # see github.com/scipy/scipy/blob/v1.4.1/scipy/stats/_multivariate.py#L3158)
    xs <- 2:m
    last_term <- sum(lfactorial(xs) * .dbinom_sum_over_probs(xs, m, ps))

    # Compute the entropy (Eq.(7) in DOI:10.3390/e21090901)
    return(-lfactorial(m) - m * sum(ps * log(ps)) + last_term)
}


#' Maximum entropy of multinomial approximation.
.max_entropy <- function(ensemble = NULL, m = NULL, n = NULL,
                         directed = NULL, selfloops = NULL) {
    if (is.null(ensemble)) {
        mat <- matrix(0, n, n)
    } else {
        mat <- ensemble$xi
        m <- ensemble$m
        n <- nrow(mat)  # TODO: Is this even used?
        directed <- ensemble$directed
        selfloops <- ensemble$selfloops
    }
    ix <- ghypernet::mat2vec.ix(mat, directed, selfloops)

    # Compute the ps (according to Eq.(9), i.e. all node-pairs equally likely)
    k <- sum(ix)
    p <- 1 / k

    # Compute the entropy (using a numpy-trick as in .entropy_ghype above)
    xs <- 2:m
    binoms <- stats::dbinom(x = xs, size = m, prob = p)
    last_term <- k * sum(lfactorial(xs) * binoms)
    return(-lfactorial(m) - m * log(p) + last_term)
}


#' Computes the model's entropy normalized by the theoretical maximum.
.relative_entropy <- function(model) {
    entropy <- .entropy_ghype(model)
    if (entropy == 0)
        return(0)
    return(entropy / .max_entropy(model))
}


#' Computes the potentiality according to a MLE fit.
#'
#' The MLE fit has the property that network is preserved as the expected
#' network over the ensemble.
#'
#' For the special cases where network has no nodes or no edges, a potentiality
#' of 0 is returned because there is always only this *one* network which
#' fulfils these constraints.
#'
#' @param network igraph graph of which to compute the potentiality.
#' @param directed Whether network is directed. If omitted, this is detected
#'   from base_network.
#' @param has_selfloops Whether base_network is allowed to have self-loops. If
#'   omitted, this is detected from base_network.
#' @param full_model Whether to use a full ghype (TRUE) or bccm (FALSE).
#'   Defaults to TRUE.
#' @return The computed potentiality.
#' @export
potentiality <- function(network,
                         directed = igraph::is_directed(network),
                         has_selfloops = any(igraph::which_loop(network)),
                         full_model = TRUE) {

    # Handle empty networks
    if (igraph::vcount(network) == 0 || igraph::ecount(network) == 0)
        return(0)

    # Compute params
    if (has_selfloops)
        stop("ERROR: Potentiality not implemented for selfloops.")

    # Compute Potentiality
    if (isTRUE(full_model)) {
        # Use full ghype propensities
        ens <- ghypernet::ghype(network, directed = directed,
                                selfloops = has_selfloops, unbiased = FALSE)
    } else {
        # Use bccm with inferred blocks
        adj <- igraph::get.adjacency(network, sparse = FALSE)
        net <- igraph::graph_from_adjacency_matrix(adj, weighted = TRUE)
        labs <- igraph::membership(igraph::cluster_fast_greedy(
            graph = igraph::as.undirected(net), modularity = FALSE)
        )
        ens <- ghypernet::bccm(adj, labels = labs, directed = directed,
                               selfloops = has_selfloops, ignore_pvals = TRUE)
    }
    return(.relative_entropy(ens))
}
