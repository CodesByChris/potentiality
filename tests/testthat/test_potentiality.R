# Tests for potentiality.

library(magrittr)


################################################################################
### Test: Full Network
################################################################################

test_that("Undirected full network has potentiality 1", {
  n_nodes <- 30
  multiplicity <- 5
  ntwk <- igraph::make_empty_graph(n_nodes, directed = FALSE)
  for (i in 1:igraph::vcount(ntwk)) {
    for (j in 1:igraph::vcount(ntwk)) {
      if (i == j)
        next
      ntwk <- igraph::add_edges(ntwk, rep(c(i, j), times = multiplicity))
    }
  }
  expect_equal(potentiality(ntwk), 1)
})


test_that("Directed full network has potentiality 1", {
  n_nodes <- 30
  multiplicity <- 5
  ntwk <- igraph::make_empty_graph(n_nodes, directed = TRUE)
  for (i in 1:igraph::vcount(ntwk)) {
    for (j in 1:igraph::vcount(ntwk)) {
      if (i == j)
        next
      ntwk <- igraph::add_edges(ntwk, rep(c(i, j), times = multiplicity))
    }
  }
  expect_equal(potentiality(ntwk), 1)
})


################################################################################
### Test: Empty Network
################################################################################

test_that("Undirected empty network has potentiality 0", {
  n_nodes <- 30
  ntwk <- igraph::make_empty_graph(n_nodes, directed = FALSE)
  expect_equal(potentiality(ntwk), 0)
})


test_that("Directed empty network has potentiality 0", {
  n_nodes <- 30
  ntwk <- igraph::make_empty_graph(n_nodes, directed = TRUE)
  expect_equal(potentiality(ntwk), 0)
})


################################################################################
### Test: Single pair
################################################################################

test_that("potentiality 0 when only one pair communicates (undirected)", {
  n_nodes <- 30
  multiplicity <- 10000
  ntwk <- igraph::make_empty_graph(n_nodes, directed = FALSE)
  ntwk <- igraph::add_edges(ntwk, rep(c(1, 2), times = multiplicity))
  expect_equal(potentiality(ntwk), 0)
})


test_that("potentiality 0 when only one pair communicates (directed)", {
  n_nodes <- 30
  multiplicity <- 10000
  ntwk <- igraph::make_empty_graph(n_nodes, directed = TRUE)
  ntwk <- igraph::add_edges(ntwk, rep(c(1, 2), times = multiplicity))
  expect_equal(potentiality(ntwk), 0)
})


################################################################################
### Test: Multinomial entropy computation
################################################################################

test_that("Directed network with 2 nodes and unidirectional preference", {
  weak_multiplicity <- 6
  strong_multiplicity <- 10000
  ntwk <- igraph::make_empty_graph(2, directed = TRUE)
  ntwk <- igraph::add_edges(ntwk, rep(c(1, 2), times = weak_multiplicity))
  ntwk <- igraph::add_edges(ntwk, rep(c(2, 1), times = strong_multiplicity))
  expect_equal(potentiality(ntwk), 0.431235, tolerance = 1E-5)
  # 0.431235 is computed by applying SciPy's multinomial entropy to the p_ij
  # from Eq.(8) in DOI:10.3390/e21090901.
})


################################################################################
### Test: Undirected star example from DOI:10.3390/e21090901
################################################################################

test_that("Undirected star", {
  i_center <- 1
  i_others <- (1:10)[-i_center]
  ntwk <- igraph::make_empty_graph(10, directed = FALSE)
  for (j in i_others) {
    ntwk <- igraph::add_edges(ntwk, rep(c(i_center, j), times = 10))
  }
  pot <- potentiality(ntwk)
  expect_true(0.265 <= pot && pot <= 0.275)
})


################################################################################
### Test: Karate club example from DOI:10.3390/e21090901
################################################################################

test_that("Karate club", {
  data(karate, package = "igraphdata")
  pot <- karate %>%
    igraph::get.adjacency(attr = "weight", sparse = FALSE) %>%
    igraph::graph_from_adjacency_matrix(mode = "undirected") %>%
    potentiality()
  expect_true(0.305 <= pot && pot <= 0.315)
})
