[![R-CMD-check](https://github.com/CodesByChris/potentiality/actions/workflows/check-release.yaml/badge.svg)](https://github.com/CodesByChris/potentiality/actions/workflows/check-release.yaml)
[![Lint](https://github.com/CodesByChris/potentiality/actions/workflows/lint.yaml/badge.svg)](https://github.com/CodesByChris/potentiality/actions/workflows/lint.yaml)
[![GitHub](https://img.shields.io/github/license/CodesByChris/potentiality?label=License)](LICENSE)




# Potentiality

R library (`potentiality`) to compute the entropy and potentiality of [ghypernet](https://ghyper.net/) network ensembles.

`potentiality` implements the potentiality measure introduced in the research paper
*"What Is the Entropy of a Social Organization?"* by Christian Zingg, Giona Casiraghi, Giacomo Vaccario, and Frank Schweitzer, published in the journal *Entropy (2019, 21(9), 901; doi:10.3390/e21090901)*.


## Installation

The package can be installed directly from GitHub with the following R commands:

```R
library(devtools)

# Install latest release version
install_github("CodesByChris/potentiality", ref = github_release())

# Or: Install latest development version
# install_github("CodesByChris/potentiality")
```


## Usage

`potentiality` provides two functions:
- `potentiality(...)`: computes a social organization's potentiality from a given (i) interaction network or (ii) `ghype` ensemble.
    Option (i) computes the potentiality as explained in Section 3 of *Zingg et al. (2019)* and internally fits the `ghype` ensemble for the potentiality computation.
    However, this approach may overfit networks with very many nodes.
    For such cases, option (ii) gives the users more control over the fitting process of the `ghype` ensemble, for example, enabling them to fit groups of nodes instead of individual nodes (see also [bccm](https://ghyper.net/reference/bccm.html)).
    `potentiality(...)` is computed according to equation (10) in *Zingg et al. (2019)*.
- `entropy_ghype(...)`: computes an approximate Shannon entropy of a [ghype ensemble](https://github.com/gi0na/r-ghypernet).
    The approximate Shannon entropy is computed according to equations (7) and (8) in *Zingg et al. (2019)*.
    `entropy_ghype(...)` corresponds to `potentiality(...)` with the exception that it is not normalized to the interval $[0, 1]$.


## Examples

To compute the potentiality of the Karate Club, you can proceed as follows:

```R
library(potentiality)
library(magrittr)
library(igraph)

data(karate, package = "igraphdata")

# Convert karate from weighted to multiedge
karate_multi <- karate %>%
    get.adjacency(attr = "weight", sparse = FALSE) %>%
    graph_from_adjacency_matrix(mode = "undirected")

# Compute potentiality
potentiality(karate_multi)  # 0.31, see Table 1 in Zingg et al. (2019)
```


## Contributors

The code in this repository has been developed by

- Christian Zingg
- Giona Casiraghi

at the Chair of Systems Design, ETH Zurich.


## Copyright

The `potentiality` R library is released under the GNU Affero General Public License v3.0

Copyright 2019-2023, ETH Zurich.
