\name{net.gen.from.P}
\alias{net.gen.from.P}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generates a network from the given connection probability
}
\description{
  Generates an adjacency matrix from a given probability matrix,
  according independent Bernoulli -- the so-called inhomogeneous
  Erdos-Renyi model. It is used to generate new networks from a given model.
}
\usage{
net.gen.from.P(P, mode = "undirected")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
connection probability between nodes
}
  \item{mode}{
"undirected" (default) if the network is undirected, so the
adjacency matrix will be symmetric with only upper diagonal entries
being generated as independent Bernoulli. Otherwise, the adjacency matrix gives independent Bernoulli everywhere.
}
}
\value{
An adjacency matrix
%% ...
}
\author{
Can M. Le and  Tianxi Li.\cr

Maintainer: Tianxi Li <tianxili@virginia.edu>
}