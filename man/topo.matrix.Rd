\name{topo.matrix}
\alias{topo.matrix}
\title{
Estimate the topology distance for all pairs of trees in a list
}
\description{
topo.matrix estimates the topology distance for all paris of trees in a list. It returns an object of class topo with the pairswise distances and the scaling factors. The scaling factors represent the differences in the magnitude of the rates among trees. See the help file for dist.topo for more details.
}
\usage{
topo.matrix(tree.list)
}
\arguments{
  \item{tree.list}{
A list of trees. It can be a list where each item is a tree, or an object of class multiPhylo. The trees can have names. See the example bellow.
}
}

\value{
An object of class 'topomat', which is a list with 
\item{topo.mat }{An object of class topomat. This is a list with a distance object for the pairwise tree distances}
}

\author{
Sebastian Duchene
}

\seealso{
dist.topo from package ape
}

\keyword{ tree topology distance }
