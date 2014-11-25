\name{partitions.topomat}
\alias{partitions.topomat}
\title{
Obtain optimal topology partitioning scheme
}
\description{
This function obtains the optimal number of partitions for gene trees with different topologies. The input should be an object of class 'topomat' with the pairwise distances among trees  (see topo.matrix for more details on this class).
partitions.topomat returns an object of class partitions with a matrix with the parition assignments, the range of the sBSDmin distances and an object of class cluster with the output from the clustering analysis.
}
\usage{
partitions.topomat(topomat.object, FUN = pam, find.best = T, B = 500, gap.best = "firstSEmax", kmax = "", ...)
}
\arguments{
  \item{topomat.object}{
An object of class topomat. It can be obtained witn topo.matrix or topo.matrix.para (to be implemented).
}
  \item{FUN}{
A clustering function as defined in the package cluster. It can be pam, clara, or fanny. However note that for fanny the maximum k should be the number of data subsets / 2.
}
  \item{find.best}{
logical. Select T to find the optimal number of partitions. If F, the function does not find the best number of partitions.
}
  \item{B}{
Numeric. Number of bootstrap replicates for the Gap statistic. Not used if find.best == F
}
  \item{gap.best}{
A character. The criterion to select the optimal number of partitions (k). It can be any of "firstSEmax", "Tibs2001SEmax", "globalSEmax", "firstmax", or "globalmax".
}
  \item{kmax}{
Numeric. Maximum number of partitions to test. The default is the number of data subsets - 1.
}
  \item{\dots}{
Additional arguments passed to the clustering function.
}
}
\details{
see the help for pacakge clustering for more details.
}
\value{
An object of class 'partitions'
\item{parts.mat }{A matrix with the cluster assignments. The partitions are numbered from 1:maximum k. The columns are each value of k, and the columns are the names of the data subsets if provided.}
\item{range.bsd }{A numeric vector with the range of the sBSDmin distances}
\item{best.k }{The optimal number of partitions. Only returned if find.best == T}
\item{clusterdata }{Results from the Gap statistic. Only returned if find.best == T. See ?clusGap for more details}
}
\references{

     Tibshirani, R., Walther, G. and Hastie, T. (2001).  Estimating the
     number of data clusters via the Gap statistic.  _Journal of the
     Royal Statistical Society B_, *63*, 411-423.

     Tibshirani, R., Walther, G. and Hastie, T. (2000).  Estimating the
     number of clusters in a dataset via the Gap statistic.  Technical
     Report. Stanford.

     Per Broberg (2006). SAGx: Statistical Analysis of the GeneChip.  R
     package version 1.9.7.  <URL:
     http://home.swipnet.se/pibroberg/expression_hemsida1.html>

}
\author{
Sebastian Ducnene
}
\keyword{ clustering }
\keyword{ gap statistic }
