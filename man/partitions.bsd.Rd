\name{partitions.bsd}
\alias{partitions.bsd}
\title{
Obtain optimal partitioning scheme
}
\description{
This function obtains the optimal number of partitions for gene trees. The input should be an object of class 'bsd' with the pairwise distances among tres and the scaling factors (see bsd.matrix for more details on this class).
partitions.bsd returns an object of class partitions with a matrix with the parition assignments, the range of the sBSDmin distances and an object of class cluster with the output from the clustering analysis.
}
\usage{
partitions.bsd(bsd.object, FUN = pam, find.best = T, B = 500, gap.best = "firstSEmax", kmax = "", ...)
}
\arguments{
  \item{bsd.object}{
An object of class bsd. It can be obtained witn bsd.matrix or bsd.matrix.para. The example(bsd.matrix) for details.
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
\examples{


# Create a list of trees of class multiPhylo with four patterns of among-lineage rate variation 
tr.fix <- rtree(10)

set.seed(12345)
rates1 <- abs(rnorm(18, sd = 0.1))
set.seed(123456)
rates2 <- abs(rnorm(18, sd = 0.1))
set.seed(1234567)
rates3 <- abs(rnorm(18, sd = 0.1))
set.seed(12345678)
rates4 <- abs(rnorm(18, sd = 0.1))

trees.list <- list()

for(i in 1:20){
      trees.list[[i]] <- tr.fix
      if(i <= 5){
      	   trees.list[[i]]$edge.length <- abs(rates1 + rnorm(18, 0, 0.01)) 
      }else if(i > 5 && i <= 10){
      	    trees.list[[i]]$edge.length <- abs(rates2 + rnorm(18, 0, 0.01)) 
      }else if(i >= 10 && i < 15){
      	    trees.list[[i]]$edge.length <- abs(rates3 + rnorm(18, 0, 0.01))  
      }else{
	    trees.list[[i]]$edge.length <- abs(rates4 + rnorm(18, 0, 0.01)) 
      }
}

names(trees.list) <- paste0("tree", 1:20)
class(trees.list) <- "multiPhylo"

# Estimate sBSDmin distance for all pairs of trees:
trees.bsd <- bsd.matrix(trees.list)

# Find the optimal number of partitions:
partitions.object <- partitions(trees.bsd)

# plot partitions in two graphics devices
plot(partitions.object)





## The function is currently defined as
function (bsd.object, FUN = pam, find.best = T, B = 500, gap.best = "firstSEmax", 
    kmax = "", ...) 
{
    dimat <- as.matrix(bsd.object[[1]])
    if (kmax == "") {
        kmax <- nrow(dimat) - 1
    }
    dimat.mds <- cmdscale(dimat, k = 3)
    if (find.best == T) {
        clusdat <- clusGap(dimat.mds, B = B, FUNcluster = FUN, K.max = kmax)
        npart <- maxSE(f = clusdat$Tab[, 3], SE.f = clusdat$Tab[, 
            4], method = gap.best)
    }
    parts.list <- list()
    for (i in 1:kmax) {
        clus.temp <- FUN(dimat.mds, k = i, ...)
        parts.list[[i]] <- clus.temp$clustering
    }
    parts.mat <- do.call("cbind", parts.list)
    rownames(parts.mat) <- rownames(dimat.mds)
    colnames(parts.mat) <- paste0("k=", 1:ncol(parts.mat))
    res <- list(parts.mat, range(bsd.object[[1]]))
    names(res) <- c("parts.mat", "range.bsd")
    if (exists("npart")) {
        colnames(res[[1]])[npart] <- paste0(colnames(parts.mat)[npart], 
            "BEST")
        res[[3]] <- npart
        names(res)[3] <- "best.k"
    }
    if (find.best == T) {
        res[[4]] <- clusdat
        names(res)[4] <- "clusterdata"
    }
    class(res) <- "partitions"
    return(res)
  }
}
\keyword{ clustering }
\keyword{ gap statistic }
